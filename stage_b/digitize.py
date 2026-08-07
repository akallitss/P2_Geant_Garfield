"""digitize.py — Stage B: ionization clusters -> per-pad charge and time.

Chain per cluster, per electron (SIM_CAMPAIGN_PLAN §2, §5):

    depth from the drift frame
      -> attachment survival        exp(-depth / lambda_a)
      -> drift time                 depth / v_d
      -> longitudinal diffusion     N(0, sigma_L sqrt(depth) / v_d)  on time
      -> transverse diffusion       N(0, sigma_T sqrt(depth))        on x, y
      -> mesh transparency          Bernoulli(eps)
      -> Polya avalanche gain       Gamma(theta+1, G/(theta+1))
      -> pad lookup                 (x, y) -> pad
      -> sum charge / earliest time per pad

DESIGN NOTES, and the reasons behind the non-obvious choices:

**Per electron, not per packet.** MX17 uses a per-cluster packet
approximation because a 30 mm gap and 512x512 pads make per-electron
expensive. Our numbers are different: a MIP in 3 mm of Ar/iso liberates ~30
electrons per event in total. Per-electron is therefore free, and it is more
faithful exactly where it matters — near a pad edge, where the difference
between "the packet is on pad A" and "80 % of its electrons are" is the
entire charge-sharing signal.

**`nPrimary` comes from Stage A, not from edep/W.** SteppingAction already
applied W and carried the sub-W remainder probabilistically. Re-deriving here
would double-count the conversion and throw away the remainder. (This is
MX17's note in `clusters.py`, and it applies verbatim.)

**AmpGas clusters skip the drift and get partial gain.** An electron released
a distance d into the 150 um amplification gap is amplified over only the
remaining (gap - d), so its gain is G^((gap-d)/gap) — the gain spans about
four decades across that single volume. Summing edepAmp and treating it like
edepDrift is explicitly wrong (`OUTPUT_FORMAT.md` §4); this is the correct
handling of the same fact.

**Timestamp per pad is the EARLIEST electron**, not the mean. That matches
what a leading-edge discriminator does. It is not the same as the
threshold-crossing time the VMM reports — that is Stage C's job
(`vmm/`), which needs the charge and the arrival-time distribution this
stage produces, and applies the shaper and threshold itself.

WHAT THIS STAGE DOES NOT DO YET, all of it recorded rather than hidden:
  * **Induced current is not modelled.** Each electron contributes its charge
    at its arrival time, i.e. a delta. The real signal is a fast electron
    spike plus a ~150 ns ion tail, which is P0.17 — cheap for our
    non-resistive stack (static Ramo, Riegler closed form) and the main model
    risk in the timing prediction.
  * **Mesh transparency is a constant**, not a field map (P0.13, Phase 4).
  * **No pillars.** P0.16: they are not in the Geant4 geometry either, so
    electrons that should land on a pillar and be lost are currently
    amplified normally. That biases efficiency upward by roughly the pillar
    area fraction, ~4.8 %.
  * **No electronics.** Stage C owns shaping, threshold, ENC, PDO/TDO.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .frame import DriftFrame
from .gas import Transport


@dataclass
class DigitizerConfig:
    """Everything Stage B needs that is not in the Stage A file."""

    gain: float = 1.0e4                 # mean avalanche gain
    polya_theta: float = 2.0            # Polya shape; MX17/ATLAS convention
    mesh_transparency: float = 1.0      # eps; 1.0 = placeholder, see P0.13
    mapping_revision: str = "B"
    seed: int = 12345

    # Electrons landing outside the instrumented pad area are dropped and
    # counted. If this fraction is large the beam is not where you think.
    warn_offplane_frac: float = 0.05


@dataclass
class PadHits:
    """Per-pad results for a set of events, flat arrays (one row per pad hit)."""

    event: np.ndarray
    pad: np.ndarray
    charge_e: np.ndarray        # electrons after gain
    t_first_ns: np.ndarray
    n_electrons: np.ndarray     # primaries that reached this pad (pre-gain)

    # bookkeeping
    n_events: int = 0
    n_electrons_total: int = 0
    n_attached: int = 0
    n_mesh_lost: int = 0
    n_offplane: int = 0
    provenance: dict = field(default_factory=dict)


def digitize(clusters, frame: DriftFrame, transport: Transport,
             padplane, cfg: DigitizerConfig) -> PadHits:
    """Digitize a set of clusters.

    `clusters` is a dict of equal-length arrays with keys
    eventID, x, y, z, nPrimary, time, volume.
    """
    rng = np.random.default_rng(cfg.seed)

    ev   = np.asarray(clusters["eventID"], np.int64)
    x    = np.asarray(clusters["x"], float)
    y    = np.asarray(clusters["y"], float)
    z    = np.asarray(clusters["z"], float)
    npr  = np.asarray(clusters["nPrimary"], np.int64)
    t0   = np.asarray(clusters["time"], float)
    vol  = np.asarray(clusters["volume"])

    is_amp = (vol == "AmpGas")

    # Expand clusters into individual electrons.
    keep = npr > 0
    if not np.any(keep):
        return PadHits(*(np.array([], t) for t in (np.int64, np.int64, float,
                                                   float, np.int64)))
    idx = np.repeat(np.flatnonzero(keep), npr[keep])
    n_e = idx.size

    e_ev, e_x, e_y, e_z, e_t = ev[idx], x[idx], y[idx], z[idx], t0[idx]
    e_amp = is_amp[idx]

    depth = frame.depth_mm(e_z)                   # mm, >=0 in the drift gas
    depth_drift = np.where(e_amp, 0.0, np.maximum(depth, 0.0))

    # ── attachment ────────────────────────────────────────────────────────
    if np.isfinite(transport.attach_length_mm):
        survive = rng.random(n_e) < np.exp(-depth_drift / transport.attach_length_mm)
    else:
        survive = np.ones(n_e, bool)
    n_attached = int((~survive).sum())

    # ── drift: time and transverse spread ────────────────────────────────
    v = transport.v_drift_mm_per_ns
    sqrt_d = np.sqrt(depth_drift)
    t_arr = e_t + depth_drift / v + rng.normal(
        0.0, np.maximum(transport.sigma_L_mm_per_sqrt_mm * sqrt_d / v, 1e-12))
    spread = np.maximum(transport.sigma_T_mm_per_sqrt_mm * sqrt_d, 0.0)
    x_arr = e_x + rng.normal(0.0, np.maximum(spread, 1e-12))
    y_arr = e_y + rng.normal(0.0, np.maximum(spread, 1e-12))

    # ── mesh transparency ────────────────────────────────────────────────
    # Amp-gap electrons are already past the mesh: transparency does not
    # apply to them.
    through = np.ones(n_e, bool)
    if cfg.mesh_transparency < 1.0:
        through = e_amp | (rng.random(n_e) < cfg.mesh_transparency)
    n_mesh_lost = int((~through).sum())

    alive = survive & through

    # ── avalanche gain ───────────────────────────────────────────────────
    # Polya, mean `gain`. Amp-gap electrons are amplified over only the
    # remaining fraction of the gap: G_eff = G^f with f = (gap - d)/gap.
    k = cfg.polya_theta + 1.0
    gain_mean = np.full(n_e, float(cfg.gain))
    if np.any(e_amp):
        d_in = frame.amp_depth_mm(e_z[e_amp])
        frac = np.clip(1.0 - d_in / frame.amp_gap_mm, 0.0, 1.0)
        gain_mean[e_amp] = np.power(float(cfg.gain), frac)
    q = rng.gamma(k, gain_mean / k)

    # ── pad assignment ───────────────────────────────────────────────────
    pad = padplane.locate(x_arr, y_arr)
    on = alive & (pad >= 0)
    n_offplane = int((alive & (pad < 0)).sum())

    # ── sum per (event, pad) ─────────────────────────────────────────────
    if not np.any(on):
        hits = PadHits(np.array([], np.int64), np.array([], np.int64),
                       np.array([], float), np.array([], float),
                       np.array([], np.int64))
    else:
        key_ev, key_pad = e_ev[on], pad[on]
        order = np.lexsort((key_pad, key_ev))
        key_ev, key_pad = key_ev[order], key_pad[order]
        qs, ts = q[on][order], t_arr[on][order]

        new = np.empty(key_ev.size, bool)
        new[0] = True
        new[1:] = (key_ev[1:] != key_ev[:-1]) | (key_pad[1:] != key_pad[:-1])
        starts = np.flatnonzero(new)

        charge = np.add.reduceat(qs, starts)
        counts = np.diff(np.append(starts, key_ev.size))
        tmin = np.minimum.reduceat(ts, starts)

        hits = PadHits(event=key_ev[starts], pad=key_pad[starts],
                       charge_e=charge, t_first_ns=tmin,
                       n_electrons=counts.astype(np.int64))

    hits.n_events = int(np.unique(ev).size)
    hits.n_electrons_total = n_e
    hits.n_attached = n_attached
    hits.n_mesh_lost = n_mesh_lost
    hits.n_offplane = n_offplane
    hits.provenance = dict(
        gas=frame.gas, gas_label=frame.gas_label,
        drift_gap_mm=frame.drift_gap_mm, amp_gap_um=frame.amp_gap_um,
        geometry_hash=frame.geometry_hash, git_hash=frame.git_hash,
        git_dirty=frame.git_dirty, thrown=frame.thrown,
        beam_spread_mm=frame.beam_spread_mm,
        gain=cfg.gain, polya_theta=cfg.polya_theta,
        mesh_transparency=cfg.mesh_transparency,
        mapping_revision=cfg.mapping_revision, seed=cfg.seed,
        transport=transport.to_dict(),
        stage_b_caveats=[
            "induced current not modelled (delta charge per electron) — P0.17",
            "mesh transparency is a constant, not a field map — P0.13",
            "no pillars: efficiency biased high by ~4.8% — P0.16",
            "gas transport parameters are placeholders — P0.10",
        ],
    )
    return hits
