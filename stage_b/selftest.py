#!/usr/bin/env python3
"""selftest.py — Stage B checks that do not need a Geant4 file.

    python3 -m stage_b.selftest

Synthetic clusters and a fake transport table, so the physics can be checked
against closed-form answers. This is what P0.9 means by "unit-test with fake
tables so Magboltz becomes a drop-in": every test here passes identically
whichever TransportTable is plugged in, because each one fixes the parameters
it depends on.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from stage_b.digitize import DigitizerConfig, digitize          # noqa: E402
from stage_b.frame import DriftFrame                            # noqa: E402
from stage_b.gas import MagboltzTable, PlaceholderTable, Transport  # noqa: E402
from stage_b.padplane import PadPlane                           # noqa: E402

FAILURES = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def make_frame(gap=3.0, mesh_z=8.0241):
    return DriftFrame(
        mesh_z_mm=mesh_z, drift_entry_z_mm=mesh_z - gap,
        amp_entry_z_mm=mesh_z + 0.038, pad_plane_z_mm=mesh_z + 0.038 + 0.150,
        drift_sign=1.0, drift_gap_mm=gap, amp_gap_um=150.0,
        gas="ArIso", gas_label="Ar/iC4H10 95/5", w_value_eV=26.36,
        particle="muon", energy_MeV=200000.0, thrown=1000,
        geometry_hash="test", git_hash="test", git_dirty=False,
        beam_spread_mm=0.0)


def flat_transport(v=0.05, sL=0.0, sT=0.0, attach=np.inf):
    """A transport with diffusion switched off, so tests are exact."""
    return Transport(gas="ArIso", v_drift_mm_per_ns=v,
                     sigma_L_mm_per_sqrt_mm=sL, sigma_T_mm_per_sqrt_mm=sT,
                     attach_length_mm=attach, n_p_per_cm=28.3,
                     n_tot_per_cm=103.2, v_scale=1.0, provenance="test")


def main():
    pp = PadPlane.load("B")
    frame = make_frame()
    print(f"Pad plane: {pp.n_pads} pads, {len(pp.ring_centres)} rings, "
          f"pitch {pp.ring_pitch_mm:.4f} mm")

    # ── 1. Frame sign convention ─────────────────────────────────────────
    d_cath = frame.depth_mm(frame.drift_entry_z_mm)
    d_mesh = frame.depth_mm(frame.mesh_z_mm)
    check("depth is +gap at the cathode and 0 at the mesh",
          abs(d_cath - 3.0) < 1e-9 and abs(d_mesh) < 1e-9,
          f"cathode {d_cath:.3f} mm, mesh {d_mesh:.3f} mm")
    d_amp = frame.depth_mm(frame.amp_entry_z_mm + 0.05)
    check("amp-gap clusters get NEGATIVE drift depth", d_amp < 0,
          f"{d_amp:.4f} mm")

    # ── 2. Pad lookup round-trips through pad centres ────────────────────
    x, y = pp.pad_r * np.cos(pp.pad_phi), pp.pad_r * np.sin(pp.pad_phi)
    got = pp.locate(x, y)
    check("every pad centre locates to its own pad",
          np.array_equal(got, np.arange(pp.n_pads)),
          f"{int((got == np.arange(pp.n_pads)).sum())}/{pp.n_pads}")
    off = pp.locate([0.0, 5000.0], [0.0, 5000.0])
    check("points off the pad plane return -1, not a nearest-pad guess",
          np.all(off == -1), f"{off}")

    # ── 3. Drift time is depth / v_d ─────────────────────────────────────
    tr = flat_transport(v=0.05)
    n = 4000
    z = np.full(n, frame.mesh_z_mm - 2.0)          # exactly 2 mm of drift
    cl = dict(eventID=np.arange(n), x=x[:1].repeat(n), y=y[:1].repeat(n),
              z=z, nPrimary=np.ones(n, int), time=np.zeros(n),
              volume=np.array(["DriftGas"] * n))
    hits = digitize(cl, frame, tr, pp, DigitizerConfig(seed=1))
    expect = 2.0 / 0.05
    check("drift time = depth / v_d with diffusion off",
          np.allclose(hits.t_first_ns, expect, atol=1e-9),
          f"{hits.t_first_ns[0]:.3f} ns vs {expect:.3f}")

    # ── 4. Longitudinal diffusion has the right width ────────────────────
    tr = flat_transport(v=0.05, sL=0.25 / np.sqrt(10.0))   # 250 um/sqrt(cm)
    hits = digitize(cl, frame, tr, pp, DigitizerConfig(seed=2))
    sigma_expect = tr.sigma_L_mm_per_sqrt_mm * np.sqrt(2.0) / 0.05
    got_sigma = hits.t_first_ns.std()
    check("longitudinal diffusion width matches sigma_L sqrt(z) / v_d",
          abs(got_sigma - sigma_expect) / sigma_expect < 0.08,
          f"{got_sigma:.3f} ns vs {sigma_expect:.3f}")

    # ── 5. Polya gain: mean and relative variance ────────────────────────
    tr = flat_transport()
    hits = digitize(cl, frame, tr, pp, DigitizerConfig(gain=1e4, polya_theta=2.0,
                                                       seed=3))
    q = hits.charge_e
    check("Polya mean gain is as configured",
          abs(q.mean() - 1e4) / 1e4 < 0.05, f"{q.mean():.0f} vs 1e4")
    # For Polya with shape k = theta+1, var/mean^2 = 1/k
    check("Polya relative variance is 1/(theta+1)",
          abs(q.var() / q.mean()**2 - 1/3.0) < 0.05,
          f"{q.var()/q.mean()**2:.3f} vs {1/3.0:.3f}")

    # ── 6. Attachment survival is exponential ────────────────────────────
    lam = 2.0
    tr = flat_transport(attach=lam)
    hits = digitize(cl, frame, tr, pp, DigitizerConfig(seed=4))
    surv = 1.0 - hits.n_attached / n
    check("attachment survival = exp(-z/lambda)",
          abs(surv - np.exp(-2.0 / lam)) < 0.03,
          f"{surv:.3f} vs {np.exp(-2.0/lam):.3f}")

    # ── 7. Mesh transparency ─────────────────────────────────────────────
    tr = flat_transport()
    hits = digitize(cl, frame, tr, pp,
                    DigitizerConfig(mesh_transparency=0.6, seed=5))
    frac = 1.0 - hits.n_mesh_lost / n
    check("mesh transparency passes the configured fraction",
          abs(frac - 0.6) < 0.03, f"{frac:.3f} vs 0.600")

    # ── 8. Amp-gap clusters skip the drift and get partial gain ──────────
    z_amp = np.full(n, frame.amp_entry_z_mm + 0.5 * frame.amp_gap_mm)
    cl_amp = dict(cl, z=z_amp, volume=np.array(["AmpGas"] * n))
    hits = digitize(cl_amp, frame, flat_transport(), pp,
                    DigitizerConfig(gain=1e4, seed=6))
    check("amp-gap clusters have no drift time",
          np.allclose(hits.t_first_ns, 0.0, atol=1e-9),
          f"{hits.t_first_ns[0]:.4f} ns")
    # Released halfway across the gap -> amplified over half -> G^0.5 = 100
    check("amp-gap gain is G^((gap-d)/gap), not G",
          abs(hits.charge_e.mean() - 100.0) / 100.0 < 0.15,
          f"{hits.charge_e.mean():.1f} vs 100 (= 1e4^0.5)")

    # ── 9. Charge sharing appears only near an edge ──────────────────────
    # A cluster at a pad centre must give exactly one pad; one placed on the
    # azimuthal boundary must share. This is the mechanism behind the ~1.1
    # pads/hit figure and behind the aim-point trap.
    tr = flat_transport(sT=0.25 / np.sqrt(10.0))
    ring = 20
    pads = pp._pad_by_ring[ring]
    p0, p1 = pads[5], pads[6]
    for label, (cx, cy), want_one in (
            ("pad centre", (pp.pad_r[p0]*np.cos(pp.pad_phi[p0]),
                            pp.pad_r[p0]*np.sin(pp.pad_phi[p0])), True),
            ("azimuthal boundary",
             (pp.pad_r[p0]*np.cos(0.5*(pp.pad_phi[p0]+pp.pad_phi[p1])),
              pp.pad_r[p0]*np.sin(0.5*(pp.pad_phi[p0]+pp.pad_phi[p1]))), False)):
        m = 2000
        c = dict(eventID=np.arange(m), x=np.full(m, cx), y=np.full(m, cy),
                 z=np.full(m, frame.mesh_z_mm - 3.0),
                 nPrimary=np.full(m, 20), time=np.zeros(m),
                 volume=np.array(["DriftGas"] * m))
        h = digitize(c, frame, tr, pp, DigitizerConfig(seed=7))
        per_ev = len(h.pad) / m
        if want_one:
            check(f"{label}: ~1 pad per event", per_ev < 1.02,
                  f"{per_ev:.3f} pads/event")
        else:
            check(f"{label}: charge shares across pads", per_ev > 1.5,
                  f"{per_ev:.3f} pads/event")

    # ── 10. Placeholder table is flagged as such ─────────────────────────
    t = PlaceholderTable().for_gas("NeIso9010")
    check("placeholder transport says so in its provenance",
          "PLACEHOLDER" in t.provenance, t.provenance[:48])
    t2 = PlaceholderTable().for_gas("ArCF4")
    check("non-campaign gas is flagged, not silently substituted",
          "NOT a campaign gas" in t2.provenance, t2.provenance[:48])

    # ── 11. Magboltz table: nodes, interpolation, guards ─────────────────
    mt = MagboltzTable()
    want = {"ArIso", "ArIso9010", "NeIso", "NeIso8515", "ArCO2Iso9352", "ArCF4Iso"}
    check("Magboltz table covers the six campaign / SPS gases",
          want <= set(mt.gases()), ", ".join(mt.gases()))
    raw = mt._d["ArIso"]
    k = 20
    node = mt.for_gas("ArIso", raw["E_V_per_cm"][k])
    check("Magboltz value at a table node is the table value",
          abs(node.v_drift_mm_per_ns * 1e3 - raw["v_um_per_ns"][k]) < 1e-6 * raw["v_um_per_ns"][k],
          f"{node.v_drift_mm_per_ns*1e3:.3f} um/ns")
    lo = mt.for_gas("ArIso", raw["E_V_per_cm"][k]).v_drift_mm_per_ns
    hi = mt.for_gas("ArIso", raw["E_V_per_cm"][k + 1]).v_drift_mm_per_ns
    mid = mt.for_gas("ArIso", np.sqrt(raw["E_V_per_cm"][k] * raw["E_V_per_cm"][k + 1])).v_drift_mm_per_ns
    check("Magboltz interpolation lies between the neighbouring nodes",
          min(lo, hi) <= mid <= max(lo, hi))
    def raises(f):
        try:
            f()
        except (ValueError, KeyError):
            return True
        return False
    check("Magboltz refuses a missing drift field", raises(lambda: mt.for_gas("ArIso")))
    check("Magboltz refuses a field outside the table",
          raises(lambda: mt.for_gas("ArIso", 10.0)))
    check("Magboltz refuses a gas it has no table for",
          raises(lambda: mt.for_gas("NeIso9010", 625.0)))
    check("Magboltz transport says so in its provenance",
          "Magboltz" in node.provenance, node.provenance[:48])

    print()
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s): {', '.join(FAILURES)}")
        return 1
    print("All Stage B self-tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
