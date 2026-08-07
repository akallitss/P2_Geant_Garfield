"""VMM3a channel emulator for the P2 wedge pads — Stage C of the campaign.

Consumes per-pad electron arrival lists produced by the Stage B response
sampler (drift + diffusion + gain), applies the VMM front-end model, and
returns digitized hits. See README.md for the model tiers, parameter
provenance and the open calibration questions (flagged NEEDS-DATA).

Tiers implemented here:
  T0 'charge'  — ballistic-scaled charge sum per pad vs threshold.
                 No timing. Fast; the campaign workhorse.
  T1 'shaper'  — full Athena closed-form shaper: first-peak amplitude
                 (PDO proxy) and threshold-crossing time (TDO proxy),
                 neighbor logic, time-window acceptance.

Not implemented (Phase 4, only if needed): TAC/BC timing quantization,
ion-tail suppression stages, per-channel trim spread, dead-time pileup.

Units: charge in ELECTRONS everywhere (Athena convention); helpers convert
to fC (1 fC = 6241.5 e-). Time in ns.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from vmm_shaper import ELECTRONS_PER_FC, ShaperOutput, VMMShaper


@dataclass
class VMMConfig:
    """One VMM operating point.

    Defaults follow the ATLAS Athena MM digitization where P2/BASKET values
    are not yet known; every NEEDS-DATA field should eventually come from
    the BASKET prototype configuration.
    """

    peak_time_ns: float = 200.0        # VMM3a: 25/50/100/200
    gain_mv_per_fc: float = 9.0        # VMM3a: 0.5/1/3/4.5/6/9/12/16 (NEEDS-DATA)
    threshold_electrons: float = 15000.0  # Athena default global threshold
    enc_electrons: float = 0.0         # Gaussian noise added to amplitude;
                                       # 0 = noiseless. (NEEDS-DATA: ENC at
                                       # 75-170 pF pad capacitance)
    # VMM3a neighbor logic. Default OFF: P2 does not plan to run with NL, and
    # on this pad plane NL reads *chip-channel* neighbors, which are only
    # partly the physical pad neighbors (vmm/nl_map.py, README §3.1). Turn on
    # only as a scan variant, with a channel-based neighbor map.
    neighbor_logic: bool = False
    nl_readout_threshold_electrons: float = 1.0  # Athena: NL strips re-run
                                                 # with threshold 1 e-
    # acceptance window for peak/threshold time (Athena: -12.5 .. +187.5 ns
    # around the BC for NSW; P2 streaming wants a wide-open window)
    time_window_ns: Tuple[float, float] = (-50.0, 500.0)
    vmm_deadtime_ns: float = 200.0     # search-window extension below
    vmm_upper_graze_ns: float = 150.0  # search-window extension above
    # ADC model (None = no quantization). PDO full scale in mV; the shaped
    # amplitude in mV is amplitude_fC * gain_mv_per_fc. (NEEDS-DATA:
    # measured PDO calibration; ~1000-1200 mV typical, quote-only for now)
    adc_bits: int = 10
    adc_full_scale_mv: Optional[float] = None
    rng: Optional[np.random.Generator] = None

    def shaper(self) -> VMMShaper:
        lo, hi = self.time_window_ns
        return VMMShaper(self.peak_time_ns,
                         lo - self.vmm_deadtime_ns,
                         hi + self.vmm_upper_graze_ns)

    @property
    def full_scale_electrons(self) -> Optional[float]:
        if self.adc_full_scale_mv is None:
            return None
        return (self.adc_full_scale_mv / self.gain_mv_per_fc) * ELECTRONS_PER_FC


@dataclass
class PadHit:
    """One digitized pad hit."""

    pad: int
    charge_electrons: float      # PDO proxy (after noise/saturation/ADC)
    time_ns: float               # TDO proxy (threshold-crossing time); nan in T0
    above_threshold: bool        # False = read out via neighbor logic only
    saturated: bool = False

    @property
    def charge_fc(self) -> float:
        return self.charge_electrons / ELECTRONS_PER_FC


# Input: {pad_id: (times_ns, charges_electrons)} for one event.
PadPulses = Dict[int, Tuple[Sequence[float], Sequence[float]]]
# Pad adjacency: pad_id -> iterable of neighbor pad ids.
#
# For NEIGHBOR LOGIC this must be the *channel* neighbor set — VMM NL fires on
# chan +-1 within a chip, and on the P2 pad plane those are only partly the
# physical pad neighbors (74 % touching, 25 % of the true neighborhood, in the
# 9-column mapping revision; 8 % / 3 % in the 7-column one). Use
# `vmm.nl_map.PadMap.channel_neighbors`. Passing geometric neighbors here
# models a detector we do not have.
#
# For CLUSTERING of pads that each passed threshold, use the geometric
# neighbors (`PadMap.geometric_neighbors`) — that is a separate, offline step
# and is unaffected by the NL setting.
NeighborMap = Callable[[int], Sequence[int]]


class PadDigitizer:
    def __init__(self, config: VMMConfig,
                 neighbor_map: Optional[NeighborMap] = None):
        self.cfg = config
        self.neighbors = neighbor_map
        self._shaper = config.shaper()
        self._rng = config.rng or np.random.default_rng()

    # ------------------------------------------------------------------
    def digitize(self, pulses: PadPulses, tier: str = "shaper") -> List[PadHit]:
        if tier == "charge":
            return self._digitize_charge(pulses)
        if tier == "shaper":
            return self._digitize_shaper(pulses)
        raise ValueError(f"unknown tier: {tier}")

    # ---------------- T0 ------------------------------------------------
    def _digitize_charge(self, pulses: PadPulses) -> List[PadHit]:
        hits: List[PadHit] = []
        scale = self._shaper.peak_time_charge_scaling
        for pad, (times, charges) in pulses.items():
            q = float(np.sum(charges)) * scale
            q = self._add_noise(q)
            if q < self.cfg.threshold_electrons:
                continue
            q, sat = self._quantize(q)
            hits.append(PadHit(pad, q, float("nan"), True, sat))
        return hits

    # ---------------- T1 ------------------------------------------------
    def _digitize_shaper(self, pulses: PadPulses) -> List[PadHit]:
        thr = self.cfg.threshold_electrons
        fired: Dict[int, bool] = {}
        for pad, (times, charges) in pulses.items():
            fired[pad] = self._shaper.has_charge_above_threshold(
                charges, times, thr)

        hits: List[PadHit] = []
        lo, hi = self.cfg.time_window_ns
        for pad, (times, charges) in pulses.items():
            read_out = fired[pad]
            via_nl = False
            if not read_out and self.cfg.neighbor_logic and self.neighbors:
                via_nl = any(fired.get(n, False) for n in self.neighbors(pad))
                read_out = via_nl
            if not read_out:
                continue
            # NL-only pads are re-evaluated with a ~zero threshold (Athena)
            eff_thr = thr if not via_nl else self.cfg.nl_readout_threshold_electrons
            out: ShaperOutput = self._shaper.threshold_response(
                charges, times, eff_thr)
            if not out.fired:
                continue
            if not (lo <= out.peak_time <= hi):
                continue  # Athena time-window acceptance
            q = self._add_noise(out.amplitude)
            q, sat = self._quantize(q)
            hits.append(PadHit(pad, q, out.threshold_time, not via_nl, sat))
        return hits

    # ------------------------------------------------------------------
    def _add_noise(self, q: float) -> float:
        if self.cfg.enc_electrons > 0:
            q += self._rng.normal(0.0, self.cfg.enc_electrons)
        return q

    def _quantize(self, q: float) -> Tuple[float, bool]:
        fs = self.cfg.full_scale_electrons
        if fs is None:
            return q, False
        sat = q >= fs
        q = min(max(q, 0.0), fs)
        lsb = fs / (2 ** self.cfg.adc_bits)
        return np.floor(q / lsb) * lsb, sat


# ----------------------------------------------------------------------
if __name__ == "__main__":
    # smoke test: 3 pads, middle pad big signal, right pad small (NL-only)
    cfg = VMMConfig(threshold_electrons=15000, enc_electrons=0,
                    neighbor_logic=True)
    adjacency = {1: [2], 2: [1, 3], 3: [2]}
    dig = PadDigitizer(cfg, neighbor_map=lambda p: adjacency.get(p, []))
    ev: PadPulses = {
        2: ([50.0, 60.0, 70.0], [30000.0, 25000.0, 20000.0]),
        3: ([55.0], [4000.0]),
    }
    for tier in ("charge", "shaper"):
        hits = dig.digitize(ev, tier=tier)
        print(f"tier={tier}:")
        for h in sorted(hits, key=lambda h: h.pad):
            t = "nan" if np.isnan(h.time_ns) else f"{h.time_ns:.1f}"
            print(f"  pad {h.pad}: q={h.charge_electrons:.0f} e- "
                  f"({h.charge_fc:.2f} fC), t={t} ns, "
                  f"above_thr={h.above_threshold}")
    assert any(h.pad == 3 for h in dig.digitize(ev, tier="shaper")), \
        "neighbor logic should read out pad 3"
    assert not any(h.pad == 3 for h in dig.digitize(ev, tier="charge")), \
        "T0 has no neighbor logic"
    print("vmm_emulator smoke test OK")
