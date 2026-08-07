"""Expected per-pad time resolution of the P2 Micromegas + VMM3a.

Standalone toy Monte Carlo (no Geant4, no Stage B) that predicts the width of
the VMM **TDO** (threshold-crossing time) for a through-going MIP, so the
prediction can be compared with the SPS muon test-beam timing distributions
(docs/TESTBEAM_PLAN.md §2.6). Full write-up, numbers and caveats:
docs/research/TIME_RESOLUTION_NOTES.md.

Chain per event:
  1. primary clusters along the track: Poisson(n_p * d / cos(theta)), uniform
     in z; electrons per cluster from p(n) ~ n^-2 truncated to match n_tot/n_p;
  2. drift to the mesh: t = z / v_d + N(0, sigma_L * sqrt(z) / v_d);
  3. avalanche: Polya(theta) gain per electron, mean G;
  4. all charge on one pad (pad multiplicity ~1.1 -- see §NL discussion in
     vmm/README.md); optional charge sharing fraction for the 2-pad case;
  5. Athena VMM shaper -> threshold crossing t_thr and peak amplitude (PDO);
  6. electronics jitter added analytically: sigma_t = ENC / |dV/dt|_thr;
  7. TDO quantization: TAC ramp / 2^8, plus the coarse BC clock.

Deliberate simplifications (all conservative-to-optimistic, flagged):
  - each arriving electron is a delta charge into the shaper (Athena's model).
    The real induced current is a fast electron spike plus a ~150 ns ion tail;
    a delta gives a slightly *faster* leading edge, so the electronics-jitter
    term here is a lower bound (the Athena ballistic-deficit factor s_bd
    handles the amplitude side, not the slope).
  - noise is white at the shaper output and enters only through the slope;
    no correlated baseline wander, no ion-tail suppression stages.
  - no time-walk from threshold-DAC trim spread between channels (a fixed
    per-channel offset, calibrated out in data).

Usage:
    python vmm/time_resolution.py                      # full table
    python vmm/time_resolution.py --gas ArCO2Iso_93_5_2 --tp 200 -n 20000
    python vmm/time_resolution.py --breakdown          # contribution budget
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from vmm_shaper import ELECTRONS_PER_FC, VMMShaper


# ----------------------------------------------------------------------
# Gas parameters.
#
# n_p / n_tot: volume-weighted from PDG 2024 Table 35.5 (same source as
# SIM_CAMPAIGN_PLAN.md §1).
# t_drift_3mm / sigma_L: PLACEHOLDERS pending the Magboltz tables (campaign
# P0.10). The drift times are the ones quoted in the campaign docs (Ar-based
# 60-75 ns, Ne/CO2 ~110 ns over 3 mm); sigma_L is a generic 230-300 um/sqrt(cm).
# Every result below scales as 1/(n_p * v_d), so these two numbers are the
# dominant input uncertainty -- re-run once real tables exist.
# ----------------------------------------------------------------------
@dataclass
class Gas:
    name: str
    n_p_per_cm: float          # primary clusters / cm
    n_tot_per_cm: float        # total ionization electrons / cm
    t_drift_3mm_ns: float      # full-gap drift time at 3 mm (placeholder)
    sigma_L_um_sqrtcm: float   # longitudinal diffusion (placeholder)

    @property
    def v_drift_mm_per_ns(self) -> float:
        return 3.0 / self.t_drift_3mm_ns

    @property
    def mean_cluster_size(self) -> float:
        return self.n_tot_per_cm / self.n_p_per_cm


GASES: Dict[str, Gas] = {
    "ArIso_95_5":       Gas("Ar/iC4H10 95/5",      28.3, 103.2,  60.0, 250.0),
    "ArCO2Iso_93_5_2":  Gas("Ar/CO2/iC4H10 93/5/2", 26.8,  99.6,  65.0, 230.0),
    "NeIso_90_10":      Gas("Ne/iC4H10 90/10",     20.7,  58.0,  75.0, 250.0),
    "NeCO2Iso_95_3_2":  Gas("Ne/CO2/iC4H10 95/3/2", 15.2,  45.4, 110.0, 300.0),
    "NeCH4_93_7":       Gas("Ne/CH4 93/7",         14.1,  41.0,  70.0, 280.0),
}

# VMM3a TAC ramps (ns) -> TDO LSB = ramp / 2^8
TAC_RAMPS_NS = (60.0, 100.0, 350.0, 650.0)


# ----------------------------------------------------------------------
def cluster_size_sampler(mean_size: float, n_max_cap: int = 400):
    """p(n) ~ n^-2 on [1, N], N tuned so <n> == mean_size (Fischle-style tail)."""
    best, best_err = None, 1e9
    for N in range(2, n_max_cap + 1):
        n = np.arange(1, N + 1, dtype=float)
        w = n ** -2.0
        m = (n * w).sum() / w.sum()
        err = abs(m - mean_size)
        if err < best_err:
            best, best_err = (n, w / w.sum()), err
        if m > mean_size:
            break
    n, p = best
    return n.astype(int), p


def polya_gain(rng: np.random.Generator, n: int, mean_gain: float,
               theta: float) -> np.ndarray:
    """Polya (Gamma) avalanche gain, shape k = theta+1, mean = mean_gain."""
    k = theta + 1.0
    return rng.gamma(k, mean_gain / k, size=n)


# ----------------------------------------------------------------------
@dataclass
class TimingConfig:
    gas: Gas
    gap_mm: float = 3.0
    angle_deg: float = 0.0
    gain: float = 1.0e4
    polya_theta: float = 2.0
    peak_time_ns: float = 200.0
    thresholds_fc: Sequence[float] = (1.0, 2.0, 5.0, 10.0)
    enc_electrons: Sequence[float] = (0.0, 1500.0, 3000.0, 6000.0)
    tac_ramp_ns: Optional[float] = 100.0     # None = no TDO quantization
    tdo_bits: int = 8
    charge_share: float = 1.0                # fraction of charge on the pad
    diffusion: bool = True
    gain_fluct: bool = True
    grid_step_ns: float = 0.2
    grid_span_ns: Tuple[float, float] = (-20.0, 500.0)


def simulate(cfg: TimingConfig, n_events: int = 5000,
             seed: int = 12345) -> Dict:
    """Return per-(threshold, ENC) timing statistics for one configuration."""
    rng = np.random.default_rng(seed)
    g = cfg.gas
    path_mm = cfg.gap_mm / math.cos(math.radians(cfg.angle_deg))
    mean_clusters = g.n_p_per_cm * path_mm / 10.0
    sizes, size_p = cluster_size_sampler(g.mean_cluster_size)
    v = g.v_drift_mm_per_ns
    sigma_L_mm_sqrtmm = g.sigma_L_um_sqrtcm * 1e-3 / math.sqrt(10.0)

    shaper = VMMShaper(cfg.peak_time_ns, cfg.grid_span_ns[0], cfg.grid_span_ns[1],
                       time_step=cfg.grid_step_ns)
    grid = np.arange(cfg.grid_span_ns[0], cfg.grid_span_ns[1], cfg.grid_step_ns)

    thr_e = np.array([t * ELECTRONS_PER_FC for t in cfg.thresholds_fc])
    n_thr = len(thr_e)
    t_cross = [[] for _ in range(n_thr)]     # noiseless crossing time
    slope = [[] for _ in range(n_thr)]       # dV/dt at crossing (e-/ns)
    amp = [[] for _ in range(n_thr)]         # PDO proxy (peak amplitude)
    t_peak = [[] for _ in range(n_thr)]      # time-at-peak (the other TDO mode)
    n_fired = np.zeros(n_thr, dtype=int)

    for _ in range(n_events):
        n_cl = rng.poisson(mean_clusters)
        if n_cl == 0:
            continue
        z = rng.uniform(0.0, cfg.gap_mm, size=n_cl)       # 0 = at the mesh
        mult = rng.choice(sizes, size=n_cl, p=size_p)
        z_e = np.repeat(z, mult)
        t_e = z_e / v
        if cfg.diffusion:
            t_e = t_e + rng.normal(0.0, sigma_L_mm_sqrtmm * np.sqrt(np.maximum(z_e, 1e-6)) / v)
        q_e = (polya_gain(rng, len(z_e), cfg.gain, cfg.polya_theta)
               if cfg.gain_fluct else np.full(len(z_e), cfg.gain))
        q_e = q_e * cfg.charge_share

        # one shaped waveform per event, reused for every threshold
        wf = shaper.response(q_e, t_e, grid)
        ipk = int(np.argmax(wf))
        peak = float(wf[ipk])
        tpk = float(grid[ipk])
        for k in range(n_thr):
            above = np.flatnonzero(wf >= thr_e[k])
            if above.size == 0 or peak < thr_e[k]:
                continue
            i = int(above[0])
            if i == 0:
                continue
            # linear interpolation of the crossing between grid samples
            y0, y1 = wf[i - 1], wf[i]
            frac = (thr_e[k] - y0) / max(y1 - y0, 1e-9)
            t_c = grid[i - 1] + frac * cfg.grid_step_ns
            s = (y1 - y0) / cfg.grid_step_ns
            if s <= 0:
                continue
            n_fired[k] += 1
            t_cross[k].append(t_c)
            slope[k].append(s)
            amp[k].append(peak)
            t_peak[k].append(tpk)

    out = {"config": _cfg_summary(cfg), "n_events": n_events, "points": []}
    for k, thr in enumerate(cfg.thresholds_fc):
        tc = np.array(t_cross[k]); sl = np.array(slope[k]); am = np.array(amp[k])
        tp = np.array(t_peak[k])
        if tc.size < 50:
            out["points"].append({"threshold_fc": thr, "efficiency": n_fired[k] / n_events,
                                  "n": int(tc.size), "insufficient_stats": True})
            continue
        for enc in cfg.enc_electrons:
            jitter = rng.normal(0.0, enc / sl) if enc > 0 else np.zeros_like(tc)
            t_meas = tc + jitter
            if cfg.tac_ramp_ns:
                lsb = cfg.tac_ramp_ns / (2 ** cfg.tdo_bits)
                t_meas = np.floor(t_meas / lsb) * lsb
            out["points"].append({
                "threshold_fc": thr,
                "enc_electrons": enc,
                "efficiency": float(n_fired[k] / n_events),
                "n": int(tc.size),
                "mean_t_ns": float(np.mean(t_meas)),
                "sigma_raw_ns": _sigma68(t_meas),
                "rms_raw_ns": float(np.std(t_meas)),
                "sigma_walkcorr_ns": _sigma68(_walk_correct(t_meas, am)),
                "rms_walkcorr_ns": float(np.std(_walk_correct(t_meas, am))),
                "median_slope_e_per_ns": float(np.median(sl)),
                "median_pdo_fc": float(np.median(am) / ELECTRONS_PER_FC),
                # the other VMM TDO mode: timestamp at the peak instead of at
                # the threshold crossing (per-chip configuration)
                "sigma_peakmode_ns": _sigma68(tp),
                "sigma_peakmode_walkcorr_ns": _sigma68(_walk_correct(tp, am)),
            })
    return out


def _cfg_summary(cfg: TimingConfig) -> Dict:
    return {"gas": cfg.gas.name, "gap_mm": cfg.gap_mm, "angle_deg": cfg.angle_deg,
            "gain": cfg.gain, "polya_theta": cfg.polya_theta,
            "peak_time_ns": cfg.peak_time_ns, "t_drift_ns": cfg.gas.t_drift_3mm_ns,
            "n_p_per_cm": cfg.gas.n_p_per_cm, "tac_ramp_ns": cfg.tac_ramp_ns,
            "charge_share": cfg.charge_share, "diffusion": cfg.diffusion,
            "gain_fluct": cfg.gain_fluct}


def _sigma68(x: np.ndarray) -> float:
    """Robust Gaussian-core width: half the central 68 % interval."""
    lo, hi = np.percentile(x, [15.865, 84.135])
    return float((hi - lo) / 2.0)


def _walk_correct(t: np.ndarray, amplitude: np.ndarray, n_bins: int = 12) -> np.ndarray:
    """Subtract the median t vs PDO trend -- what an offline analysis does with
    the streamed (PDO, TDO) pair. Non-parametric: quantile bins in log(PDO)."""
    if t.size < 5 * n_bins:
        n_bins = max(2, t.size // 20)
    edges = np.percentile(amplitude, np.linspace(0, 100, n_bins + 1))
    edges[0] -= 1e-6; edges[-1] += 1e-6
    idx = np.clip(np.digitize(amplitude, edges) - 1, 0, n_bins - 1)
    corr = t.copy()
    for b in range(n_bins):
        m = idx == b
        if m.sum() >= 3:
            corr[m] -= np.median(t[m])
    return corr


# ----------------------------------------------------------------------
def analytic_first_cluster_sigma(gas: Gas) -> float:
    """sigma of the arrival time of the first primary cluster = 1/(n_p * v_d).

    The classic Micromegas leading-edge limit: cluster positions are Poisson
    along z, so the distance from the mesh to the first cluster is exponential
    with mean 1/n_p, and its arrival time inherits sigma = mean."""
    return (10.0 / gas.n_p_per_cm) / gas.v_drift_mm_per_ns   # mm/ (mm/ns)


def print_table(results: List[Dict], enc: float, threshold_fc: float) -> None:
    print(f"\n{'gas':<24}{'t_p':>5}{'T_drift':>9}{'1/(np*vd)':>11}"
          f"{'eff':>7}{'<t>':>8}{'sig_raw':>9}{'sig_corr':>10}")
    print(f"{'':<24}{'[ns]':>5}{'[ns]':>9}{'[ns]':>11}{'':>7}{'[ns]':>8}"
          f"{'[ns]':>9}{'[ns]':>10}")
    print("-" * 83)
    for r in results:
        c = r["config"]
        pts = [p for p in r["points"]
               if p.get("enc_electrons") == enc and p["threshold_fc"] == threshold_fc]
        if not pts:
            continue
        p = pts[0]
        gas = next(g for g in GASES.values() if g.name == c["gas"])
        print(f"{c['gas']:<24}{c['peak_time_ns']:>5.0f}{c['t_drift_ns']:>9.0f}"
              f"{analytic_first_cluster_sigma(gas):>11.1f}{p['efficiency']:>7.3f}"
              f"{p['mean_t_ns']:>8.1f}{p['sigma_raw_ns']:>9.2f}"
              f"{p['sigma_walkcorr_ns']:>10.2f}")


def breakdown(gas_key: str, n_events: int, seed: int = 7) -> None:
    """Contribution budget: switch effects on one at a time."""
    gas = GASES[gas_key]
    thr, enc_full = 2.0, 3000.0
    stages = [
        ("ionization statistics only",     dict(diffusion=False, gain_fluct=False,
                                                enc_electrons=(0.0,), tac_ramp_ns=None)),
        ("+ longitudinal diffusion",       dict(diffusion=True, gain_fluct=False,
                                                enc_electrons=(0.0,), tac_ramp_ns=None)),
        ("+ Polya gain fluctuations",      dict(diffusion=True, gain_fluct=True,
                                                enc_electrons=(0.0,), tac_ramp_ns=None)),
        ("+ electronics noise (3 ke- ENC)", dict(diffusion=True, gain_fluct=True,
                                                enc_electrons=(enc_full,), tac_ramp_ns=None)),
        ("+ TDO quantization (100 ns TAC)", dict(diffusion=True, gain_fluct=True,
                                                enc_electrons=(enc_full,), tac_ramp_ns=100.0)),
    ]
    print(f"\n=== contribution budget: {gas.name}, 3 mm, t_p = 200 ns, "
          f"threshold {thr} fC, gain 1e4 ===")
    print(f"{'stage':<34}{'sigma_raw [ns]':>16}{'sigma_walkcorr [ns]':>21}")
    print("-" * 71)
    for label, kw in stages:
        cfg = TimingConfig(gas=gas, thresholds_fc=(thr,), **kw)
        r = simulate(cfg, n_events=n_events, seed=seed)
        p = r["points"][0]
        print(f"{label:<34}{p['sigma_raw_ns']:>16.2f}{p['sigma_walkcorr_ns']:>21.2f}")
    print(f"\nanalytic first-cluster limit 1/(n_p*v_d) = "
          f"{analytic_first_cluster_sigma(gas):.1f} ns")


# ----------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gas", default=None, choices=sorted(GASES), help="single gas")
    ap.add_argument("--tp", type=float, nargs="*", default=[25.0, 50.0, 100.0, 200.0])
    ap.add_argument("--gap", type=float, default=3.0)
    ap.add_argument("--angle", type=float, default=0.0)
    ap.add_argument("--gain", type=float, default=1.0e4)
    ap.add_argument("-n", "--n-events", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--enc", type=float, nargs="*", default=[0.0, 1500.0, 3000.0, 6000.0])
    ap.add_argument("--thresholds", type=float, nargs="*", default=[1.0, 2.0, 5.0, 10.0])
    ap.add_argument("--breakdown", action="store_true", help="contribution budget only")
    ap.add_argument("--json", default=None, help="write full results here")
    args = ap.parse_args()

    if args.breakdown:
        breakdown(args.gas or "ArCO2Iso_93_5_2", args.n_events, args.seed)
        return

    keys = [args.gas] if args.gas else sorted(GASES)
    results = []
    for key in keys:
        for tp in args.tp:
            cfg = TimingConfig(gas=GASES[key], gap_mm=args.gap, angle_deg=args.angle,
                               gain=args.gain, peak_time_ns=tp,
                               thresholds_fc=tuple(args.thresholds),
                               enc_electrons=tuple(args.enc))
            results.append(simulate(cfg, n_events=args.n_events, seed=args.seed))
            print(f"  done: {key} t_p={tp:.0f} ns", flush=True)

    for enc in args.enc:
        for thr in args.thresholds:
            print(f"\n##### ENC = {enc:.0f} e-, threshold = {thr:.1f} fC #####")
            print_table(results, enc, thr)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(results, fh, indent=1)
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
