#!/usr/bin/env python3
"""Does the Magboltz drift velocity reproduce the SPS hit-time shift? (2026-10-15)

Validation step 1 of the chain. In drift_mesh_scan_1 (Ar/CO2/iso 93/5/2,
mesh 450 V, drift 500-900 V, 4 mm) P2_IN (HV fixed) stays at 280 +- 0.5 ns,
so any change in P2_MID / P2_OUT is the drift of the ionisation column. The
DREAM `time` is a 30 % leading edge on 60 ns samples (WaveformAnalyzer.cpp,
timingPercentMax = 0.3). Data = largest good hit per event, 200-3500 ADC,
from sps_hit_times_extract.py (run on lxplus):
  median time vs field      -> shift of the charge arrival (needs v_d)
  (MID - OUT) 68 % width/V2 -> single-chamber spread, trigger cancelled

Toy forward model per drift field E = (V - 450) / 0.4 V/cm:
  MIP clusters        Poisson, N_CL /cm, uniform over the gap, cluster size
                      from a 1/n^2-tailed distribution
  drift               t = z / v_d(E) + N(0, D_L(E) sqrt(z) / v_d)  (Magboltz)
  avalanche           Polya (theta 0.5) per electron
  induced current     electron part (1/ln G) instantaneous + ion part: box of
                      duration (g - s) / v_ion, s ~ Exp(1/alpha) from the anode
  DREAM               CR-RC^n, peaking TP_NS, 16 samples x 60 ns, random phase,
                      noise, then the same 30 % leading-edge algorithm.
The data's absolute time offset is free (one constant per chamber); the
shape of t(E) is the test. v_d is then scaled by k and k fitted.

    python3 scripts/analysis/sps_drift_velocity_check.py [--n 4000]
"""
import argparse
import csv
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
DATA = pathlib.Path("/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/"
                    "mpgd26_workspace/products/analysis")
OUT = REPO / "docs" / "results" / "drift_velocity"
GAS_JSON = REPO / "stage_b" / "gas_tables" / "magboltz_dry_20261008.json"
GAS = "ArCO2Iso9352"
D_UM = 3964.0            # drift gas of the 4 mm frame (same as the contamination fit)
G_UM = 150.0             # amplification gap
SAMPLE_NS = 60.0
N_SAMPLES = 16
VOLTS = [500, 550, 600, 650, 700, 750, 800, 850, 900]


def load_gas():
    t = json.load(open(GAS_JSON))[GAS]
    E = np.array(t["E_V_per_cm"])
    return lambda key, x: np.interp(np.log(x), np.log(E), np.array(t[key], float))


def load_data():
    """Per drift point: MID / OUT median time and median amplitude, and the
    single-chamber width from the MID - OUT difference."""
    d = json.load(open(OUT / "hit_times_drift_mesh_scan_1.json"))
    out = {}
    for det in ("P2_MID", "P2_OUT"):
        out[det] = np.array([d[f"drift_{v}"][det]["time"]["p50"] for v in VOLTS])
        out[det + "_amp"] = np.array([d[f"drift_{v}"][det]["amp"]["p50"] for v in VOLTS])
    m = [d[f"drift_{v}"]["MID_minus_OUT"] for v in VOLTS]
    out["width"] = np.array([(x["p84"] - x["p16"]) / 2 / np.sqrt(2) for x in m])
    return out


def cluster_sizes(rng, n):
    # ~65 % single electrons, 1/n^2 tail to 50 (mean ~2.6): Heed-like for Ar
    k = np.arange(1, 51)
    p = np.where(k == 1, 0.65, 0.35 / k ** 2 / np.sum(1 / k[1:] ** 2))
    p /= p.sum()
    return rng.choice(k, size=n, p=p)


def response(t, tp, n_shape, t_ion, f_e):
    """Shaped response to one avalanche at t = 0 (arbitrary units)."""
    dt = t[1] - t[0]
    tau = tp / n_shape
    h = (t / tau) ** n_shape * np.exp(-t / tau)
    # ion current: average over the start depth s (from the anode), box of
    # length (1 - s/g) * t_ion; s ~ Exp with mean g / ln G (~ 15 um)
    s = np.random.default_rng(1).exponential(G_UM / 9.8, 4000)
    s = s[s < G_UM]
    dur = (1 - s / G_UM) * t_ion
    i_ion = np.mean((t[None, :] < dur[:, None]) / dur[:, None], axis=0) * dt
    i = (1 - f_e) * i_ion
    i[0] += f_e
    return np.convolve(i, h)[: len(t)]


def frac_time(w, frac=0.3):
    """WaveformAnalyzer: parabola peak, walk left to frac*peak, interpolate."""
    n, m = w.shape
    imax = w.argmax(axis=1)
    out = np.full(n, np.nan)
    for k in range(n):
        i = imax[k]
        a = w[k, i]
        if 0 < i < m - 1:
            y1, y2, y3 = w[k, i - 1], w[k, i], w[k, i + 1]
            den = y1 - 2 * y2 + y3
            if abs(den) > 1e-9:
                off = 0.5 * (y1 - y3) / den
                a = y2 - 0.25 * (y1 - y3) * off
        lv = frac * a
        j = i
        while j > 0 and w[k, j] > lv:
            j -= 1
        if j + 1 < m and abs(w[k, j + 1] - w[k, j]) > 1e-9:
            out[k] = j + (lv - w[k, j]) / (w[k, j + 1] - w[k, j])
        else:
            out[k] = j + 0.5
    return out * SAMPLE_NS


def simulate(E, vd, dl, n_ev, rng, *, n_cl, tp, n_shape, t_ion, noise, t_start,
             amp_med=300.0, amp_cut=(200.0, 3500.0)):
    """Mean and robust sigma of the 30 % time for n_ev MIPs at field E."""
    dt = 1.0
    t = np.arange(0, 3000, dt)
    h = response(t, tp, n_shape, t_ion, f_e=1 / 9.8)
    times = []
    amps = []
    for _ in range(n_ev):
        ncl = rng.poisson(n_cl * D_UM * 1e-4)
        z = rng.uniform(0, D_UM, ncl)
        sz = cluster_sizes(rng, ncl)
        z = np.repeat(z, sz)
        ta = z / vd + rng.normal(0, 1, len(z)) * dl * np.sqrt(z * 1e-4) / vd
        q = rng.gamma(1.5, 1 / 1.5, len(z))
        hist = np.bincount(np.clip((ta / dt).astype(int), 0, len(t) - 1),
                           weights=q, minlength=len(t))
        wf = np.convolve(hist, h)[: len(t)]
        phase = rng.uniform(0, SAMPLE_NS)
        idx = ((t_start + phase + np.arange(N_SAMPLES) * SAMPLE_NS) / dt).astype(int)
        s = np.where((idx >= 0) & (idx < len(t)), wf[np.clip(idx, 0, len(t) - 1)], 0)
        times.append(phase)
        amps.append(s)
    w = np.array(amps)
    # scale so the median pulse height matches the data at this field, then
    # apply the data's amplitude window (noise ~10 ADC)
    scale = amp_med / np.median(w.max(axis=1))
    w = w * scale + rng.normal(0, noise, w.shape)
    pk = w.max(axis=1)
    keep = (pk >= amp_cut[0]) & (pk <= amp_cut[1])
    # sample 0 sits at t_start + phase on the physical (track-crossing) axis
    tf = frac_time(w[keep]) + t_start + np.array(times)[keep]
    tf = tf[np.isfinite(tf)]
    med = np.median(tf)
    q16, q84 = np.percentile(tf, [16, 84])
    return np.mean(tf), med, 0.5 * (q84 - q16)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=3000)
    a = ap.parse_args()
    gas = load_gas()
    dat = load_data()
    E = (np.array(VOLTS) - 450.0) / 0.4
    vd = gas("v_um_per_ns", E)
    dl = gas("DL_um_sqrtcm", E)
    amp = 0.5 * (dat["P2_MID_amp"] + dat["P2_OUT_amp"])
    nominal = dict(n_cl=25.0, tp=180.0, n_shape=4, t_ion=300.0, noise=10.0, t_start=-200.0)
    variants = {
        "nominal": {},
        "N_CL 18": dict(n_cl=18.0),
        "N_CL 32": dict(n_cl=32.0),
        "t_ion 150 ns": dict(t_ion=150.0),
        "t_ion 500 ns": dict(t_ion=500.0),
        "CR-RC^2": dict(n_shape=2),
        "noise 20 ADC": dict(noise=20.0),
    }
    scales = [0.6, 0.8, 1.0, 1.2]
    rows, res = [], {}
    for name, ov in variants.items():
        par = {**nominal, **ov}
        for k in (scales if name == "nominal" else [1.0]):
            rng = np.random.default_rng(7)
            m = np.array([simulate(e, k * v, d, a.n, rng, amp_med=am, **par)
                          for e, v, d, am in zip(E, vd, dl, amp)])
            res[(name, k)] = m
            for i, v in enumerate(VOLTS):
                rows.append(dict(variant=name, vd_scale=k, drift_v=v, E=E[i],
                                 vd_um_ns=k * vd[i], mean_ns=m[i, 0],
                                 median_ns=m[i, 1], sigma68_ns=m[i, 2]))
            print(f"{name:13s} k={k}  median " + " ".join(f"{x:.0f}" for x in m[:, 1])
                  + "  | w68 " + " ".join(f"{x:.0f}" for x in m[:, 2]), flush=True)
    with open(OUT / "drift_time_sim.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    ref = VOLTS.index(850)
    print("data  median MID " + " ".join(f"{x:.0f}" for x in dat["P2_MID"])
          + "  OUT " + " ".join(f"{x:.0f}" for x in dat["P2_OUT"])
          + "  | w68 " + " ".join(f"{x:.0f}" for x in dat["width"]))
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for det, mfc in (("P2_MID", "k"), ("P2_OUT", "none")):
        axs[0].plot(E, dat[det] - dat[det][ref], "o", color="k", mfc=mfc, label=f"SPS {det}")
    axs[1].plot(E, dat["width"], "o", color="k", label="SPS (MID−OUT)/√2")
    for (name, k), m in res.items():
        lab = f"sim v_d × {k}" if name == "nominal" else name
        lw = 2.2 if (name == "nominal" and k == 1.0) else 0.9
        ls = "-" if name == "nominal" else "--"
        axs[0].plot(E, m[:, 1] - m[ref, 1], ls, lw=lw, label=lab)
        axs[1].plot(E, m[:, 2], ls, lw=lw, label=lab)
    axs[0].set_ylabel("median hit time − (850 V) [ns]")
    axs[1].set_ylabel("single-chamber 68 % half-width [ns]")
    axs[1].set_yscale("log")
    for ax in axs:
        ax.set_xlabel("drift field [V/cm]")
    axs[1].legend(fontsize=7, ncol=2)
    axs[0].legend(fontsize=7)
    fig.suptitle("SPS drift_mesh_scan_1, Ar/CO2/iso 93/5/2: hit time vs Magboltz v_d (toy DREAM model)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fig_drift_time_vs_field.png", dpi=130)


if __name__ == "__main__":
    main()
