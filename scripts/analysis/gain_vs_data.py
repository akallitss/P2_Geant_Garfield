#!/usr/bin/env python3
"""Simulated avalanche gain (geometry DAG, mesh_gain.cc) vs the measured gain
(2026-10-14).

Simulation: docs/results/gain/geom_gain_lines.txt (Townsend alpha integrated
along drift lines through the COMSOL cell, Penning r = 0-0.7).
Data (docs/results/gain/dream_calibration.md): DREAM 200 fC range, 18 ADC/fC;
the 180 ns peak sees 0.46-0.76 of the ion charge (central 0.6).
  SPS P2_OUT beam_nominal_meshscan_1, Ar/CO2/iso 93/5/2: Landau MPV of the
    cluster charge / (22 +- 5 primaries x T 0.975)
  Cosmic det2 Fe55 mesh scan, Ar/iso 95/5: photopeak / 227 primaries

    python3 scripts/analysis/gain_vs_data.py
"""
import csv
import pathlib
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = REPO / "docs" / "results" / "gain"
E_PER_ADC = 6241.5 / 18.0
COLL = (0.46, 0.6, 0.76)   # ion-charge fraction at the 180 ns peak: low, central, high
SPS = pathlib.Path("/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/mpgd26_workspace/"
                   "products/analysis/P2_OUT/beam_nominal_meshscan_1/scan/20_beam_spectra/"
                   "beam_mpv_vs_hv_spark_vetoed.csv")
FE55 = pathlib.Path("/local/home/ak271430/Documents/PostDocSaclay/data/Cosmic_Bench/Analysis/det2/"
                    "p2_fe55_det2_det3_mesh_scan_7-18-26/fe55_scan/18_fe55_spectra/fe55_gain_vs_hv_spark_vetoed.csv")


def sim():
    """{(tag, gas): {r: G}}"""
    d = {}
    for line in open(OUT / "geom_gain_lines.txt"):
        m = re.match(r"gain_(\S+)_(ArCO2Iso9352|ArIso)\.out:GAIN .* r=([\d.]+) .*mean_gain=(\S+)", line)
        if m:
            d.setdefault((m[1], m[2]), {})[float(m[3])] = float(m[4])
    return d


def tag_of(gap, vm):
    return f"g{gap}" + ("" if vm == 450 else f"_vm{vm}")


def data(path, vcol, qcol, n0):
    V, Q = [], []
    for r in csv.DictReader(open(path)):
        if r[qcol] and (path != FE55 or float(r[vcol]) >= 385):
            V.append(float(r[vcol]))
            Q.append(float(r[qcol]) * E_PER_ADC / n0)
    return np.array(V), np.array(Q)


def main():
    s = sim()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6), constrained_layout=True)
    col = {120: "#D62728", 135: "#FF7F0E", 150: "#1F77B4", 165: "#2CA02C"}

    # (a) Ar/CO2/iso: G vs mesh V, gaps at r = 0.3 / 0.4, SPS data
    ax = axes[0]
    for gap, c in col.items():
        for r, ls in ((0.3, "--"), (0.4, "-")):
            V = [430, 450, 470]
            G = [s[(tag_of(gap, v), "ArCO2Iso9352")][r] for v in V]
            ax.plot(V, G, ls, color=c, marker="o", ms=4,
                    label=f"{gap} µm, r = {r}" if True else None)
    V, Q = data(SPS, "hv", "mpv_adc", 22 * 0.975)
    ax.errorbar(V, Q / COLL[1], yerr=[Q / COLL[1] - Q / COLL[2] * 22 / 27, Q / COLL[0] * 22 / 17 - Q / COLL[1]],
                fmt="ks", ms=7, label="SPS P2_OUT (MIP MPV), band: primaries 17-27, ion fraction 0.46-0.76")
    ax.set_yscale("log")
    ax.set_xlabel("mesh voltage [V]")
    ax.set_ylabel("gain")
    ax.set_title("Ar/CO2/iso 93/5/2: simulated gain vs SPS", fontsize=11)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3, which="both")

    # (b) r scan at 450 V for each gap vs the SPS band at 450 V
    ax = axes[1]
    q450 = Q[V == 450].mean()
    lo, hi = q450 / COLL[2] * 22 / 27, q450 / COLL[0] * 22 / 17
    ax.axhspan(lo, hi, color="#BBBBBB", alpha=0.5, label="SPS 450 V (range)")
    rows = []
    for gap, c in col.items():
        d = s[(tag_of(gap, 450), "ArCO2Iso9352")]
        r = np.array(sorted(d))
        G = np.array([d[x] for x in r])
        ax.plot(r, G, "-o", color=c, ms=4, label=f"{gap} µm")
        lr = lambda g: float(np.interp(np.log(g), np.log(G), r))  # noqa: E731
        rows.append(("ArCO2Iso9352 450 V", gap, lr(lo), lr(np.sqrt(lo * hi)), lr(hi)))
    ax.set_yscale("log")
    ax.set_ylim(10, 3e5)
    ax.set_xlabel("Penning transfer rate r")
    ax.set_ylabel("gain at 450 V")
    ax.set_title("Ar/CO2/iso: which r matches the SPS gain, per gap", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    # (c) Ar/iso 95/5: Fe55 (det2) vs the 150 um cell at 415 / 450 V
    ax = axes[2]
    V, Q = data(FE55, "hv", "peak_adc", 227)
    ax.errorbar(V, Q / COLL[1], yerr=[Q / COLL[1] - Q / COLL[2], Q / COLL[0] - Q / COLL[1]],
                fmt="ks", ms=7, label="cosmic det2 Fe55 (band: ion fraction)")
    for r, c in ((0.2, "#9467BD"), (0.3, "#1F77B4"), (0.4, "#D62728"), (0.5, "#8C564B")):
        Vs = [415, 430, 450, 470]
        Gs = [s[(t, "ArIso")][r] for t in ("vm415", "g150_vm430", "g150", "g150_vm470")]
        ax.plot(Vs, Gs, "-o", color=c, ms=4, label=f"sim 150 µm, r = {r}")
    d = s[("vm415", "ArIso")]
    r = np.array(sorted(d))
    G = np.array([d[x] for x in r])
    q415 = Q[V == 415].mean()
    for lab, g in (("lo", q415 / COLL[2]), ("mid", q415 / COLL[1]), ("hi", q415 / COLL[0])):
        pass
    rows.append(("ArIso 415 V (Fe55)", 150, *[float(np.interp(np.log(q415 / f), np.log(G), r)) for f in COLL[::-1]]))
    ax.set_yscale("log")
    ax.set_xlabel("mesh voltage [V]")
    ax.set_ylabel("gain")
    ax.set_title("Ar/iso 95/5: simulated gain vs Fe55 (det2)", fontsize=11)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, which="both")
    fig.savefig(OUT / "fig_gain_vs_data.png", dpi=150)

    with open(OUT / "penning_r_from_gain.csv", "w") as f:
        f.write("case,gap_um,r_low,r_central,r_high\n")
        for row in rows:
            f.write(f"{row[0]},{row[1]},{row[2]:.3f},{row[3]:.3f},{row[4]:.3f}\n")
            print(row[0], row[1], f"r = {row[3]:.2f} ({row[2]:.2f}-{row[4]:.2f})")
    # simulated e-folding 430-470 V at r = 0.3 / 0.4
    for gap in col:
        for r in (0.3, 0.4):
            g = [s[(tag_of(gap, v), "ArCO2Iso9352")][r] for v in (430, 470)]
            print(f"e-folding {gap} um r={r}: {40 / np.log(g[1] / g[0]):.1f} V")


if __name__ == "__main__":
    main()
