#!/usr/bin/env python3
"""
rnd2026 first look: photons vs 100 MeV electrons (Geant4 stage only).

For each gas x drift gap of the `--scan rnd2026` batch:
  * P(>= 1 primary in the drift gas) per photon energy, split into
    conversions in the gas and electrons coming out of a wall;
  * drift-gas deposit spectra, photons (detected events) vs electrons;
  * an upper cut on the drift-gas deposit set at 99 % / 95 % electron
    efficiency, and the fraction of detected photons it removes;
  * the detected photon rate per Micromegas layer, folding P(E) with
    Matthieu's P2Sim flux (BkgHistograms.root, hEnergy_i, MHz).

Observable: edepDrift (eV) -- docs/OUTPUT_FORMAT.md §4 says to prefer edep
over nPrimary until Stage B; "detected" = nPrimDrift >= 1. No gain, no
electronics: this is the primary-ionization layer of the answer only.

Usage (lxplus, LCG view python):
    python3 rnd2026_photons_electrons.py --indir /eos/.../rnd2026 \
        --bkg BkgHistograms.root --out rnd2026_out
"""
import argparse
import csv
import glob
import json
import os
import re

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GASES = ["ArIso", "ArIso9010", "NeIso", "NeIso8515"]
LABEL = {"ArIso": "Ar/iso 95/5", "ArIso9010": "Ar/iso 90/10",
         "NeIso": "Ne/iso 95/5", "NeIso8515": "Ne/iso 85/15"}
COLOR = {"ArIso": "#1f5fa8", "ArIso9010": "#6aa7e0",
         "NeIso": "#c4501a", "NeIso8515": "#f0a066"}
GAPS = ["3p964mm", "5mm"]
GAPLABEL = {"3p964mm": "3.96 mm (4 mm frame)", "5mm": "5 mm"}
# Gas vs wall is split on the conversion VOLUME, not interactionClass: up
# to commit f643863 the class labelled a soft (<20 keV) beam photon absorbed
# in the gas as fluorescence, because the beam is born in the World volume.
GAS_VOLUMES = ("DriftGas", "AmpGas")
BRANCHES = ["edepDrift", "nPrimDrift", "interactionClass", "convVolume"]

TAG = re.compile(r"p2_(?P<gas>[A-Za-z0-9]+)_(?P<gap>[0-9p]+mm)_(?P<part>gamma|electron|muon)"
                 r"(?P<e>[0-9p]+)(?P<unit>keV|MeV|GeV)_th(?P<th>[0-9p]+)")


def load(path):
    f = uproot.open(path)
    thrown = int(f["RunMeta"]["thrown"].array(library="np").sum())
    ev = f["EventTree"].arrays(BRANCHES, library="np")
    return thrown, ev


def scan(indir):
    runs = {}
    for p in sorted(glob.glob(os.path.join(indir, "p2_*.root"))):
        m = TAG.search(os.path.basename(p))
        if not m:
            continue
        e = float(m["e"].replace("p", "."))
        e_keV = e * {"keV": 1, "MeV": 1e3, "GeV": 1e6}[m["unit"]]
        key = (m["gas"], m["gap"], m["part"], e_keV, float(m["th"].replace("p", ".")))
        runs[key] = p
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--bkg", required=True)
    ap.add_argument("--out", default="rnd2026_out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    runs = scan(a.indir)

    # ── electrons: deposit spectra and the upper cuts ──────────────────────
    elec = {}
    for (gas, gap, part, e, th), p in runs.items():
        if part != "electron":
            continue
        thrown, ev = load(p)
        ed = ev["edepDrift"][ev["nPrimDrift"] >= 1] / 1e3      # keV
        elec[(gas, gap, th)] = dict(thrown=thrown, edep=ed,
                                    eff=len(ed) / thrown,
                                    median=float(np.median(ed)),
                                    cut99=float(np.percentile(ed, 99)),
                                    cut95=float(np.percentile(ed, 95)))
        h, b = np.histogram(ed, bins=np.linspace(0, 10, 201))
        elec[(gas, gap, th)]["mpv"] = float(0.5 * (b[h.argmax()] + b[h.argmax() + 1]))

    # ── photons ────────────────────────────────────────────────────────────
    rows, phot = [], {}
    for (gas, gap, part, e, th), p in sorted(runs.items()):
        if part != "gamma":
            continue
        thrown, ev = load(p)
        det = ev["nPrimDrift"] >= 1
        ed = ev["edepDrift"][det] / 1e3
        e0 = elec[(gas, gap, 0.0)]
        r = dict(gas=gas, gap=gap, E_keV=e, thrown=thrown,
                 P_detect=det.sum() / thrown,
                 P_gas=(det & np.isin(ev["convVolume"], GAS_VOLUMES)).sum() / thrown,
                 P_wall=(det & ~np.isin(ev["convVolume"], GAS_VOLUMES)).sum() / thrown,
                 n_detected=int(det.sum()),
                 median_edep_keV=float(np.median(ed)) if len(ed) else np.nan,
                 frac_above_cut99=float((ed > e0["cut99"]).mean()) if len(ed) else np.nan,
                 frac_above_cut95=float((ed > e0["cut95"]).mean()) if len(ed) else np.nan)
        rows.append(r)
        phot[(gas, gap, e)] = ed
    with open(os.path.join(a.out, "photon_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    def series(gas, gap, key):
        rr = sorted((r for r in rows if r["gas"] == gas and r["gap"] == gap),
                    key=lambda r: r["E_keV"])
        return np.array([r["E_keV"] for r in rr]), np.array([r[key] for r in rr]), rr

    # ── figure 1: P(>=1 primary) vs E ──────────────────────────────────────
    fig, axs = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, gap in zip(axs, GAPS):
        for gas in GASES:
            E, P, rr = series(gas, gap, "P_detect")
            err = np.sqrt(np.array([r["n_detected"] for r in rr])) / np.array([r["thrown"] for r in rr])
            ax.errorbar(E, P, yerr=err, marker="o", color=COLOR[gas], label=LABEL[gas])
            ax.plot(E, series(gas, gap, "P_gas")[1], ls=":", color=COLOR[gas])
        ax.axvspan(40, 100, color="0.9", zorder=0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("photon energy [keV]"); ax.set_title(f"drift gas {GAPLABEL[gap]}")
        ax.grid(alpha=0.3, which="both")
    axs[0].set_ylabel("P(≥ 1 primary in drift gas) per photon")
    axs[0].legend(title="solid: all · dotted: gas conversions only")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig1_P_detect.png"), dpi=130)

    # ── figure 2: deposit spectra, photons vs electrons ───────────────────
    fig, axs = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    bins = np.logspace(-1, np.log10(200), 45)
    for ax, gas in zip(axs.flat, GASES):
        gap = "3p964mm"
        ed = elec[(gas, gap, 0.0)]["edep"]
        ax.hist(ed, bins=bins, weights=np.full(len(ed), 1 / len(ed)), histtype="step", lw=2.5, color="k",
                label="e⁻ 100 MeV, 0°")
        for e, c in zip([20, 40, 60, 100], ["#2b8cbe", "#41ab5d", "#e6550d", "#756bb1"]):
            x = phot.get((gas, gap, float(e)))
            if x is not None and len(x):
                ax.hist(x, bins=bins, weights=np.full(len(x), 1 / len(x)), histtype="step", lw=1.5, color=c,
                        label=f"γ {e} keV (detected)")
        ax.axvline(elec[(gas, gap, 0.0)]["cut99"], color="k", ls="--", lw=1)
        ax.set_xscale("log"); ax.set_title(f"{LABEL[gas]}, {GAPLABEL[gap]}")
        ax.grid(alpha=0.3)
    for ax in axs[1]:
        ax.set_xlabel("energy deposited in the drift gas [keV]")
    for ax in axs[:, 0]:
        ax.set_ylabel("fraction of events per bin")
    axs[0, 0].legend(fontsize=9, title="dashed: 99 % e⁻ cut")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig2_edep_spectra.png"), dpi=130)

    # ── figure 3: photons removed by the upper cut ─────────────────────────
    fig, axs = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, gap in zip(axs, GAPS):
        for gas in GASES:
            E, F, _ = series(gas, gap, "frac_above_cut99")
            ax.plot(E, F, marker="o", color=COLOR[gas], label=LABEL[gas])
            ax.plot(E, series(gas, gap, "frac_above_cut95")[1], ls=":", color=COLOR[gas])
        ax.set_xscale("log"); ax.set_ylim(0, 1)
        ax.set_xlabel("photon energy [keV]"); ax.set_title(f"drift gas {GAPLABEL[gap]}")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("detected photons above the upper cut")
    axs[0].legend(title="solid: cut at 99 % e⁻ eff. · dotted: 95 %")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig3_upper_cut.png"), dpi=130)

    # ── fold with the P2Sim flux ───────────────────────────────────────────
    bk = uproot.open(a.bkg)
    fold = {}
    for i in range(3):
        h = bk[f"hEnergy_{i}"]
        c_keV = h.axis().centers() * 1e3
        flux = h.values()                       # MHz per 5 keV bin
        sel = (c_keV >= 10) & (c_keV <= 150)
        for gas in GASES:
            for gap in GAPS:
                E, P, rr = series(gas, gap, "P_detect")
                F99 = series(gas, gap, "frac_above_cut99")[1]
                Pi = np.exp(np.interp(np.log(c_keV[sel]), np.log(E), np.log(P)))
                Fi = np.interp(np.log(c_keV[sel]), np.log(E), F99)
                det = float((flux[sel] * Pi).sum())
                fold[(i, gas, gap)] = dict(
                    flux_covered_MHz=float(flux[sel].sum()),
                    flux_total_MHz=float(flux.sum()),
                    detected_MHz=det,
                    detected_after_cut99_MHz=float((flux[sel] * Pi * (1 - Fi)).sum()))
    with open(os.path.join(a.out, "fold_p2sim.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["layer", "gas", "gap", "flux_10_150keV_MHz", "flux_total_MHz",
                    "detected_MHz", "detected_after_cut99_MHz"])
        for (i, gas, gap), v in sorted(fold.items()):
            w.writerow([f"MM{i}", gas, gap, f"{v['flux_covered_MHz']:.0f}",
                        f"{v['flux_total_MHz']:.0f}", f"{v['detected_MHz']:.2f}",
                        f"{v['detected_after_cut99_MHz']:.2f}"])

    fig, ax = plt.subplots(figsize=(10, 5))
    xs = np.arange(len(GASES))
    for k, gap in enumerate(GAPS):
        vals = [fold[(0, g, gap)]["detected_MHz"] for g in GASES]
        ax.bar(xs + (k - 0.5) * 0.38, vals, width=0.36,
               color=[COLOR[g] for g in GASES], alpha=1.0 if k == 0 else 0.55,
               label=GAPLABEL[gap])
    ax.set_xticks(xs); ax.set_xticklabels([LABEL[g] for g in GASES])
    ax.set_ylabel("detected photon rate, MM0 [MHz, whole layer]")
    ax.set_title("P2Sim flux (10–150 keV) × P(≥ 1 primary)  ·  solid 3.96 mm, light 5 mm")
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig4_fold_MM0.png"), dpi=130)

    # ── electrons table ────────────────────────────────────────────────────
    erows = [dict(gas=g, gap=gp, theta=th, eff=v["eff"], mpv_keV=v["mpv"],
                  median_keV=v["median"], cut99_keV=v["cut99"], cut95_keV=v["cut95"])
             for (g, gp, th), v in sorted(elec.items())]
    with open(os.path.join(a.out, "electron_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(erows[0]))
        w.writeheader(); w.writerows(erows)
    print(json.dumps({"photon_runs": len(rows), "electron_runs": len(erows)}))


if __name__ == "__main__":
    main()
