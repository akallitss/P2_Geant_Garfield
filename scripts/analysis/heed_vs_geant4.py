#!/usr/bin/env python3
"""Primary-electron spectra in the drift gas: Garfield/Heed vs Geant4
(plan item 2, 2026-10-14).

Heed: scripts/garfield/heed_mip.cc, a bare gas slab (no window or cathode
upstream). Geant4: Stage A rnd2026 files, EventTree.nPrimDrift (events with no
drift ionization are not stored: added back as zeros from RunMeta.thrown).

    python3 scripts/analysis/heed_vs_geant4.py <heed dir> <rnd2026 dir> <outdir>

Writes heed_vs_geant4.csv (mean, P(n<=k), quantiles per case) and
fig_heed_vs_geant4.png (ArIso / NeIso spectra at 0 deg, both gaps).
"""
import csv
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

GASES = ["ArIso", "ArIso9010", "NeIso", "NeIso8515", "ArCO2Iso9352", "ArCF4Iso"]
LABEL = {"ArIso": "Ar/iso 95/5", "ArIso9010": "Ar/iso 90/10", "NeIso": "Ne/iso 95/5",
         "NeIso8515": "Ne/iso 85/15", "ArCO2Iso9352": "Ar/CO2/iso 93/5/2",
         "ArCF4Iso": "Ar/CF4/iso 88/10/2"}


def heed(d, gas, gap, th, mu=False):
    tag = f"heed_{gas}_{gap}mm_th{th}" + ("_muon150GeV" if mu else "") + ".txt"
    p = d / tag
    if not p.exists():
        return None
    return np.loadtxt(p, usecols=1).astype(int)


def g4(d, gas, gapname, th, part="electron100MeV"):
    p = d / f"p2_{gas}_{gapname}_{part}_th{th}_t0.root"
    if not p.exists():
        return None
    with uproot.open(p) as f:
        thrown = int(f["RunMeta"]["thrown"].array(library="np").sum())
        n = f["EventTree"]["nPrimDrift"].array(library="np").astype(int)
    return np.concatenate([n, np.zeros(max(thrown - len(n), 0), int)])


def summ(n):
    mpv = int(np.bincount(n).argmax())
    return {"mean": round(float(n.mean()), 2), "median": float(np.median(n)), "mpv": mpv,
            "P0": round(float((n == 0).mean()), 5), "Ple5": round(float((n <= 5).mean()), 5),
            "Ple10": round(float((n <= 10).mean()), 5), "q05": float(np.quantile(n, 0.05)),
            "q95": float(np.quantile(n, 0.95)), "N": len(n)}


def main():
    hd, gd, out = (pathlib.Path(a) for a in sys.argv[1:4])
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    cases = [(g, gap, gname, th, False) for g in GASES[:4]
             for gap, gname in [("3.964", "3p964mm"), ("5", "5mm")] for th in (0, 20, 40)]
    cases += [(g, "3.964", "3p964mm", 0, True) for g in ("ArCO2Iso9352", "ArCF4Iso", "ArIso")]
    data = {}
    for g, gap, gname, th, mu in cases:
        h = heed(hd, g, gap, th, mu)
        part = "muon150GeV" if mu else "electron100MeV"
        e = g4(gd, g, gname, th, part)
        if h is None or e is None:
            continue
        data[(g, gname, th, mu)] = (h, e)
        sh, se = summ(h), summ(e)
        rows.append({"gas": g, "gap": gname, "theta": th, "particle": part,
                     **{f"heed_{k}": v for k, v in sh.items()},
                     **{f"g4_{k}": v for k, v in se.items()},
                     "mean_ratio_g4_over_heed": round(se["mean"] / sh["mean"], 3)})
        print(f"{g:13s} {gname:8s} th{th:2d} {part:15s} mean {sh['mean']:6.2f}/{se['mean']:6.2f} "
              f"median {sh['median']:4.0f}/{se['median']:4.0f} mpv {sh['mpv']:3d}/{se['mpv']:3d} "
              f"| P(n<=5) heed {sh['Ple5']:.4f} g4 {se['Ple5']:.4f} | P0 heed {sh['P0']:.4f} g4 {se['P0']:.4f}")
    with open(out / "heed_vs_geant4.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, (g, gname) in zip(axes.flat, [("ArIso", "3p964mm"), ("ArIso", "5mm"),
                                          ("NeIso", "3p964mm"), ("NeIso8515", "3p964mm")]):
        if (g, gname, 0, False) not in data:
            continue
        h, e = data[(g, gname, 0, False)]
        bins = np.arange(-0.5, np.quantile(np.concatenate([h, e]), 0.995) + 1, 1)
        ax.hist(h, bins, density=True, histtype="step", lw=1.6, color="#D62728", label=f"Heed (mean {h.mean():.1f})")
        ax.hist(e, bins, density=True, histtype="step", lw=1.6, color="#2F6DB5", label=f"Geant4 (mean {e.mean():.1f})")
        ax.set_yscale("log")
        ax.set_xlabel("primary electrons in the drift gas")
        ax.set_title(f"{LABEL[g]}, {gname.replace('p', '.')}, 100 MeV e⁻, 0°", fontsize=10)
        ax.legend(fontsize=8)
    fig.suptitle("Primary ionization per crossing: Garfield/Heed (bare gas slab) vs Geant4 (full stack)",
                 fontsize=11)
    fig.savefig(out / "fig_heed_vs_geant4.png", dpi=160)
    print("wrote", out)


if __name__ == "__main__":
    main()
