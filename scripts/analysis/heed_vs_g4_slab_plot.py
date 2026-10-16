#!/usr/bin/env python3
"""Figures for the like-for-like Heed vs Geant4 check (2026-10-16).

    python3 scripts/analysis/heed_vs_g4_slab_plot.py <work dir> <outdir>

fig_slab_vs_heed.png : Ar/iso 95/5, 3.964 mm, 100 MeV e-, bare slab: Heed vs
                       Geant4 PAI / option4 (PIXE off) and option4 + PIXE +
                       100 um step (the production settings).
fig_stack_vs_slab.png: the same physics (PAI, PIXE off) on the bare slab and
                       in the full P2 stack.
"""
import glob
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 13})


def g4(work, tag):
    out = []
    for fn in sorted(glob.glob(str(work / "out" / f"{tag}_s*_t0.root"))):
        with uproot.open(fn) as f:
            n = f["EventTree"]["nPrimDrift"].array(library="np")
            th = int(f["RunMeta"]["thrown"].array(library="np").sum())
        out.append(np.concatenate([n, np.zeros(th - len(n), int)]))
    return np.concatenate(out)


def step(ax, n, label, **kw):
    b = np.arange(-0.5, 151.5, 2)
    h, _ = np.histogram(n, b)
    ax.step(b[:-1], h / len(n), where="post", label=f"{label}: {n.mean():.1f} / {np.median(n):.0f}", **kw)


def main():
    work, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    heed = np.loadtxt(work / "heed" / "heed_ArIso_4mm.txt", usecols=1)
    fig, ax = plt.subplots(figsize=(8, 5))
    step(ax, heed, "Garfield/Heed", color="k", lw=2.2)
    step(ax, g4(work, "np_slab_ArIso_4mm_pai_st0"), "Geant4 PAI", color="#1f6fb4", lw=1.6)
    step(ax, g4(work, "np_slab_ArIso_4mm_opt4_st0"), "Geant4 option4 (Urban)", color="#e07b00", lw=1.4)
    step(ax, g4(work, "slab_ArIso_4mm_opt4_st100"), "Geant4 option4 + PIXE + 100 µm step", color="#888888",
         lw=1.4, ls="--")
    ax.set_xlim(0, 140)
    ax.set_xlabel("primary electrons in the drift gas")
    ax.set_ylabel("fraction of events / 2")
    ax.set_title("Bare gas slab: Ar/iso 95/5, 3.964 mm, 100 MeV e-", fontsize=13)
    ax.legend(title="mean / median", fontsize=11, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_slab_vs_heed.png", dpi=150)

    fig, ax = plt.subplots(figsize=(8, 5))
    step(ax, heed, "Garfield/Heed, bare slab", color="k", lw=2.2)
    step(ax, g4(work, "np_slab_ArIso_4mm_pai_st0"), "Geant4 PAI, bare slab", color="#1f6fb4", lw=1.6)
    step(ax, g4(work, "np_p2_ArIso_pai_st0"), "Geant4 PAI, full P2 stack", color="#c0392b", lw=1.6)
    ax.set_xlim(0, 140)
    ax.set_xlabel("primary electrons in the drift gas")
    ax.set_ylabel("fraction of events / 2")
    ax.set_title("Same physics (PAI, PIXE off), with and without the material", fontsize=13)
    ax.legend(title="mean / median", fontsize=11, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_stack_vs_slab.png", dpi=150)


if __name__ == "__main__":
    main()
