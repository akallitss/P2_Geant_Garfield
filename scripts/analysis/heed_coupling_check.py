#!/usr/bin/env python3
"""Validation of the Geant4 + Garfield++/Heed coupling (mm_sim --heed,
src/HeedGasModel.cc), 2026-10-16.

1. Bare slab: Geant4+Heed against Heed standalone (heed_mip.cc), same gas and
   particle -> must agree in shape, not only in mean (fraction of events per
   range of primaries, P(n<10), KS distance).
2. Full P2 stack, 100 MeV e-: Geant4+Heed vs Geant4 PAI (both PIXE off), and
   the layers downstream of the gas (mesh, pad copper): the beam particle
   must leave the gas as before.
3. 30 keV photons in the full stack: below the Heed thresholds nothing may
   change (P(drift ionization), nPrimDrift).

    python3 scripts/analysis/heed_coupling_check.py <work dir> <outdir>
"""
import glob
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from scipy.stats import ks_2samp

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 13})
RANGES = [(0, 10), (10, 15), (15, 40), (40, 100), (100, 10**9)]


def g4(work, tag, branches=("nPrimDrift",)):
    cols = {b: [] for b in branches}
    thrown = 0
    for fn in sorted(glob.glob(str(work / "out" / f"{tag}_s*_t0.root"))):
        with uproot.open(fn) as f:
            th = int(f["RunMeta"]["thrown"].array(library="np").sum())
            ev = f["EventTree"].arrays(list(branches), library="np")
        for b in branches:
            pad = np.zeros(th - len(ev[b]))
            cols[b].append(np.concatenate([ev[b], pad]))
        thrown += th
    return {b: np.concatenate(v) for b, v in cols.items()}, thrown


def heed(work, name):
    return np.loadtxt(work / "heed" / f"{name}.txt", usecols=1)


def shape(n):
    return [100 * np.mean((n >= a) & (n < b)) for a, b in RANGES]


def row(label, n, ref=None):
    s = f"{label:34s} {n.mean():7.2f} {np.median(n):6.0f} " + " ".join(f"{x:6.2f}" for x in shape(n))
    if ref is not None:
        s += f"   KS {ks_2samp(n, ref).statistic:.4f}"
    return s


def main():
    work, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    lines = ["case                                 mean median   <10  10-15  15-40 40-100   >100"]

    # 1. bare slab
    for gas, h, tg, tp in [("Ar/iso e-", "heed_ArIso_4mm", "hd_slab_ArIso_4mm_heed_st0", "np_slab_ArIso_4mm_pai_st0"),
                           ("Ne/iso e-", "heed_NeIso_4mm", "hd_slab_NeIso_4mm_heed_st0", "np_slab_NeIso_4mm_pai_st0"),
                           ("Ar/iso mu 150 GeV", "heed_ArIso_4mm_mu", "hd_slab_ArIso_4mm_mu_heed_st0", None)]:
        ref = heed(work, h)
        lines.append(f"-- bare slab, {gas}")
        lines.append(row("Heed standalone", ref))
        lines.append(row("Geant4 + Heed (--heed)", g4(work, tg)[0]["nPrimDrift"], ref))
        if tp:
            lines.append(row("Geant4 PAI (edep/W)", g4(work, tp)[0]["nPrimDrift"], ref))

    # 2. full stack
    br = ("nPrimDrift", "edepDrift", "edepMesh", "edepPadCu", "nPrimAmp")
    lines.append("-- full P2 stack, Ar/iso, 100 MeV e-, PIXE off")
    stack = {}
    for lab, tag in [("Geant4 PAI, no step limit", "np_p2_ArIso_pai_st0"),
                     ("Geant4 + Heed, no step limit", "hd_p2_ArIso_heed_st0"),
                     ("Geant4 + Heed, 100 um step", "hd_p2_ArIso_heed_st100")]:
        d, _ = g4(work, tag, br)
        stack[lab] = d
        lines.append(row(lab, d["nPrimDrift"]))
    lines.append("   downstream means [eV]: mesh, pad Cu, amp-gas primaries")
    for lab, d in stack.items():
        lines.append(f"   {lab:31s} {d['edepMesh'].mean():8.0f} {d['edepPadCu'].mean():8.0f} {d['nPrimAmp'].mean():6.2f}")

    # 3. photons
    lines.append("-- full P2 stack, 30 keV photons: P(drift ionization), mean nPrimDrift when > 0")
    for lab, tag in [("Geant4 PAI", "hd_p2_ArIso_g30_pai"), ("Geant4 + Heed", "hd_p2_ArIso_g30_heed")]:
        d, th = g4(work, tag)
        n = d["nPrimDrift"]
        k = n[n > 0]
        lines.append(f"   {lab:31s} P = {len(k) / th:.5f} +- {np.sqrt(len(k)) / th:.5f}   <n> = {k.mean():.1f}")

    txt = "\n".join(lines)
    print(txt)
    (out / "heed_coupling_check.txt").write_text(txt + "\n")

    fig, ax = plt.subplots(figsize=(8, 5))
    b = np.arange(-0.5, 151.5, 2)
    for lab, n, kw in [("Garfield/Heed standalone", heed(work, "heed_ArIso_4mm"), dict(color="k", lw=2.2)),
                       ("Geant4 + Heed (coupled)", g4(work, "hd_slab_ArIso_4mm_heed_st0")[0]["nPrimDrift"],
                        dict(color="#c0392b", lw=1.6, ls="--")),
                       ("Geant4 PAI, edep/W", g4(work, "np_slab_ArIso_4mm_pai_st0")[0]["nPrimDrift"],
                        dict(color="#1f6fb4", lw=1.4))]:
        h, _ = np.histogram(n, b)
        ax.step(b[:-1], h / len(n), where="post", label=f"{lab}: P(n<10) {100 * np.mean(n < 10):.1f} %", **kw)
    ax.set_xlim(0, 140)
    ax.set_xlabel("primary electrons in the drift gas")
    ax.set_ylabel("fraction of events / 2")
    ax.set_title("Bare gas slab: Ar/iso 95/5, 3.964 mm, 100 MeV e-", fontsize=13)
    ax.legend(fontsize=11, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fig_heed_coupling_slab.png", dpi=150)


if __name__ == "__main__":
    main()
