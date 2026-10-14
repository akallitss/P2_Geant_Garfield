#!/usr/bin/env python3
"""Check the Geant4 pillar map (include/P2Pillars.hh) against the insulation
masks in Detector_Drawings (plan item 6, 2026-10-14).

For each mask: number of small / big pillars, radii, nearest-neighbour
spacing, covered area fraction; and, for the masks the simulation uses
(saclay = V1 for det1-4, cern = P2_Mask2 for det5), a one-to-one position
match against the header. Also writes an overlay figure.

    python3 scripts/gerber/check_pillars.py
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
ROOT = os.path.dirname(REPO)
sys.path.insert(0, os.path.join(ROOT, "P2_basket_analysis", "cosmic_bench_analysis"))
from p2_mapping import load_pillars  # noqa: E402

DD = os.path.join(ROOT, "Detector_Drawings")
MASKS = {
    "V1 (det1-4, used)": os.path.join(DD, "Version_Apr26/Insulation_masks/V1/P2_BASKET-Mask_M2_V1.gbr"),
    "V2": os.path.join(DD, "Version_Apr26/Insulation_masks/V2/P2_BASKET-Mask_M2_V2.gbr"),
    "V3B_final": os.path.join(DD, "Version_Sept26/Insulation_masks/V3B_final/P2_BASKET-Mask_M2_V3B_final.gbr"),
    "CERN P2_Mask2 (det5, used)": os.path.join(REPO, "design/gerbers/bulk_masks_CERN/P2_Mask2.gbr"),
}
OUT = os.path.join(REPO, "docs", "results", "pillars")


def header(kind):
    s = open(os.path.join(REPO, "include", "P2Pillars.hh")).read()
    def arr(name):
        m = re.search(rf"{name}\[\] = \{{(.*?)\}};", s, re.S)
        return np.array([[float(a), float(b)] for a, b in
                         re.findall(r"\{([-\d.]+),\s*([-\d.]+)\}", m.group(1))])
    r = float(re.search(rf"kPillarR{kind}\s*=\s*([\d.]+)", s).group(1))
    rb = float(re.search(rf"kPillarRBig{kind}\s*=\s*([\d.]+)", s).group(1))
    return arr(f"kPillars{kind}"), arr(f"kPillarsBig{kind}"), r, rb


def stats(df):
    sm, bg = df[~df.big], df[df.big]
    xy = sm[["x", "y"]].to_numpy()
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    return {"n_small": len(sm), "r_small": np.round(np.unique(sm.r.round(4)), 4).tolist(),
            "n_big": len(bg), "r_big": np.round(np.unique(bg.r.round(4)), 4).tolist(),
            "nn_median": round(float(np.median(nn)), 3),
            "area_mm2": round(float((np.pi * sm.r ** 2).sum() + (np.pi * bg.r ** 2).sum()), 1)}


def match(xy_mask, xy_hdr):
    d, _ = cKDTree(xy_hdr).query(xy_mask)
    d2, _ = cKDTree(xy_mask).query(xy_hdr)
    return {"mask_unmatched(>1um)": int((d > 1e-3).sum()),
            "header_unmatched(>1um)": int((d2 > 1e-3).sum()),
            "max_offset_um": round(1e3 * float(d2.max()), 3)}


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    for name, path in MASKS.items():
        df = load_pillars(path)
        if len(df) == 0:
            print(f"{name}: not read ({path})")
            continue
        res[name] = (df, stats(df))
        print(f"{name:28s} {res[name][1]}")
    for name, kind in [("V1 (det1-4, used)", "Saclay"), ("CERN P2_Mask2 (det5, used)", "Cern")]:
        if name not in res:
            continue
        df = res[name][0]
        hs, hb, r, rb = header(kind)
        sm = df[~df.big][["x", "y"]].to_numpy()
        print(f"  header {kind}: {len(hs)} small r={r} + {len(hb)} big r={rb};",
              "small", match(sm, hs), "big", match(df[df.big][["x", "y"]].to_numpy(), hb))
    # overlay: one region, V1 vs V2 vs V3B_final vs header
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    box = (300, 330, 170, 200)
    hs, hb, r, rb = header("Saclay")
    sel = (hs[:, 0] > box[0]) & (hs[:, 0] < box[1]) & (hs[:, 1] > box[2]) & (hs[:, 1] < box[3])
    for (x, y) in hs[sel]:
        ax.add_patch(plt.Circle((x, y), r, color="#2F6DB5", alpha=0.5))
    for name, c, mk in [("V2", "#D62728", "x"), ("V3B_final", "#2CA02C", "+")]:
        if name in res:
            d = res[name][0]
            s = (d.x > box[0]) & (d.x < box[1]) & (d.y > box[2]) & (d.y < box[3]) & ~d.big
            ax.plot(d.x[s], d.y[s], mk, color=c, ms=5, label=f"{name} mask centres")
    for (x, y) in hb:
        ax.add_patch(plt.Circle((x, y), rb, fill=False, color="#2F6DB5", lw=2))
    ax.plot([], [], "o", color="#2F6DB5", alpha=0.5, label="Geant4 (P2Pillars.hh = V1)")
    ax.set_xlim(box[0], box[1])
    ax.set_ylim(box[2], box[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm] (wedge apex frame)")
    ax.set_ylabel("y [mm]")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("Pillars around a big pillar: Geant4 vs mask versions", fontsize=11)
    fig.savefig(os.path.join(OUT, "pillars_overlay.png"), dpi=160)
    print("wrote", os.path.join(OUT, "pillars_overlay.png"))


if __name__ == "__main__":
    main()
