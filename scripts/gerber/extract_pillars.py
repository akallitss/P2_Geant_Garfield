#!/usr/bin/env python3
"""
The P2 mesh-support pillars, read out of the bulk masks.

Two pillar patterns exist and the detectors are split between them:

    saclay  Insulation_masks/V1/P2_BASKET-Mask_M2_V1.gbr  det1-det4 (Saclay bulk)
            12 477 x Ø0.8 mm on a ~4 mm grid + 5 x Ø6.15 mm
    cern    Bulk_CERN/P2_Mask2.gbr                        det5 (CERN bulk)
            41 366 x Ø0.5 mm on a 2 mm grid  + 5 x Ø6.15 mm

The gerbers are parsed with `load_pillars` from P2_basket_analysis
(cosmic_bench_analysis/p2_mapping.py) -- the same function the cosmic-bench
efficiency maps use -- so the simulation and the analysis share ONE
definition of where the pillars are.  That repo must sit next to this one.
Both masks share the pad-map frame (wedge apex at the origin): the five big
pillars land on identical coordinates in both.

Checks before writing: every pillar lies inside the board-side frame opening
(the AmpGas outline) with its full radius, and no two pillars overlap -- the
Geant4 placements skip the per-volume overlap check, so these ARE the check.

Usage:
    python scripts/gerber/extract_pillars.py            # report + checks
    python scripts/gerber/extract_pillars.py --write    # + include/P2Pillars.hh
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
ROOT = os.path.dirname(REPO)
sys.path.insert(0, os.path.join(REPO, "scripts", "model"))
sys.path.insert(0, os.path.join(ROOT, "P2_basket_analysis", "cosmic_bench_analysis"))

from p2_model import wedge_outline, DRIFT_OPENING, ACTIVE   # noqa: E402

try:
    from p2_mapping import load_pillars                      # noqa: E402
except ImportError as e:
    sys.exit(f"needs P2_basket_analysis next to this repo ({ROOT}): {e}")

MASKS = {
    "saclay": os.path.join(REPO, "design", "gerbers", "bulk_masks_CERN",
                           "older_revisions", "P2_BASKET-Mask_M2_V1.gbr"),
    "cern":   os.path.join(REPO, "design", "gerbers", "bulk_masks_CERN",
                           "P2_Mask2.gbr"),
}
HEADER = os.path.join(REPO, "include", "P2Pillars.hh")


def edge_distance(pts, poly):
    """Distance from each point to the nearest edge of a closed polygon."""
    a = np.asarray(poly)
    b = np.roll(a, -1, axis=0)
    ab = b - a
    t = ((pts[:, None, :] - a[None]) * ab[None]).sum(-1) / (ab ** 2).sum(-1)
    t = np.clip(t, 0, 1)
    proj = a[None] + t[..., None] * ab[None]
    return np.linalg.norm(pts[:, None, :] - proj, axis=-1).min(axis=1)


def drop_covered(name, df):
    """The Saclay mask also draws small pillars under the big ones (6 of them,
    each wholly inside a Ø6.15 mm disc). They add no material; drop them so
    the placements do not overlap. A PARTIAL overlap would be a real conflict."""
    small, big = df[~df.big], df[df.big]
    d = np.hypot(small.x.to_numpy()[:, None] - big.x.to_numpy()[None],
                 small.y.to_numpy()[:, None] - big.y.to_numpy()[None])
    touch = d < small.r.to_numpy()[:, None] + big.r.to_numpy()[None]
    inside = d + small.r.to_numpy()[:, None] <= big.r.to_numpy()[None]
    if (touch & ~inside).any():
        sys.exit(f"{name}: a small pillar partially overlaps a big one")
    covered = inside.any(axis=1)
    if covered.any():
        print(f"{name:7s} dropped {covered.sum()} small pillars lying inside big ones")
    return df.drop(small.index[covered])


def check(name, df):
    from matplotlib.path import Path
    df = drop_covered(name, df)
    pts = df[["x", "y"]].to_numpy()
    r = df["r"].to_numpy()
    opening = np.array(wedge_outline(**DRIFT_OPENING))
    inside = Path(opening).contains_points(pts)
    margin = edge_distance(pts, opening) - r
    tree = cKDTree(pts)
    pairs = tree.query_pairs(2 * r.max() + 0.01, output_type="ndarray")
    gap = (np.linalg.norm(pts[pairs[:, 0]] - pts[pairs[:, 1]], axis=1)
           - r[pairs[:, 0]] - r[pairs[:, 1]]) if len(pairs) else np.array([np.inf])
    small, big = df[~df.big], df[df.big]
    act = 0.5 * np.radians(59.0 - 0.97) * (ACTIVE["r_out"] ** 2 - ACTIVE["r_in"] ** 2)
    dead = (np.pi * df.r ** 2).sum() / act
    print(f"{name:7s} {len(small):6d} small Ø{2*small.r.iloc[0]:.3f} mm + {len(big)} big"
          f" Ø{2*big.r.iloc[0]:.3f} mm | dead area {100*dead:.2f} % of the active area"
          f" | min edge margin {margin.min():.3f} mm | min pillar gap {gap.min():.3f} mm")
    ok = inside.all() and margin.min() > 0 and gap.min() > 0
    if not ok:
        sys.exit(f"{name}: pillar outside the AmpGas outline or overlapping -- not writing")
    return small, big


def write(sets):
    lines = [
        "#pragma once",
        "// P2Pillars.hh  --  GENERATED, DO NOT EDIT BY HAND",
        "//",
        "// Regenerate with:",
        "//     python scripts/gerber/extract_pillars.py --write",
        "//",
        "// Mesh-support pillar centres (mm, gerber frame: wedge apex at the",
        "// origin) read with P2_basket_analysis' load_pillars, so simulation and",
        "// cosmic-bench analysis use one pillar map.  Checked: every pillar sits",
        "// inside the AmpGas outline with its full radius, none overlap.",
        "",
        "namespace P2 {",
        "",
        "struct PillarXY { double x, y; };",
        "",
    ]
    for name, (small, big) in sets.items():
        src = os.path.relpath(MASKS[name], REPO)
        tag = name.capitalize()
        lines += [
            f"// {name}: {src}",
            f"constexpr double kPillarR{tag}    = {small.r.iloc[0]:.4f};   // small pillar radius",
            f"constexpr double kPillarRBig{tag} = {big.r.iloc[0]:.4f};   // big pillar radius",
            f"constexpr int    kNPillars{tag}    = {len(small)};",
            f"constexpr int    kNPillarsBig{tag} = {len(big)};",
            f"inline constexpr PillarXY kPillarsBig{tag}[] = {{",
        ]
        lines += [f"    {{{x:.4f}, {y:.4f}}}," for x, y in big[["x", "y"]].to_numpy()]
        lines += ["};", f"inline constexpr PillarXY kPillars{tag}[] = {{"]
        xy = small[["x", "y"]].to_numpy()
        for i in range(0, len(xy), 4):
            lines.append("    " + " ".join(f"{{{x:.4f},{y:.4f}}}," for x, y in xy[i:i + 4]))
        lines += ["};", ""]
    lines += ["}  // namespace P2", ""]
    with open(HEADER, "w") as fh:
        fh.write("\n".join(lines))
    print("wrote", HEADER, f"({os.path.getsize(HEADER)/1e6:.2f} MB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true", help="write include/P2Pillars.hh")
    args = ap.parse_args()
    sets = {}
    for name, path in MASKS.items():
        df = load_pillars(path)
        if df.empty:
            sys.exit(f"{name}: no pillars read from {path}")
        sets[name] = check(name, df)
    if args.write:
        write(sets)


if __name__ == "__main__":
    main()
