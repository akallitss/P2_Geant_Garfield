#!/usr/bin/env python3
"""
P2 readout pad map, from the Maarten/Irakli channel-mapping files.

The mapping ships in two flavours under design/mapping/:

  connector_<i>.txt          7 columns, no cartesian coordinates, pads ordered
                             from the OUTER radius inward
      Connector Channel PadIndex AbsIndex Radius Phi DeltaPhi

  Mapping/connector_<i>.txt  9 columns, adds X, Y and PadName, pads ordered
                             from the INNER radius outward
      Connector PinName Channel PadIndex PadName X Y Radius Phi DeltaPhi

*** Phi in both files is in RADIANS, despite sitting next to a Radius in mm. ***
Verified against the X/Y columns: atan2(Y, X) reproduces Phi exactly.

The pads are "equalized": the radial pitch is a constant 11.429 mm while the
number of pads per ring grows with radius (9 at the inner edge, 51 at the
outer), keeping each pad close to 11.4 x 12 mm everywhere.

Usage:
    python scripts/gerber/analyze_p2_readout.py
    python scripts/gerber/analyze_p2_readout.py --plot padmap.pdf
"""

from __future__ import annotations

import argparse
import collections
import glob
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
MAPDIR = os.path.join(REPO, "design", "mapping")


def load(with_xy: bool = True) -> np.ndarray:
    """Return a structured array of the pad map."""
    pattern = os.path.join(MAPDIR, "Mapping" if with_xy else "", "connector_*.txt")
    rows = []
    for path in sorted(glob.glob(pattern)):
        for line in open(path).read().splitlines()[1:]:
            p = line.split()
            if with_xy and len(p) >= 9:
                rows.append((int(p[0]), int(p[2]), int(p[3]), -1,
                             float(p[5]), float(p[6]), float(p[7]), float(p[8])))
            elif not with_xy and len(p) >= 6:
                r, phi = float(p[4]), float(p[5])
                rows.append((int(p[0]), int(p[1]), int(p[2]), int(p[3]),
                             r * math.cos(phi), r * math.sin(phi), r, phi))
    return np.array(rows, dtype=[("conn", "i4"), ("chan", "i4"),
                                 ("pad", "i4"), ("abs", "i4"),
                                 ("x", "f8"), ("y", "f8"),
                                 ("r", "f8"), ("phi", "f8")])


def report(a: np.ndarray, label: str):
    print("\n" + "=" * 74)
    print(f"READOUT PAD MAP  --  {label}")
    print("=" * 74)
    print(f"  pads            : {len(a)}")
    print(f"  connectors      : {len(np.unique(a['conn']))} "
          f"x {np.bincount(a['conn']).max()} channels")
    print(f"  radius          : {a['r'].min():.3f} .. {a['r'].max():.3f} mm")
    print(f"  phi             : {math.degrees(a['phi'].min()):.4f} .. "
          f"{math.degrees(a['phi'].max()):.4f} deg  "
          f"({a['phi'].min():.5f} .. {a['phi'].max():.5f} rad)")

    rings = np.unique(np.round(a["r"], 3))
    dr = np.diff(rings)
    print(f"  radial rings    : {len(rings)}")
    print(f"  radial pitch    : {np.median(dr):.4f} mm "
          f"(min {dr.min():.4f}, max {dr.max():.4f})")

    cnt = collections.Counter(np.round(a["r"], 3))
    print(f"  pads per ring   : {min(cnt.values())} (inner) .. "
          f"{max(cnt.values())} (outer)")

    print(f"\n  {'ring r [mm]':>12} {'n pads':>7} {'dphi [deg]':>11} "
          f"{'arc pitch [mm]':>15}")
    for rr in [rings[0], rings[len(rings) // 4], rings[len(rings) // 2],
               rings[3 * len(rings) // 4], rings[-1]]:
        ph = np.sort(a["phi"][np.abs(a["r"] - rr) < 1e-3])
        dp = np.median(np.diff(ph)) if len(ph) > 2 else float("nan")
        print(f"  {rr:>12.3f} {len(ph):>7} {math.degrees(dp):>11.4f} "
              f"{dp * rr:>15.3f}")

    print(f"\n  => pad cell is ~{np.median(dr):.2f} mm (radial) x "
          f"~{np.median([np.median(np.diff(np.sort(a['phi'][np.abs(a['r']-rr)<1e-3])))*rr for rr in rings[1:-1]]):.2f} mm "
          f"(azimuthal), i.e. roughly area-equalized")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    a = load(with_xy=True)
    report(a, "design/mapping/Mapping/connector_*.txt  (9-column, has X/Y)")

    b = load(with_xy=False)
    if len(b):
        report(b, "design/mapping/connector_*.txt  (7-column, has AbsIndex)")
        # Do the two revisions describe the same pads? Match by nearest
        # neighbour, NOT by rounded coordinates: the 7-column file has no X/Y,
        # so x,y are recomputed from r,phi and differ from the printed 9-column
        # values in the last digit. (An earlier rounded-key comparison here
        # reported a spurious "79/1280 pads disagree", which propagated into
        # the campaign docs as a risk item -- the pad planes are identical.)
        PA = np.c_[a["x"], a["y"]]
        PB = np.c_[b["x"], b["y"]]
        dmin = np.sqrt(((PA[:, None, :] - PB[None, :, :]) ** 2).sum(-1)).min(axis=1)
        n_off = int((dmin > 0.01).sum())
        print(f"\n  pad POSITIONS: max nearest-pad mismatch "
              f"{dmin.max()*1000:.1f} um; {n_off}/{len(a)} pads off by >10 um")
        if n_off == 0:
            print("  => the two revisions describe the SAME pad plane.")

        # ... but do they assign the same channel to the same pad?
        ib = {(int(c), int(ch)): i for i, (c, ch) in enumerate(zip(b["conn"], b["chan"]))}
        same = sum(1 for j, (c, ch) in enumerate(zip(a["conn"], a["chan"]))
                   if (int(c), int(ch)) in ib
                   and abs(b["r"][ib[(int(c), int(ch))]] - a["r"][j]) < 0.01
                   and abs(b["phi"][ib[(int(c), int(ch))]] - a["phi"][j]) < 1e-4)
        print(f"  pad CHANNEL ASSIGNMENT: (connector, channel) -> same pad for "
              f"{same}/{len(a)} channels ({100*same/len(a):.1f}%)")
        if same < len(a):
            print("  NOTE: the revisions disagree on the READOUT ORDER, not the\n"
                  "  geometry. Everything channel-level -- VMM neighbour logic,\n"
                  "  dead-channel masks, per-channel test-beam comparisons --\n"
                  "  depends on which one the DAQ uses. See vmm/nl_map.py.")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 9))
        sc = ax.scatter(a["x"], a["y"], c=a["conn"], s=6, cmap="tab10")
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title("P2 readout pads, coloured by connector (0-9)")
        plt.colorbar(sc, ax=ax, label="connector", shrink=0.7)
        ax.grid(alpha=0.25, lw=0.4)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"\nWrote {args.plot}")


if __name__ == "__main__":
    main()
