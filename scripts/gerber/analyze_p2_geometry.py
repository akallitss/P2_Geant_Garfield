#!/usr/bin/env python3
"""
Extract P2 wedge detector dimensions from the fabrication gerbers.

Reports the board outline (bounding box, wedge arc radii and opening angle),
the readout strip/pad pitch from the copper apertures, and the drill map --
i.e. everything needed to write the Geant4 solid.

Usage:
    python scripts/gerber/analyze_p2_geometry.py                 # default gerber set
    python scripts/gerber/analyze_p2_geometry.py --plot out.pdf  # + outline drawing
    python scripts/gerber/analyze_p2_geometry.py --dir <gerberdir>
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gerber_outline import (parse, build_contours, polygon_area, bbox_of,
                            fit_circle)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DEFAULT_DIR = os.path.join(REPO, "design", "gerbers", "P2_BASKET_Apr26", "Gerber")


def hr(title):
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def analyze_outline(path, plot_ax=None):
    hr(f"OUTLINE  --  {os.path.basename(path)}")
    gf = parse(path)
    print(f"  units             : {gf.unit}")
    print(f"  stroked segments  : {len(gf.segments)}")
    print(f"  regions           : {len(gf.regions)}")
    x0, y0, x1, y1 = gf.bbox()
    print(f"  overall bbox      : X [{x0:.3f}, {x1:.3f}] mm   "
          f"Y [{y0:.3f}, {y1:.3f}] mm")
    print(f"  overall size      : {x1 - x0:.3f} x {y1 - y0:.3f} mm")

    contours = build_contours(gf.segments)
    contours.sort(key=polygon_area, reverse=True)
    print(f"  chained contours  : {len(contours)}")

    print("\n  Largest contours (candidate board profile):")
    print(f"    {'#':>3} {'pts':>6} {'area/mm^2':>12} {'w x h [mm]':>26}")
    for i, c in enumerate(contours[:6]):
        bx0, by0, bx1, by1 = bbox_of(c)
        print(f"    {i:>3} {len(c):>6} {polygon_area(c):>12.1f} "
              f"{f'{bx1-bx0:.2f} x {by1-by0:.2f}':>26}")

    main = contours[0]

    # --- wedge characterisation -------------------------------------------
    # Fit a circle to the arcs of the main contour; a wedge (annular sector)
    # should show two dominant radii about a common apex.
    arcs = [s for s in gf.segments if s.kind == "arc" and s.radius > 5.0]
    if arcs:
        print("\n  Long arcs in the profile (r > 5 mm) -- wedge radii:")
        radii = Counter()
        for s in arcs:
            radii[round(s.radius, 2)] += 1
        for r, n in sorted(radii.items(), reverse=True)[:12]:
            print(f"    r = {r:>10.2f} mm   ({n} segment(s))")

        big = max(arcs, key=lambda s: s.length())
        print(f"\n  Longest arc: r = {big.radius:.3f} mm, "
              f"centre = ({big.cx:.3f}, {big.cy:.3f}), "
              f"arc length = {big.length():.2f} mm")
        a0 = math.degrees(math.atan2(big.y0 - big.cy, big.x0 - big.cx))
        a1 = math.degrees(math.atan2(big.y1 - big.cy, big.x1 - big.cx))
        span = abs((a1 - a0 + 180) % 360 - 180)
        print(f"  Longest arc angular span: {span:.3f} deg "
              f"({a0:.2f} -> {a1:.2f})")

    cx, cy, r, rms = fit_circle(main)
    print(f"\n  Circle fit to main contour: centre=({cx:.2f}, {cy:.2f}) "
          f"r={r:.2f} mm  rms={rms:.2f} mm"
          f"{'  <-- poor fit, not a pure annulus' if rms > 5 else ''}")

    if plot_ax is not None:
        for c in contours[:200]:
            xs = [p[0] for p in c]
            ys = [p[1] for p in c]
            plot_ax.plot(xs, ys, "-", lw=0.6, color="k")
        plot_ax.set_aspect("equal")
        plot_ax.set_xlabel("x [mm]")
        plot_ax.set_ylabel("y [mm]")
        plot_ax.set_title(os.path.basename(path))

    return gf, contours


def analyze_copper(path):
    hr(f"COPPER  --  {os.path.basename(path)}")
    gf = parse(path)
    x0, y0, x1, y1 = gf.bbox()
    print(f"  bbox              : {x1-x0:.3f} x {y1-y0:.3f} mm  "
          f"(X {x0:.2f}..{x1:.2f}, Y {y0:.2f}..{y1:.2f})")
    print(f"  segments={len(gf.segments)}  flashes={len(gf.flashes)}  "
          f"regions={len(gf.regions)}  apertures={len(gf.apertures)}")

    by_fn = {}
    for ap in gf.apertures.values():
        by_fn.setdefault(ap.function or "(none)", []).append(ap)
    print("\n  Apertures by function:")
    for fn, aps in sorted(by_fn.items()):
        sizes = Counter(round(a.size, 4) for a in aps)
        top = ", ".join(f"{s}mm x{n}" for s, n in sizes.most_common(6))
        print(f"    {fn:<20} n={len(aps):<5} {top}")

    # Strip pads on this design are RotRect apertures; their rotation angle is
    # the third macro parameter and tells us the strip fan-out angles.
    rot = [a for a in gf.apertures.values() if a.template == "RotRect"]
    if rot:
        print(f"\n  RotRect (radial strip pad) apertures: {len(rot)}")
        dims = Counter((round(a.params[0], 3), round(a.params[1], 3))
                       for a in rot if len(a.params) >= 2)
        for (w, h), n in dims.most_common():
            print(f"    {w} x {h} mm   x{n}")
        angles = sorted({round(a.params[2], 3) for a in rot if len(a.params) >= 3})
        print(f"  distinct rotation angles: {len(angles)}")
        if len(angles) > 1:
            print(f"    range {min(angles):.3f} .. {max(angles):.3f} deg")
            deltas = [round(b - a, 4) for a, b in zip(angles, angles[1:])]
            print(f"    step  {min(deltas):.4f} .. {max(deltas):.4f} deg "
                  f"(median {sorted(deltas)[len(deltas)//2]:.4f})")

    # Conductor traces = readout strips; their aperture width is the strip width
    cond = [a for a in gf.apertures.values() if a.function == "Conductor"]
    if cond:
        print("\n  Conductor (strip) widths: "
              + ", ".join(f"{round(a.size, 4)}mm" for a in
                          sorted(cond, key=lambda a: a.size)[:12]))
    return gf


def analyze_drill(path):
    hr(f"DRILL  --  {os.path.basename(path)}")
    tools, counts, cur = {}, Counter(), None
    unit_mm = True
    with open(path, errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("INCH"):
                unit_mm = False
            if s.startswith("T") and "C" in s:
                try:
                    tid = int(s[1:s.index("C")])
                    tools[tid] = float(s[s.index("C") + 1:])
                except ValueError:
                    pass
            elif s.startswith("T") and s[1:].isdigit():
                cur = int(s[1:])
            elif s.startswith("X") and cur is not None:
                counts[cur] += 1
    print(f"  units: {'mm' if unit_mm else 'inch'}")
    print(f"  {'tool':>6} {'dia':>10} {'holes':>8}")
    total = 0
    for tid in sorted(tools):
        print(f"  {'T'+str(tid):>6} {tools[tid]:>10.4f} {counts.get(tid, 0):>8}")
        total += counts.get(tid, 0)
    print(f"  total holes: {total}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=DEFAULT_DIR, help="gerber directory")
    ap.add_argument("--drill", default=None, help="NC drill file")
    ap.add_argument("--plot", default=None, help="write outline plot to this PDF/PNG")
    ap.add_argument("--skip-copper", action="store_true",
                    help="skip the (slow, multi-MB) copper layers")
    args = ap.parse_args()

    files = sorted(os.listdir(args.dir))
    edge = [f for f in files if "Edge_Cuts" in f and f.endswith(".gbr")]
    cu = [f for f in files if f.endswith(".gbr") and "_Cu" in f and "Edge" not in f]

    fig = axes = None
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, max(1, len(edge)), figsize=(9 * max(1, len(edge)), 9))
        if len(edge) <= 1:
            axes = [axes]

    for i, f in enumerate(edge):
        analyze_outline(os.path.join(args.dir, f),
                        plot_ax=axes[i] if axes is not None else None)

    if not args.skip_copper:
        for f in cu:
            analyze_copper(os.path.join(args.dir, f))

    drill = args.drill
    if drill is None:
        ncd = os.path.join(os.path.dirname(args.dir.rstrip("/")), "NCDrill")
        if os.path.isdir(ncd):
            for f in sorted(os.listdir(ncd)):
                if f.endswith(".drl"):
                    drill = os.path.join(ncd, f)
                    break
    if drill and os.path.exists(drill):
        analyze_drill(drill)

    if args.plot:
        import matplotlib.pyplot as plt
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"\nWrote {args.plot}")


if __name__ == "__main__":
    main()
