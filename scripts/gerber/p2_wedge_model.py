#!/usr/bin/env python3
"""
Idealised P2 wedge outline, and a check of it against the fabrication gerber.

The P2_BASKET board profile recovered from P2_BASKET-Edge_Cuts.gbr is:

    apex        : gerber origin (0, 0)  -- the MESA beam axis
    r_inner     :  95.0 mm
    r_outer     : 650.0 mm
    sector      : phi = 0 deg .. 60 deg
    edge offset : both radial edges pushed OUTWARD by 10.0 mm (perpendicular),
                  so the board is slightly wider than a true 60 deg sector
    top cut     : y <= 540.0 mm  (flat chord across the wide end)

Run this file to regenerate the comparison figure:

    python scripts/gerber/p2_wedge_model.py --out wedge_outline.pdf
"""

from __future__ import annotations

import math

# ── Recovered parameters (mm / deg) ─────────────────────────────────────────
R_INNER = 95.0
R_OUTER = 650.0
PHI0 = 0.0
PHI1 = 60.0
EDGE_OFFSET = 10.0
Y_TOP_CUT = 540.01
ROUTE_WIDTH = 2.0        # separation between board edge and panel route

# Fabrication panel the board sits in (not part of the detector)
PANEL = (38.594, -39.996, 648.595, 570.004)


def corner_angles() -> dict:
    """Angular positions where each boundary element meets the next."""
    return {
        # inner arc endpoints (offset edges cut it wider than 60 deg)
        "in_lo": -math.degrees(math.asin(EDGE_OFFSET / R_INNER)),
        "in_hi": PHI1 + math.degrees(math.asin(EDGE_OFFSET / R_INNER)),
        # outer arc endpoints
        "out_lo": -math.degrees(math.asin(EDGE_OFFSET / R_OUTER)),
        "out_hi": math.degrees(math.asin(Y_TOP_CUT / R_OUTER)),   # stopped by top cut
    }


def outline(n_arc: int = 400) -> list[tuple[float, float]]:
    """Closed polygon of the board profile, counter-clockwise, in mm.

    Traversal: inner-arc/bottom corner -> along the phi=0 offset edge -> outer
    arc -> top chord -> along the phi=60 offset edge -> back down the inner arc.
    """
    c = corner_angles()

    def arc_pts(r, a0, a1, n):
        return [(r * math.cos(math.radians(a)), r * math.sin(math.radians(a)))
                for a in [a0 + (a1 - a0) * i / n for i in range(n + 1)]]

    pts = []
    # bottom edge: from the inner arc out to the outer arc, along y = -EDGE_OFFSET
    pts.append((R_INNER * math.cos(math.radians(c["in_lo"])), -EDGE_OFFSET))
    pts.append((R_OUTER * math.cos(math.radians(c["out_lo"])), -EDGE_OFFSET))
    # outer arc, CCW, until the top chord cuts it
    pts += arc_pts(R_OUTER, c["out_lo"], c["out_hi"], n_arc)
    # top chord across to the phi=60 offset edge
    phi1 = math.radians(PHI1)
    x_at_top = (Y_TOP_CUT * math.cos(phi1) - EDGE_OFFSET) / math.sin(phi1)
    pts.append((x_at_top, Y_TOP_CUT))
    # down the phi=60 offset edge to the inner arc
    a = math.radians(c["in_hi"])
    pts.append((R_INNER * math.cos(a), R_INNER * math.sin(a)))
    # inner arc, CW, back to the start
    pts += arc_pts(R_INNER, c["in_hi"], c["in_lo"], n_arc)
    return pts


def area_mm2() -> float:
    pts = outline(4000)
    a = 0.0
    for i in range(len(pts)):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % len(pts)]
        a += x0 * y1 - x1 * y0
    return abs(a) / 2


def summary() -> str:
    da_in = math.degrees(math.asin(EDGE_OFFSET / R_INNER))
    da_out = math.degrees(math.asin(EDGE_OFFSET / R_OUTER))
    return f"""P2 wedge (P2_BASKET, Version_Apr26/FABRICATION_P2)
  apex / beam axis      : gerber (0, 0)
  inner radius          : {R_INNER} mm
  outer radius          : {R_OUTER} mm
  nominal sector        : {PHI1 - PHI0:.0f} deg  (phi {PHI0} .. {PHI1})
  radial edges offset   : {EDGE_OFFSET} mm outward (perpendicular)
  -> inner arc spans    : {PHI0 - da_in:+.3f} .. {PHI1 + da_in:+.3f} deg  ({PHI1 - PHI0 + 2*da_in:.3f} deg)
  -> outer arc spans    : {PHI0 - da_out:+.3f} .. {PHI1 + da_out:+.3f} deg  ({PHI1 - PHI0 + 2*da_out:.3f} deg)
  top chord cut         : y <= {Y_TOP_CUT} mm  (outer arc truncated between
                          x = {math.sqrt(R_OUTER**2 - Y_TOP_CUT**2):.1f} mm and the 60 deg edge)
  radial extent         : {R_OUTER - R_INNER} mm
  board area            : {area_mm2()/1e2:.1f} cm^2
  panel (fab only)      : x {PANEL[0]}..{PANEL[2]}, y {PANEL[1]}..{PANEL[3]} mm
  route width           : {ROUTE_WIDTH} mm between board edge and panel
"""


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="wedge_outline.pdf")
    ap.add_argument("--gerber", default=None,
                    help="Edge_Cuts.gbr to overlay (default: repo copy)")
    args = ap.parse_args()

    print(summary())

    import os
    import sys
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    repo = os.path.abspath(os.path.join(here, "..", ".."))
    gpath = args.gerber or os.path.join(
        repo, "design/gerbers/P2_BASKET_Apr26/Gerber/P2_BASKET-Edge_Cuts.gbr")

    fig, ax = plt.subplots(figsize=(10, 10))

    if os.path.exists(gpath):
        from gerber_outline import parse
        gf = parse(gpath)
        for s in gf.segments:
            p = s.sample(24)
            ax.plot([q[0] for q in p], [q[1] for q in p],
                    "-", color="0.55", lw=0.7, zorder=1)
        ax.plot([], [], "-", color="0.55", lw=0.7, label="Edge_Cuts gerber")

    pts = outline(600)
    ax.plot([p[0] for p in pts] + [pts[0][0]],
            [p[1] for p in pts] + [pts[0][1]],
            "-", color="crimson", lw=1.8, zorder=3,
            label="idealised model")

    for r in (R_INNER, R_OUTER):
        a = [math.radians(t) for t in
             [i * (PHI1 - PHI0 + 20) / 200 + PHI0 - 10 for i in range(201)]]
        ax.plot([r * math.cos(t) for t in a], [r * math.sin(t) for t in a],
                ":", color="steelblue", lw=0.9, zorder=2)
    for phi in (PHI0, PHI1):
        t = math.radians(phi)
        ax.plot([R_INNER * math.cos(t) * 0, R_OUTER * math.cos(t)],
                [0, R_OUTER * math.sin(t)], ":", color="steelblue", lw=0.9, zorder=2)
    ax.plot([], [], ":", color="steelblue", lw=0.9,
            label=f"true {PHI1-PHI0:.0f}° sector, r={R_INNER}/{R_OUTER} mm")

    ax.plot(0, 0, "+", color="k", ms=12, mew=1.5, zorder=4)
    ax.annotate("apex = beam axis", (0, 0), textcoords="offset points",
                xytext=(12, -14), fontsize=9)

    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title("P2 wedge outline — gerber vs. idealised model")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.25, lw=0.4)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
