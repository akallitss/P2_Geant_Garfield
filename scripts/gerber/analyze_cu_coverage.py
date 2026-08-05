#!/usr/bin/env python3
"""
Copper coverage of the P2 wedge copper layers — exact, via shapely.

The B.Cu layer is *signal lines*, not a ground plane (Alexandra 2026-08-05),
so the Geant4 model must not use a solid 18 µm sheet. This script builds the
exact copper geometry from the gerbers (strokes buffered to their aperture
width with round caps, flashes as discs/rectangles, G36 regions as polygons —
both layers are purely additive, no LPC), unions it, and reports the covered
area fraction over

  * the board outline        (what the model's PCB_Cu_* prisms span), and
  * the active area          (r 119.87..589.80 mm, the pad field — the
                              region sampled by tracks in the acceptance),

plus a radial profile over the active area. The active-area numbers are
folded into the simulation as effective density scales
(SimConfig p2_bcu_coverage / p2_fcu_coverage; docs/P2_MODEL.md).

Result (2026-08-05, this script):
    B_Cu: 0.181 board / 0.174 active   (profile 0.03 -> 0.26 inner->outer)
    F_Cu: 0.937 board / 0.983 active

Usage:
    python scripts/gerber/analyze_cu_coverage.py
"""

from __future__ import annotations

import math
import os
import sys

from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "model"))
from gerber_outline import parse                      # noqa: E402
from p2_model import wedge_outline, BOARD, ACTIVE     # noqa: E402

REPO = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
GDIR = os.path.join(REPO, "design", "gerbers", "P2_BASKET_Apr26", "Gerber")

# Active area of the pad field: radii from the bulk mask (exact), sector
# edges from the pad map (phi 0.97..59.03 deg, docs/P2_GEOMETRY.md).
ACT_PHI0, ACT_PHI1 = 0.97, 59.03


def sector(r_in: float, r_out: float, phi0: float, phi1: float,
           n: int = 256) -> Polygon:
    a0, a1 = math.radians(phi0), math.radians(phi1)
    pts = [(r_out * math.cos(a0 + (a1 - a0) * i / n),
            r_out * math.sin(a0 + (a1 - a0) * i / n)) for i in range(n + 1)]
    pts += [(r_in * math.cos(a1 - (a1 - a0) * i / n),
             r_in * math.sin(a1 - (a1 - a0) * i / n)) for i in range(n + 1)]
    return Polygon(pts)


def copper_union(gbr):
    parts = []
    for s in gbr.segments:
        w = gbr.apertures[s.aperture].size
        pts = s.sample(48) if s.kind == "arc" else [(s.x0, s.y0), (s.x1, s.y1)]
        ln = LineString(pts)
        if ln.length > 1e-9:
            parts.append(ln.buffer(w / 2))          # round caps/joins = gerber
        else:
            parts.append(Point(pts[0]).buffer(w / 2))
    for f in gbr.flashes:
        ap = gbr.apertures[f.aperture]
        if ap.template == "R" and len(ap.params) >= 2:
            parts.append(box(f.x - ap.params[0] / 2, f.y - ap.params[1] / 2,
                             f.x + ap.params[0] / 2, f.y + ap.params[1] / 2))
        else:                                       # C/O/P/macros -> disc
            parts.append(Point(f.x, f.y).buffer(ap.size / 2))
    for reg in gbr.regions:
        if len(reg) >= 3:
            p = Polygon(reg)
            if not p.is_valid:
                p = p.buffer(0)
            parts.append(p)
    return unary_union(parts)


def main():
    board = Polygon(wedge_outline(**BOARD))
    active = sector(ACTIVE["r_in"], ACTIVE["r_out"], ACT_PHI0, ACT_PHI1)

    for layer in ("B_Cu", "F_Cu"):
        g = parse(os.path.join(GDIR, f"P2_BASKET-{layer}.gbr"))
        cu = copper_union(g)
        f_board = cu.intersection(board).area / board.area
        f_act = cu.intersection(active).area / active.area
        print(f"\n{layer}: {len(g.segments)} segments, {len(g.flashes)} "
              f"flashes, {len(g.regions)} regions; union {cu.area/100:.1f} cm2")
        print(f"  coverage over board outline : {f_board:.4f}")
        print(f"  coverage over ACTIVE area   : {f_act:.4f}   <-- sim input")
        print("  radial profile (active area):")
        for i in range(10):
            r_lo = ACTIVE["r_in"] + (ACTIVE["r_out"] - ACTIVE["r_in"]) * i / 10
            r_hi = ACTIVE["r_in"] + (ACTIVE["r_out"] - ACTIVE["r_in"]) * (i + 1) / 10
            ring = sector(r_lo, r_hi, ACT_PHI0, ACT_PHI1)
            print(f"    r {r_lo:5.1f}..{r_hi:5.1f} mm : "
                  f"{cu.intersection(ring).area / ring.area:.4f}")


if __name__ == "__main__":
    main()
