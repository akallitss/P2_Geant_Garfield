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

Result (2026-08-05, this script, BEFORE the aperture fix below):
    B_Cu: 0.181 board / 0.174 active   (profile 0.03 -> 0.26 inner->outer)
    F_Cu: 0.937 board / 0.983 active   <-- BOTH WRONG, see below

**Aperture handling fixed 2026-08-07.** Flashes were buffered as discs of
radius `ApertureDef.size / 2`, and `size` is `max(params)` -- which for
KiCad's `RotRect` macro apertures is the ROTATION ANGLE in degrees
(`%ADD19RotRect,0.300000X1.800000X302.763000*%` -> "size" 302.763 mm). That
painted a ~151 mm-radius disc at each of the 1420 connector-pin flashes on
F_Cu. Apertures tagged `NonConductor` and `Profile` were also counted as
copper. B_Cu was barely affected (it has no macro apertures); F_Cu was
inflated badly. Corrected values, cross-checked independently against an
exact rasterization in `extract_readout_pattern.py`:

    B_Cu: 0.180 board / 0.171 over the pad field
    F_Cu: 0.787 board / 0.913 over the pad field

**Prefer `scripts/gerber/extract_readout_pattern.py`.** It is the script the
Geant4 geometry is actually generated from: it emits `include/P2PadMap.hh`
with the exact pad grid plus per-radial-band coverage, needs no shapely, and
the numbers above come from it. This file is kept as the independent
cross-check -- a different algorithm (exact vector unions rather than a
raster) on the same input.

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


# Apertures that are drawn on a copper layer but are not copper. KiCad tags
# them via %TA.AperFunction%; `Profile` is the board outline, which this file
# used to buffer into the union as 0.1 mm-wide "copper" all round the panel.
NON_COPPER_FUNCTIONS = ("NonConductor", "Profile")


def _flash_shape(f, ap):
    """One flashed aperture as a shapely polygon, honouring its TEMPLATE.

    `ap.size` is `max(params)`, which is only a size for a plain circle. For
    R/O it is the longer side; for P it may be the vertex count or a rotation;
    and for KiCad's `RotRect` macro the third parameter is a ROTATION ANGLE
    IN DEGREES:

        %ADD19RotRect,0.300000X1.800000X302.763000*%

    Buffering that as a disc of radius `size/2` painted a 151 mm-radius blob
    at each of the 1420 connector-pin flashes, which is what used to lift
    F_Cu to 0.937 board / 0.983 active. Both numbers were wrong.
    """
    p = ap.params
    t = ap.template
    if t == "R" and len(p) >= 2:
        return box(f.x - p[0] / 2, f.y - p[1] / 2, f.x + p[0] / 2, f.y + p[1] / 2)
    if t == "O" and len(p) >= 2:                    # obround = stadium
        w, h = p[0], p[1]
        if w >= h:
            return LineString([(f.x - (w - h) / 2, f.y),
                               (f.x + (w - h) / 2, f.y)]).buffer(h / 2)
        return LineString([(f.x, f.y - (h - w) / 2),
                           (f.x, f.y + (h - w) / 2)]).buffer(w / 2)
    if t == "P" and len(p) >= 2:                    # regular polygon
        n = int(p[1])
        rot = math.radians(p[2]) if len(p) >= 3 else 0.0
        r = p[0] / 2
        return Polygon([(f.x + r * math.cos(rot + 2 * math.pi * k / n),
                         f.y + r * math.sin(rot + 2 * math.pi * k / n))
                        for k in range(n)])
    if t == "RotRect" and len(p) >= 3:              # KiCad macro: L x W x angle
        from shapely.affinity import rotate
        return rotate(box(f.x - p[0] / 2, f.y - p[1] / 2,
                          f.x + p[0] / 2, f.y + p[1] / 2), p[2],
                      origin=(f.x, f.y))
    # C, and any unknown macro: a disc on the FIRST parameter (the leading
    # dimension in every KiCad macro), never on max(params).
    return Point(f.x, f.y).buffer((p[0] / 2) if p else 0.0)


def copper_union(gbr):
    parts = []
    for s in gbr.segments:
        ap = gbr.apertures[s.aperture]
        if ap.function in NON_COPPER_FUNCTIONS:
            continue
        w = ap.params[0] if ap.params else 0.0      # stroking width, not size
        pts = s.sample(48) if s.kind == "arc" else [(s.x0, s.y0), (s.x1, s.y1)]
        ln = LineString(pts)
        if ln.length > 1e-9:
            parts.append(ln.buffer(w / 2))          # round caps/joins = gerber
        else:
            parts.append(Point(pts[0]).buffer(w / 2))
    for f in gbr.flashes:
        ap = gbr.apertures[f.aperture]
        if ap.function in NON_COPPER_FUNCTIONS:
            continue
        parts.append(_flash_shape(f, ap))
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
