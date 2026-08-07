#!/usr/bin/env python3
"""
The P2 readout copper, measured out of the production gerbers.

Two jobs, both feeding the Geant4 geometry:

1. **The pad pattern.**  The F.Cu pad field is an *exactly regular polar grid*
   — 42 rings on a constant 11.4286 mm radial pitch, each ring a set of
   uniformly spaced annular-sector pads.  This script fits that grid to the
   1280 pad polygons in the gerber, checks the regularity claim ring by ring,
   and writes the table out as `include/P2PadMap.hh` so `DetectorConstruction`
   can build the real copper instead of a density-scaled sheet.

2. **The copper coverage everywhere else.**  Outside the pad field the artwork
   is irregular (fan-out, connector footprints, via fields, and on B.Cu 30789
   individual 0.125 mm traces), so it stays homogenized — but per *zone*, not
   as one board-wide average.  Coverage is measured here by stratified Monte
   Carlo against the exact vector geometry.

Why not `analyze_cu_coverage.py`: that script has two aperture bugs that
inflate F.Cu badly (see `--compare`).  It buffers every non-"R" flash as a disc
of radius `ap.size/2`, and `ApertureDef.size` is `max(params)` — so for KiCad's
`RotRect` macro apertures it picks up the **rotation angle**:

    %ADD19RotRect,0.300000X1.800000X302.763000*%   ->  "size" = 302.763 mm

which paints a 151 mm-radius disc at each of the 1420 connector-pin flashes.
It also counts apertures tagged `NonConductor` and `Profile` as copper.  This
script honours the aperture template and the `.AperFunction` attribute.

Usage:
    python scripts/gerber/extract_readout_pattern.py              # report
    python scripts/gerber/extract_readout_pattern.py --write      # + header
    python scripts/gerber/extract_readout_pattern.py --compare    # + old-script bug
    python scripts/gerber/extract_readout_pattern.py --plot p.png
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "scripts", "model"))

from gerber_outline import parse                        # noqa: E402
from p2_model import wedge_outline, BOARD, SECTOR_DEG as SECTOR_DEG_  # noqa: E402

GDIR = os.path.join(REPO, "design", "gerbers", "P2_BASKET_Apr26", "Gerber")
HEADER = os.path.join(REPO, "include", "P2PadMap.hh")

# A pad is a G36 region of roughly 11.3 x 11.8 mm.  The band is wide enough to
# admit the fat inner ring (144 mm2) and narrow enough to exclude everything
# else on the layer: the next-largest non-pad region is 25 mm2, the next
# smallest above the band is 1580 mm2.
PAD_AREA_LO, PAD_AREA_HI = 100.0, 200.0

# Apertures whose flashes/strokes are not copper.  Both attributes are written
# by KiCad; `Profile` is the board outline drawn onto the copper layer.
NON_COPPER_FUNCTIONS = ("NonConductor", "Profile")


# ─────────────────────────────────────────────────────────────────────────────
# Exact vector geometry of one copper layer
# ─────────────────────────────────────────────────────────────────────────────
#
# Every feature becomes one of three primitives with an exact point-inside
# test.  Nothing is rasterized, so there is no pixel-edge convention to get
# wrong; the only approximation is the Monte-Carlo sampling of the zones,
# whose statistical error is reported alongside every number.

class Capsule:
    """Stroked path element: all points within `r` of the segment p0->p1."""
    __slots__ = ("x0", "y0", "dx", "dy", "len2", "r", "bbox")

    def __init__(self, x0, y0, x1, y1, r):
        self.x0, self.y0 = x0, y0
        self.dx, self.dy = x1 - x0, y1 - y0
        self.len2 = self.dx * self.dx + self.dy * self.dy
        self.r = r
        self.bbox = (min(x0, x1) - r, min(y0, y1) - r,
                     max(x0, x1) + r, max(y0, y1) + r)

    def contains(self, x, y):
        px, py = x - self.x0, y - self.y0
        if self.len2 > 0.0:
            t = np.clip((px * self.dx + py * self.dy) / self.len2, 0.0, 1.0)
            px = px - t * self.dx
            py = py - t * self.dy
        return px * px + py * py <= self.r * self.r


class Rect:
    """Flashed rectangle, optionally rotated (KiCad's RotRect macro)."""
    __slots__ = ("cx", "cy", "hw", "hh", "c", "s", "bbox")

    def __init__(self, cx, cy, w, h, rot_deg=0.0):
        self.cx, self.cy = cx, cy
        self.hw, self.hh = w / 2.0, h / 2.0
        a = math.radians(rot_deg)
        self.c, self.s = math.cos(a), math.sin(a)
        ex = abs(self.hw * self.c) + abs(self.hh * self.s)
        ey = abs(self.hw * self.s) + abs(self.hh * self.c)
        self.bbox = (cx - ex, cy - ey, cx + ex, cy + ey)

    def contains(self, x, y):
        px, py = x - self.cx, y - self.cy
        u = px * self.c + py * self.s          # rotate into the pad frame
        v = -px * self.s + py * self.c
        return (np.abs(u) <= self.hw) & (np.abs(v) <= self.hh)


class Poly:
    """G36 region, filled by scanline.

    A point-in-polygon test per pixel is the obvious thing and it is far too
    slow here: the fan-out ground pour spans the whole 650 mm outer band, so
    at a 50 um pitch its bounding box alone is ~1.7e8 pixels and testing each
    against ~40 edges is ~7e9 operations.  Scanline instead costs
    (rows x edges) to find the crossings and then fills whole spans, which is
    the same answer several thousand times faster.
    """
    __slots__ = ("ax", "ay", "bx", "by", "ylo", "yhi", "bbox")

    # The board-edge copper pours are single KiCad zones with every cutout
    # stitched into one contour: 155979 vertices on F.Cu, 91471 on B.Cu, each
    # spanning the whole board.  Crossing every edge with every scanline in
    # one shot is a (rows x edges) array -- ~10^8 doubles for those two -- so
    # rows are processed in chunks and, within a chunk, only the edges whose
    # y-extent reaches it are considered.  The edges of a stitched pour are
    # individually tiny, so that selection is what makes this tractable.
    ROW_CHUNK = 256
    BIG_EDGES = 2_000

    def __init__(self, pts):
        p = np.asarray(pts, dtype=float)
        q = np.roll(p, -1, axis=0)
        keep = p[:, 1] != q[:, 1]              # horizontal edges never cross
        ax, ay = p[keep, 0], p[keep, 1]
        bx, by = q[keep, 0], q[keep, 1]
        order = np.argsort(np.minimum(ay, by))          # for range selection
        self.ax, self.ay = ax[order], ay[order]
        self.bx, self.by = bx[order], by[order]
        self.ylo = np.minimum(self.ay, self.by)
        self.yhi = np.maximum(self.ay, self.by)
        self.bbox = (p[:, 0].min(), p[:, 1].min(), p[:, 0].max(), p[:, 1].max())

    def mark(self, img, x0, y0, pitch, i0, i1, j0, j1):
        if len(self.ax) == 0:
            return
        if len(self.ax) <= self.BIG_EDGES:
            self._mark_rows(img, x0, y0, pitch, i0, i1, j0, j1, None)
            return
        for jj in range(j0, j1, self.ROW_CHUNK):
            jn = min(jj + self.ROW_CHUNK, j1)
            ylo = y0 + (jj + 0.5) * pitch
            yhi = y0 + (jn - 0.5) * pitch
            n = np.searchsorted(self.ylo, yhi, side="right")
            sel = np.flatnonzero(self.yhi[:n] >= ylo)
            if len(sel):
                self._mark_rows(img, x0, y0, pitch, i0, i1, jj, jn, sel)

    def _mark_rows(self, img, x0, y0, pitch, i0, i1, j0, j1, sel):
        ax, ay = (self.ax, self.ay) if sel is None else (self.ax[sel], self.ay[sel])
        bx, by = (self.bx, self.by) if sel is None else (self.bx[sel], self.by[sel])
        ys = y0 + (np.arange(j0, j1) + 0.5) * pitch
        # crossings of every edge with every scanline, +inf where there is none
        straddle = (ay[None, :] > ys[:, None]) != (by[None, :] > ys[:, None])
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = ax + (ys[:, None] - ay) * (bx - ax) / (by - ay)
        xs = np.sort(np.where(straddle, xint, np.inf), axis=1)

        # even-odd: consecutive sorted crossings bound an interior span.
        # Accumulate span starts/ends into a difference image, one pair index
        # at a time, then a single cumsum turns it into the filled mask.
        w = i1 - i0
        diff = np.zeros((j1 - j0, w + 1), dtype=np.int16)
        rows = np.arange(j1 - j0)
        n_pairs = xs.shape[1] // 2
        for p in range(n_pairs):
            lo, hi = xs[:, 2 * p], xs[:, 2 * p + 1]
            ok = np.isfinite(hi)
            if not ok.any():
                break
            # pixel i is interior iff lo <= centre(i) < hi
            a = np.ceil((lo[ok] - x0) / pitch - 0.5).astype(np.int64) - i0
            b = np.ceil((hi[ok] - x0) / pitch - 0.5).astype(np.int64) - i0
            a = np.clip(a, 0, w)
            b = np.clip(b, 0, w)
            np.add.at(diff, (rows[ok], a), 1)
            np.add.at(diff, (rows[ok], b), -1)
        img[j0:j1, i0:i1] |= np.cumsum(diff, axis=1)[:, :w] > 0


def copper_primitives(gbr, arc_steps: int = 48):
    """Exact copper geometry of a parsed gerber, as primitives."""
    prims = []
    skipped = 0

    for s in gbr.segments:
        ap = gbr.apertures.get(s.aperture)
        if ap is None or ap.function in NON_COPPER_FUNCTIONS:
            skipped += 1
            continue
        r = _stroke_radius(ap)
        pts = s.sample(arc_steps) if s.kind == "arc" else [(s.x0, s.y0),
                                                           (s.x1, s.y1)]
        for (ax, ay), (bx, by) in zip(pts[:-1], pts[1:]):
            prims.append(Capsule(ax, ay, bx, by, r))

    for f in gbr.flashes:
        ap = gbr.apertures.get(f.aperture)
        if ap is None or ap.function in NON_COPPER_FUNCTIONS:
            skipped += 1
            continue
        prims.extend(_flash_primitives(f, ap))

    for reg in gbr.regions:
        if len(reg) >= 3:
            prims.append(Poly(reg))

    return prims, skipped


def _stroke_radius(ap) -> float:
    """Half-width of a stroking aperture.  Only round/rect apertures stroke."""
    if ap.template in ("C", "O", "R") and ap.params:
        return ap.params[0] / 2.0
    return (ap.params[0] / 2.0) if ap.params else 0.0


def _flash_primitives(f, ap):
    """One flashed aperture -> primitives.  Honours the aperture TEMPLATE.

    The old script's `ap.size` shortcut (`max(params)`) is wrong for every
    template whose parameter list is not purely a size: R/O carry two sides,
    P carries a vertex count and a rotation, and the RotRect macro carries a
    rotation angle in degrees that dwarfs the two real dimensions.
    """
    p = ap.params
    t = ap.template

    if t == "C":
        return [Capsule(f.x, f.y, f.x, f.y, p[0] / 2.0)]
    if t == "R" and len(p) >= 2:
        return [Rect(f.x, f.y, p[0], p[1])]
    if t == "O" and len(p) >= 2:                    # obround = capsule
        w, h = p[0], p[1]
        if w >= h:
            return [Capsule(f.x - (w - h) / 2, f.y, f.x + (w - h) / 2, f.y, h / 2)]
        return [Capsule(f.x, f.y - (h - w) / 2, f.x, f.y + (h - w) / 2, w / 2)]
    if t == "P" and len(p) >= 2:                    # regular polygon
        n = int(p[1])
        rot = math.radians(p[2]) if len(p) >= 3 else 0.0
        r = p[0] / 2.0
        return [Poly([(f.x + r * math.cos(rot + 2 * math.pi * k / n),
                       f.y + r * math.sin(rot + 2 * math.pi * k / n))
                      for k in range(n)])]
    if t == "RotRect" and len(p) >= 3:              # KiCad macro: L x W x angle
        return [Rect(f.x, f.y, p[0], p[1], p[2])]
    # Unknown macro: fall back to a disc on the FIRST parameter, which is the
    # leading dimension in every KiCad macro, rather than on max(params).
    return [Capsule(f.x, f.y, f.x, f.y, (p[0] / 2.0) if p else 0.0)]


# ─────────────────────────────────────────────────────────────────────────────
# Monte-Carlo coverage over a zone
# ─────────────────────────────────────────────────────────────────────────────

class CopperRaster:
    """The union of the copper primitives, sampled on a regular grid.

    Each primitive is stamped into a shared boolean image, so overlapping
    traces, via pads and pours are unioned rather than double-counted — which
    is the whole reason the old script needed shapely.

    The image is sampled at PIXEL CENTRES.  That is the midpoint rule, and it
    is what keeps this honest: a pixel counts as copper iff its centre is
    inside a feature, so a feature narrower than the pitch is neither
    guaranteed to be caught nor guaranteed to be missed, and the error
    averages out over a zone rather than accumulating.  (MX17's rasterizer bug
    was the opposite convention -- endpoint-inclusive sampling, which biases
    every feature up by one pixel row and column.)  At the default 50 um
    pitch, the narrowest real feature on either layer is a 125 um trace, i.e.
    2.5 pixels wide.
    """

    # Raster window: the board's bounding box, nothing more.  The gerbers also
    # carry the 610 x 610 mm fabrication panel and its routing channel; those
    # features are real copper but they are cut away before the detector
    # exists, and rastering them would triple the image for no answer.
    XLO, XHI = 0.0, BOARD["r_out"] + 2.0
    YLO, YHI = -BOARD["edge_off"] - 2.0, BOARD["top_cut"] + 2.0
    CHUNK_PX = 4_000_000          # cap on one primitive's working set

    def __init__(self, prims, pitch=0.05):
        self.pitch = pitch
        self.x0, self.y0 = self.XLO, self.YLO
        self.nx = int(math.ceil((self.XHI - self.XLO) / pitch))
        self.ny = int(math.ceil((self.YHI - self.YLO) / pitch))
        self.img = np.zeros((self.ny, self.nx), dtype=bool)
        self.n_marked = 0

        for pr in prims:
            bx0, by0, bx1, by1 = pr.bbox
            i0 = max(int((bx0 - self.x0) / pitch) - 1, 0)
            i1 = min(int((bx1 - self.x0) / pitch) + 2, self.nx)
            j0 = max(int((by0 - self.y0) / pitch) - 1, 0)
            j1 = min(int((by1 - self.y0) / pitch) + 2, self.ny)
            if i1 <= i0 or j1 <= j0:
                continue                       # entirely off the board window
            self.n_marked += 1

            if hasattr(pr, "mark"):            # polygons fill by scanline
                pr.mark(self.img, self.x0, self.y0, pitch, i0, i1, j0, j1)
                continue

            # Capsules and rectangles evaluate a closed-form inside test on
            # the sub-grid, in row chunks so a long diagonal trace's bounding
            # box cannot blow up the working set.
            xs = self.x0 + (np.arange(i0, i1) + 0.5) * pitch
            step = max(1, self.CHUNK_PX // max(i1 - i0, 1))
            for jj in range(j0, j1, step):
                jn = min(jj + step, j1)
                ys = self.y0 + (np.arange(jj, jn) + 0.5) * pitch
                gx, gy = np.meshgrid(xs, ys)
                self.img[jj:jn, i0:i1] |= pr.contains(gx, gy)

    def rows(self, block=512):
        """Yield (pixel x, pixel y, copper) for horizontal blocks of the image."""
        xs = self.x0 + (np.arange(self.nx) + 0.5) * self.pitch
        for j0 in range(0, self.ny, block):
            j1 = min(j0 + block, self.ny)
            ys = self.y0 + (np.arange(j0, j1) + 0.5) * self.pitch
            gx, gy = np.meshgrid(xs, ys)
            yield gx, gy, self.img[j0:j1]


def board_mask(x, y):
    """Exact analytic form of P2::WedgeOutline(BOARD) -- see p2_model.py.

    The board is the 60 deg annular sector with BOTH radial edges pushed
    `edge_off` mm outward, truncated by the y = top_cut chord.  Written as
    half-planes it is cheap enough to evaluate per pixel, which a polygon
    containment test on 10^8 points is not.
    """
    phi1 = math.radians(SECTOR_DEG_)
    r2 = x * x + y * y
    n_x, n_y = -math.sin(phi1), math.cos(phi1)      # outward normal, 60 deg edge
    return ((r2 >= BOARD["r_in"] ** 2) & (r2 <= BOARD["r_out"] ** 2)
            & (y >= -BOARD["edge_off"])
            & (n_x * x + n_y * y <= BOARD["edge_off"])
            & (y <= BOARD["top_cut"]))


# ─────────────────────────────────────────────────────────────────────────────
# The pad grid
# ─────────────────────────────────────────────────────────────────────────────

def _poly_area(p):
    x, y = p[:, 0], p[:, 1]
    return abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2.0


def pad_rings(gbr):
    """Fit the polar pad grid to the F.Cu pad polygons.

    Returns (rings, stats).  Each ring is a dict with the exact gerber values:
      k, n, r_in, r_out, dphi_pitch, phi_first (centre of the first pad),
      dphi_pad, area.
    """
    recs = []
    for reg in gbr.regions:
        p = np.asarray(reg, dtype=float)
        a = _poly_area(p)
        if not (PAD_AREA_LO < a < PAD_AREA_HI):
            continue
        r = np.hypot(p[:, 0], p[:, 1])
        ph = np.arctan2(p[:, 1], p[:, 0])
        recs.append((r.min(), r.max(), ph.min(), ph.max(), a))
    recs = np.array(recs)
    if len(recs) == 0:
        raise SystemExit("no pad polygons found")

    r_in0 = recs[:, 0].min()
    pitch0 = 11.43                       # seed; refined below
    k = np.rint((recs[:, 0] - r_in0) / pitch0).astype(int)

    rings, resid = [], []
    for kk in sorted(set(k)):
        s = k == kk
        phc = np.sort((recs[s, 2] + recs[s, 3]) / 2.0)
        n = int(s.sum())
        d = np.diff(phc)
        rings.append(dict(
            k=int(kk), n=n,
            r_in=float(recs[s, 0].mean()), r_out=float(recs[s, 1].mean()),
            dphi_pitch=float(np.median(d)) if n > 1 else float("nan"),
            phi_first=float(phc[0]),
            dphi_pad=float((recs[s, 3] - recs[s, 2]).mean()),
            area=float(recs[s, 4].mean()),
        ))
        if n > 2:
            resid.append(float(np.ptp(d)))          # phi-pitch spread in ring

    r_ins = np.array([r["r_in"] for r in rings])
    fit = np.polyfit(np.arange(len(rings)), r_ins, 1)
    stats = dict(
        n_pads=int(len(recs)),
        n_rings=len(rings),
        pitch=float(fit[0]),
        r_in0=float(fit[1]),
        pitch_resid=float(np.abs(r_ins - np.polyval(fit, np.arange(len(rings)))).max()),
        pad_dr=float((recs[:, 1] - recs[:, 0]).mean()),
        pad_dr_spread=float((recs[:, 1] - recs[:, 0]).std()),
        phi_pitch_spread=float(max(resid)) if resid else 0.0,
        cu_area=float(sum(r["n"] * r["area"] for r in rings)),
        # A pad polygon vs the perfect annular sector it is meant to be
        rect_fill=float(np.mean(recs[:, 4] /
                        (((recs[:, 0] + recs[:, 1]) / 2) *
                         (recs[:, 3] - recs[:, 2]) *
                         (recs[:, 1] - recs[:, 0])))),
    )
    return rings, stats


def ring_zone_area(rings):
    """Area of the pad-field envelope: each ring's full replicated wedge."""
    return sum(r["n"] * r["dphi_pitch"] * (r["r_out"] ** 2 - r["r_in"] ** 2) / 2
               for r in rings)


# ─────────────────────────────────────────────────────────────────────────────
# Report
# ─────────────────────────────────────────────────────────────────────────────

def zone_edges(rings, n_pad_bands=8):
    """Radial band edges for the homogenized copper zones.

    Aligned to the pad field rather than to round numbers: the first band is
    the inner margin the artwork never reaches, the last is the fan-out, and
    the pad field is split into `n_pad_bands` so B.Cu's routing-density
    gradient is carried by the model instead of averaged away.
    """
    r0, r1 = rings[0]["r_in"], rings[-1]["r_out"]
    return np.concatenate([[BOARD["r_in"]],
                           np.linspace(r0, r1, n_pad_bands + 1),
                           [BOARD["r_out"]]])


def measure(layer, rings, pitch, verbose=True):
    """Coverage of one layer, over the board and over the three readout zones.

    Returns (zones, profile).  `zones` maps label -> (coverage, area_cm2);
    `profile` is the radial coverage gradient in `nbins` equal-radius bands.
    """
    gbr = parse(os.path.join(GDIR, f"P2_BASKET-{layer}.gbr"))
    prims, skipped = copper_primitives(gbr)
    ras = CopperRaster(prims, pitch=pitch)

    r_pad_in, r_pad_out = rings[0]["r_in"], rings[-1]["r_out"]
    zones = [
        (f"inner margin  r < {r_pad_in:.1f}", BOARD["r_in"], r_pad_in),
        (f"pad field     {r_pad_in:.1f}-{r_pad_out:.1f}", r_pad_in, r_pad_out),
        (f"fan-out       r > {r_pad_out:.1f}", r_pad_out, BOARD["r_out"]),
        ("whole board outline", BOARD["r_in"], BOARD["r_out"]),
    ]
    edges = zone_edges(rings)
    nbins = len(edges) - 1

    if verbose:
        print(f"\n{layer}: {len(gbr.segments)} segments, {len(gbr.flashes)} "
              f"flashes, {len(gbr.regions)} regions -> {len(prims)} primitives"
              f"  ({skipped} non-copper apertures skipped)")
        print(f"    raster {ras.nx} x {ras.ny} px at {pitch*1000:.0f} um; "
              f"{ras.n_marked} primitives inside the board window")

    n_tot = np.zeros(len(zones), dtype=np.int64)
    n_cu = np.zeros(len(zones), dtype=np.int64)
    p_tot = np.zeros(nbins, dtype=np.int64)
    p_cu = np.zeros(nbins, dtype=np.int64)

    for gx, gy, cu in ras.rows():
        on = board_mask(gx, gy)
        if not on.any():
            continue
        r = np.hypot(gx[on], gy[on])
        c = cu[on]
        for i, (_, r0, r1) in enumerate(zones):
            sel = (r >= r0) & (r < r1)
            n_tot[i] += int(sel.sum())
            n_cu[i] += int(c[sel].sum())
        b = np.clip(np.searchsorted(edges, r, side="right") - 1, 0, nbins - 1)
        p_tot += np.bincount(b, minlength=nbins)
        p_cu += np.bincount(b, weights=c, minlength=nbins).astype(np.int64)

    px_cm2 = pitch * pitch / 100.0
    out = {}
    for i, (label, _, _) in enumerate(zones):
        f = n_cu[i] / n_tot[i] if n_tot[i] else float("nan")
        out[label] = (f, n_tot[i] * px_cm2, n_cu[i] * px_cm2)
        if verbose:
            print(f"    {label:28s}  coverage {f:.4f}   "
                  f"({n_cu[i]*px_cm2:8.2f} of {n_tot[i]*px_cm2:8.2f} cm2)")

    profile = [(edges[i], edges[i + 1],
                p_cu[i] / p_tot[i] if p_tot[i] else float("nan"))
               for i in range(nbins)]
    return out, profile


HEADER_TEMPLATE = '''#pragma once
// P2PadMap.hh  --  GENERATED, DO NOT EDIT BY HAND
//
// Regenerate with:
//     python scripts/gerber/extract_readout_pattern.py --write
//
// The P2 F.Cu readout pad field, measured feature-by-feature out of
//     design/gerbers/P2_BASKET_Apr26/Gerber/P2_BASKET-F_Cu.gbr
// (KiCad 9.0.2, 2026-01-20).  {n_pads} pad polygons, no exceptions.
//
// The artwork is an exactly regular POLAR grid, which is what makes it
// cheap to build as real geometry:
//
//   * {n_rings} rings on a constant radial pitch of {pitch:.6f} mm
//     (max deviation of any ring from the fitted line: {pitch_resid:.4f} mm)
//   * every pad is {pad_dr:.5f} mm tall radially (spread {pad_dr_spread:.5f} mm),
//     leaving a {rad_gap:.5f} mm radial gap between rings
//   * within a ring the pads are uniformly spaced in phi to
//     {phi_pitch_spread:.2e} rad -- i.e. exactly, at the 1e-6 mm coordinate
//     resolution of the gerber.  Each ring is therefore ONE G4PVReplica.
//   * each pad is a true annular sector: its polygon area is {rect_fill:.4f} of
//     r_c * dphi * dr, so a G4Tubs segment is not an approximation.
//
// Pads per ring grows with radius to keep the cell area constant
// ("equalized" pads): {n_first} on the inner ring, {n_last} on the outer,
// {n_pads} channels total = 10 connectors x 128.
//
// Angles are RADIANS, lengths are mm.  phiFirst is the CENTRE of the first
// pad of the ring; the replica envelope therefore starts at
// phiFirst - dPhiPitch/2.

namespace P2 {{

struct PadRing {{
    int    n;           // pads in this ring
    double rIn;         // pad inner radius  [mm]
    double rOut;        // pad outer radius  [mm]
    double dPhiPitch;   // pad-to-pad angular pitch [rad]
    double phiFirst;    // centre of the first pad  [rad]
    double dPhiPad;     // angular width of the copper pad itself [rad]
}};

// Total copper in the pad field: {cu_area_cm2:.2f} cm^2 over a {zone_area_cm2:.2f} cm^2
// replica envelope, i.e. {pad_cov:.4f} of it (the azimuthal gaps).  Folding in
// the radial gaps between rings gives {pad_cov_pitch:.4f} of the full pad-field
// annulus -- that is the number a homogenized F.Cu sheet should carry there.
constexpr int kNPadRings = {n_rings};
constexpr int kNPads     = {n_pads};

constexpr PadRing kPadRings[kNPadRings] = {{
{rows}
}};

// ── Homogenized copper, by radial band ──────────────────────────────────────
//
// Everything that is NOT the pad field is irregular artwork and stays
// homogenized -- but per band, not as one board-wide average.  Coverage is
// the copper area fraction of the band, measured on the exact vector
// geometry (see the script; rasterized at pixel centres).
//
// Three facts these numbers encode that a single average destroys:
//   * Inside r = {r_pad_in:.1f} mm there is no signal copper at all on either
//     layer -- only a ~1.7 mm board-edge guard band and four mounting pads,
//     which is why the first band reads the SAME on F.Cu and B.Cu (the two
//     layers carry identical artwork there).  A board-wide F.Cu sheet instead
//     puts ~98 % copper across that 20 mm annulus.
//   * The fan-out band beyond r = {r_pad_out:.1f} mm is only ~18 % copper on
//     F.Cu: this is a 2-layer board, every pad drops through two vias
//     immediately, and ALL the routing happens on B.Cu.  A board-wide sheet
//     gets this zone wrong by a factor of five.
//   * B.Cu is not a ground plane -- it is {n_bcu_traces} stroked 0.125 mm
//     signal traces, and their density climbs steadily with radius as more
//     channels are gathered toward the connectors.
//
// fCu is used only when the pattern is switched off (--homogenized-readout);
// with the pattern on, the bands inside the pad field are replaced by the
// real pads and only the first and last bands are built from this table.

struct CuBand {{
    double rIn;    // [mm]
    double rOut;   // [mm]
    double fCu;    // F.Cu area coverage in this band
    double bCu;    // B.Cu area coverage in this band
}};

constexpr int kNCuBands = {n_bands};
constexpr CuBand kCuBands[kNCuBands] = {{
{band_rows}
}};

// Pad-field boundaries, so the geometry code does not have to re-derive them.
constexpr double kPadFieldRIn  = {r_pad_in:.6f};
constexpr double kPadFieldROut = {r_pad_out:.6f};

}}  // namespace P2
'''


def write_header(rings, stats, zone_area, bands):
    rows = "\n".join(
        "    {{ {n:2d}, {r_in:11.6f}, {r_out:11.6f}, {dphi_pitch:.9f}, "
        "{phi_first:.9f}, {dphi_pad:.9f} }},".format(**r) for r in rings)
    band_rows = "\n".join(
        f"    {{ {r0:10.6f}, {r1:10.6f}, {ff:.6f}, {fb:.6f} }},"
        for r0, r1, ff, fb in bands)
    pad_cov = stats["cu_area"] / zone_area
    txt = HEADER_TEMPLATE.format(
        rows=rows, band_rows=band_rows,
        n_pads=stats["n_pads"], n_rings=stats["n_rings"],
        pitch=stats["pitch"], pitch_resid=stats["pitch_resid"],
        pad_dr=stats["pad_dr"], pad_dr_spread=stats["pad_dr_spread"],
        rad_gap=stats["pitch"] - stats["pad_dr"],
        phi_pitch_spread=stats["phi_pitch_spread"],
        rect_fill=stats["rect_fill"],
        n_first=rings[0]["n"], n_last=rings[-1]["n"],
        cu_area_cm2=stats["cu_area"] / 100.0,
        zone_area_cm2=zone_area / 100.0,
        pad_cov=pad_cov,
        pad_cov_pitch=pad_cov * stats["pad_dr"] / stats["pitch"],
        n_bands=len(bands),
        n_bcu_traces=stats["n_bcu_traces"],
        r_pad_in=rings[0]["r_in"], r_pad_out=rings[-1]["r_out"],
    )
    with open(HEADER, "w") as fh:
        fh.write(txt)
    print(f"\nWrote {os.path.relpath(HEADER, REPO)}  "
          f"({stats['n_rings']} rings, {stats['n_pads']} pads, "
          f"{len(bands)} copper bands)")


def compare_old(rings):
    """Reproduce the aperture bug in analyze_cu_coverage.py, to size it."""
    gbr = parse(os.path.join(GDIR, "P2_BASKET-F_Cu.gbr"))
    bad = 0
    worst = 0.0
    for f in gbr.flashes:
        ap = gbr.apertures.get(f.aperture)
        if ap is None:
            continue
        size = max(ap.params) if ap.params else 0.0     # the old `ap.size`
        real = _flash_primitives(f, ap)
        rb = real[0].bbox
        if size > 2 * max(rb[2] - rb[0], rb[3] - rb[1]) + 1e-6:
            bad += 1
            worst = max(worst, size)
    print("\n--- analyze_cu_coverage.py aperture bug ---")
    print(f"  flashes whose old buffered radius is oversized : {bad}")
    print(f"  worst old 'aperture size'                      : {worst:.3f} mm"
          f"   (a {worst/2:.1f} mm-radius disc per flash)")
    print("  cause: ApertureDef.size = max(params); RotRect's 3rd parameter")
    print("         is a ROTATION ANGLE in degrees, not a length.")
    print("  effect: the 10 connector footprints get painted over as ~150 mm")
    print("         discs, which is what lifts F_Cu to the 0.937 board /")
    print("         0.983 active figures quoted in SimConfig and the docs.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true",
                    help="write include/P2PadMap.hh")
    ap.add_argument("--compare", action="store_true",
                    help="quantify the analyze_cu_coverage.py aperture bug")
    ap.add_argument("--pitch", type=float, default=0.05,
                    help="raster pitch in mm (default 0.05; the narrowest "
                         "real feature is a 0.125 mm trace)")
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    fcu = parse(os.path.join(GDIR, "P2_BASKET-F_Cu.gbr"))
    rings, stats = pad_rings(fcu)
    zone_area = ring_zone_area(rings)
    bcu = parse(os.path.join(GDIR, "P2_BASKET-B_Cu.gbr"))
    stats["n_bcu_traces"] = sum(
        1 for s in bcu.segments
        if bcu.apertures.get(s.aperture) is not None
        and bcu.apertures[s.aperture].function not in NON_COPPER_FUNCTIONS)
    del bcu

    print("=" * 74)
    print("P2 READOUT PAD PATTERN  --  F_Cu, exact vector geometry")
    print("=" * 74)
    print(f"  pad polygons        : {stats['n_pads']}")
    print(f"  rings               : {stats['n_rings']}"
          f"   ({rings[0]['n']} pads inner .. {rings[-1]['n']} outer)")
    print(f"  radial pitch        : {stats['pitch']:.6f} mm"
          f"   (max residual {stats['pitch_resid']:.4f} mm)")
    print(f"  pad radial height   : {stats['pad_dr']:.5f} mm"
          f"   (spread {stats['pad_dr_spread']:.5f} mm)")
    print(f"  radial gap          : {stats['pitch'] - stats['pad_dr']:.5f} mm")
    print(f"  phi pitch spread    : {stats['phi_pitch_spread']:.2e} rad"
          f"   <- uniform => one G4PVReplica per ring")
    print(f"  pad vs annular sector: area / (r*dphi*dr) = {stats['rect_fill']:.5f}"
          f"   <- G4Tubs segment is exact")
    print(f"  radial span         : {rings[0]['r_in']:.3f} .. "
          f"{rings[-1]['r_out']:.3f} mm")
    print(f"  copper in pad field : {stats['cu_area']/100:.2f} cm2 over "
          f"{zone_area/100:.2f} cm2 envelope = {stats['cu_area']/zone_area:.4f}")

    print(f"\n  {'k':>3} {'n':>4} {'r_in':>10} {'r_out':>10} "
          f"{'dphi_pitch':>12} {'phi_first':>11} {'dphi_pad':>11} "
          f"{'az gap':>8}")
    for r in rings:
        rc = (r["r_in"] + r["r_out"]) / 2
        print(f"  {r['k']:>3} {r['n']:>4} {r['r_in']:>10.4f} {r['r_out']:>10.4f} "
              f"{r['dphi_pitch']:>12.8f} {r['phi_first']:>11.8f} "
              f"{r['dphi_pad']:>11.8f} "
              f"{(r['dphi_pitch']-r['dphi_pad'])*rc:>8.4f}")

    print("\n" + "=" * 74)
    print("COPPER COVERAGE BY ZONE  --  exact vector geometry, rasterized")
    print("=" * 74)
    profiles = {}
    for layer in ("F_Cu", "B_Cu"):
        _, profiles[layer] = measure(layer, rings, args.pitch)

    print("\n  radial coverage bands (what the homogenized zones are built from):")
    print(f"    {'r [mm]':>18} {'F_Cu':>8} {'B_Cu':>8}")
    bands = []
    for (r0, r1, ff), (_, _, fb) in zip(profiles["F_Cu"], profiles["B_Cu"]):
        bands.append((r0, r1, ff, fb))
        print(f"    {r0:7.3f} .. {r1:7.3f} {ff:>8.4f} {fb:>8.4f}")

    if args.compare:
        compare_old(rings)

    if args.write:
        write_header(rings, stats, zone_area, bands)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge
        fig, ax = plt.subplots(figsize=(9, 9))
        for r in rings:
            for i in range(r["n"]):
                pc = r["phi_first"] + i * r["dphi_pitch"]
                ax.add_patch(Wedge((0, 0), r["r_out"],
                                   math.degrees(pc - r["dphi_pad"] / 2),
                                   math.degrees(pc + r["dphi_pad"] / 2),
                                   width=r["r_out"] - r["r_in"],
                                   facecolor="#b5651d", edgecolor="none"))
        b = np.asarray(wedge_outline(**BOARD))
        ax.plot(b[:, 0], b[:, 1], "k-", lw=0.8)
        ax.set_aspect("equal")
        ax.set_xlim(0, 680)
        ax.set_ylim(-60, 580)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_title("P2 F.Cu pad copper as built in Geant4 "
                     f"({stats['n_pads']} pads, {stats['n_rings']} rings)")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"\nWrote {args.plot}")


if __name__ == "__main__":
    main()
