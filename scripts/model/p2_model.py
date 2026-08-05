#!/usr/bin/env python3
"""
Python mirror of the P2 Geant4 geometry (include/P2Wedge.hh +
DetectorConstruction::ConstructP2). Used by plot_p2_model.py to draw the
as-built model without needing Geant4.

KEEP IN SYNC with include/P2Wedge.hh / SimConfig.hh — both carry the same
parameter table; docs/P2_MODEL.md documents which numbers are measured and
which are guesses.

Coordinates: gerber frame, apex (= MESA beam axis) at (0,0), mm everywhere.
z = 0 at the front window plane (top of the gas frame), +z downstream.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

# ── Board outline (Edge_Cuts gerber, exact) ─────────────────────────────────
BOARD = dict(r_in=95.0, r_out=650.0, edge_off=10.0, top_cut=540.01)
SECTOR_DEG = 60.0

# ── Readout PCB stack (Stack_Up_P2.txt, exact) ──────────────────────────────
T_CU_F, T_FR4, T_CU_B = 0.018, 0.200, 0.018

# ── Gas frame (P2_Frame_V1_3mm.stp bbox analysis; opening ≈, see docs) ──────
FRAME = dict(r_in=104.5, r_out=605.5, edge_off=7.5, top_cut=537.5)
OPENING = dict(r_in=107.0, r_out=603.0, edge_off=0.0, top_cut=535.0)

# ── Active area (bulk mask, exact; informational) ───────────────────────────
ACTIVE = dict(r_in=119.87, r_out=589.80)

# ── Defaults from SimConfig.hh (GUESS flags in docs/P2_MODEL.md) ────────────
T_DRIFT      = 3.0      # confirmed baseline; campaign scans 1..4 mm
T_AMP        = 0.150    # confirmed 2026-08-05
MESH_WIRE    = 0.019    # woven SS mesh: wire diameter — confirmed 2026-08-05
MESH_OPEN    = 0.048    # opening ("48x19"); pitch = wire + opening = 67 um
T_MESH       = 2 * MESH_WIRE          # weave height, effective-density slab
MESH_FILL    = math.pi * MESH_WIRE / (4 * (MESH_WIRE + MESH_OPEN))  # ~0.22
FRONT_GAP    = 4.0      # window -> first drift foil (window->mesh still 8 mm,
                        # frame STEP: ledge z=3 -> top z=8)
CATH_GAP     = 1.0      # between the two drift-cathode foils ("maybe 1 mm")
BACK_GAP     = 1.0      # carbon back-frame depth; 1 mm normal, <=3 this prod.
BULGE_FRONT  = 10.0     # GUESS
BULGE_BACK   = 5.0      # GUESS
T_WINDOW     = 0.040    # MX17-like
T_CATH_MYLAR = 0.012    # each foil — thickness GUESS
T_CATH_AL    = 0.0001   # on the drift-gas side of the downstream foil
FCU_COVERAGE = 0.983    # Cu area fraction over the active area, from gerbers
BCU_COVERAGE = 0.174    #   (B.Cu = signal lines; radial 0.03->0.26);
                        #   scripts/gerber/analyze_cu_coverage.py

N_TERRACE = 6


def wedge_outline(r_in: float, r_out: float, edge_off: float, top_cut: float,
                  n_arc: int = 96) -> list[tuple[float, float]]:
    """Port of P2::WedgeOutline (itself a port of p2_wedge_model.outline)."""
    phi1 = math.radians(SECTOR_DEG)
    a_in_lo = -math.asin(edge_off / r_in)
    a_in_hi = phi1 - a_in_lo
    a_out_lo = -math.asin(edge_off / r_out)
    a_out_edge = phi1 + math.asin(edge_off / r_out)
    a_out_cut = math.asin(top_cut / r_out)

    # The y = top_cut chord only clips the outer arc if it gets there before
    # the 60 deg radial edge does. True for the board profile, false for the
    # gas frame and the frame opening — forcing a chord corner in anyway
    # folded those polygons over themselves and killed G4ExtrudedSolid.
    # Fixed 2026-08-05; keep in sync with src/P2Wedge.cc.
    clipped = a_out_cut < a_out_edge
    a_out_hi = a_out_cut if clipped else a_out_edge

    pts: list[tuple[float, float]] = []
    arc = lambda r, a0, a1: pts.extend(
        (r * math.cos(a0 + (a1 - a0) * i / n_arc),
         r * math.sin(a0 + (a1 - a0) * i / n_arc)) for i in range(n_arc + 1))

    pts.append((r_in * math.cos(a_in_lo), -edge_off))
    pts.append((r_out * math.cos(a_out_lo), -edge_off))
    arc(r_out, a_out_lo, a_out_hi)
    if clipped:
        x_top = (top_cut * math.cos(phi1) - edge_off) / math.sin(phi1)
        pts.append((x_top, top_cut))
    pts.append((r_in * math.cos(a_in_hi), r_in * math.sin(a_in_hi)))
    arc(r_in, a_in_hi, a_in_lo)

    ded = []
    for p in pts:
        if ded and (p[0]-ded[-1][0])**2 + (p[1]-ded[-1][1])**2 < 1e-12:
            continue
        ded.append(p)
    if len(ded) > 1 and (ded[0][0]-ded[-1][0])**2 + (ded[0][1]-ded[-1][1])**2 < 1e-12:
        ded.pop()
    return ded


def centroid(poly) -> tuple[float, float]:
    a2 = cx = cy = 0.0
    for i in range(len(poly)):
        (x0, y0), (x1, y1) = poly[i], poly[(i + 1) % len(poly)]
        w = x0 * y1 - x1 * y0
        a2 += w
        cx += (x0 + x1) * w
        cy += (y0 + y1) * w
    return cx / (3 * a2), cy / (3 * a2)


def scale_about(poly, c, s):
    return [(c[0] + s * (x - c[0]), c[1] + s * (y - c[1])) for x, y in poly]


@dataclass
class Layer:
    name: str
    z0: float               # upstream face
    t: float
    material: str
    poly: list = field(repr=False, default=None)
    color: str = "0.5"
    alpha: float = 0.8


def terrace_profile(H: float, n: int = N_TERRACE):
    """(sigma_k, h_k) k=0..n of the spherical-cap terraced dome, as in C++."""
    sig = [max(math.cos(k * math.pi / (2 * n)), 0.04) for k in range(n + 1)]
    h = [H * math.sin(k * math.pi / (2 * n)) for k in range(n + 1)]
    return sig, h


def build_window(tag: str, z_base: float, sign: int, H: float,
                 opening_poly, open_c) -> list[Layer]:
    """Terraced dome: gas steps + mylar rings/cap, mirroring BuildWindow()."""
    sig, h = terrace_profile(H)
    out = []
    for k in range(N_TERRACE):
        poly = scale_about(opening_poly, open_c, sig[k])
        dz = h[k + 1] - h[k]
        z0 = z_base + sign * h[k] - (dz if sign < 0 else 0)
        out.append(Layer(f"{tag}_Gas{k}", z0, dz, "gas", poly,
                         color="#8dd8f0", alpha=0.15))
        zm = z_base + sign * h[k + 1] - (T_WINDOW if sign < 0 else 0)
        if k < N_TERRACE - 1:
            inner = scale_about(opening_poly, open_c, sig[k + 1])
            lay = Layer(f"{tag}_Mylar{k}", zm, T_WINDOW, "mylar", poly,
                        color="#63b8d8", alpha=0.9)
            lay.hole = inner  # annular ring
        else:
            lay = Layer(f"{tag}_Mylar{k}", zm, T_WINDOW, "mylar", poly,
                        color="#63b8d8", alpha=0.9)
            lay.hole = None
        out.append(lay)
    return out


def build_stack():
    """Full detector: returns (layers, frames, windows, key_z) — all Layer."""
    board = wedge_outline(**BOARD)
    frame = wedge_outline(**FRAME)
    opening = wedge_outline(**OPENING)
    open_c = centroid(opening)

    z = 0.0
    layers = []

    def add(name, t, mat, poly, color, alpha=0.85):
        nonlocal z
        layers.append(Layer(name, z, t, mat, poly, color, alpha))
        z += t

    add("FrontGas",            FRONT_GAP,    "gas",   opening, "#8dd8f0", 0.15)
    add("DriftCathode_Mylar1", T_CATH_MYLAR, "mylar", opening, "#a8e0a8")
    add("DriftCathode_Gas",    CATH_GAP,     "gas",   opening, "#8dd8f0", 0.15)
    add("DriftCathode_Mylar2", T_CATH_MYLAR, "mylar", opening, "#a8e0a8")
    add("DriftCathode_Al",     T_CATH_AL,    "Al",    opening, "#b0b0b0")
    add("DriftGas",            T_DRIFT,      "gas",   opening, "#3380ff", 0.30)
    add("Micromesh",           T_MESH,       "steel(eff)", opening, "#808080")
    add("AmpGas",              T_AMP,        "gas",   opening, "#ff4d4d", 0.35)
    add("PCB_Cu_F",            T_CU_F,       "Cu(x0.98)", board, "#cc6619")
    add("PCB_FR4",             T_FR4,        "FR4",   board,   "#339933")
    add("PCB_Cu_B",            T_CU_B,       "Cu(x0.17)", board, "#cc6619")
    frame_h = (FRONT_GAP + T_CATH_MYLAR + CATH_GAP + T_CATH_MYLAR
               + T_CATH_AL + T_DRIFT + T_MESH + T_AMP)
    z_pcb_end = z
    add("BackGas",            BACK_GAP,     "gas",   opening, "#8dd8f0", 0.15)
    z_back_win = z

    frames = [Layer("GasFrame", 0.0, frame_h, "plastic", frame,
                    "#e6e6d9", 0.9),
              Layer("GasFrameBack", z_pcb_end, BACK_GAP, "carbon", frame,
                    "#4d4d4d", 0.9)]
    for f in frames:
        f.hole = wedge_outline(**OPENING)

    windows = (build_window("FrontWindow", 0.0, -1, BULGE_FRONT, opening, open_c)
               + build_window("BackWindow", z_back_win, +1, BULGE_BACK, opening, open_c))

    key_z = dict(front_window=0.0, pad_plane=frame_h, pcb_end=z_pcb_end,
                 back_window=z_back_win,
                 z_min=-(BULGE_FRONT + T_WINDOW),
                 z_max=z_back_win + BULGE_BACK + T_WINDOW)
    return layers, frames, windows, key_z


if __name__ == "__main__":
    layers, frames, windows, key_z = build_stack()
    print(f"{'layer':22s} {'z0 [mm]':>10s} {'t [mm]':>9s}  material")
    for L in layers:
        print(f"{L.name:22s} {L.z0:10.4f} {L.t:9.4f}  {L.material}")
    for F in frames:
        print(f"{F.name:22s} {F.z0:10.4f} {F.t:9.4f}  {F.material} (ring)")
    print(f"\nwindows: {len(windows)} terrace volumes "
          f"(2 x {N_TERRACE} gas + {N_TERRACE} mylar)")
    for k, v in key_z.items():
        print(f"  {k:14s} z = {v:9.4f} mm")
