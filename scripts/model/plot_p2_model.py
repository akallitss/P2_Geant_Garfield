#!/usr/bin/env python3
"""
Plots of the as-built P2 Geant4 model (scripts/model/p2_model.py) and of the
real board copper (design gerbers).

    python scripts/model/plot_p2_model.py            # all figures
    python scripts/model/plot_p2_model.py --only 3d  # board | xsec | 3d

Outputs to docs/figures/ by default:
    p2_board_copper.png    F.Cu / B.Cu rendered from the production gerbers
    p2_stack_xsec.png      model cross-section at phi=30 (true scale + zooms)
    p2_3d_overview.png     assembled 3D views (z exaggerated where noted)
    p2_3d_exploded.png     exploded 3D view (no axes, for slides)
    p2_3d_exploded_lr.png  same, rolled 90 deg: stack reads left to right,
                           labels on leader arrows
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Circle, Rectangle, Polygon as MplPolygon

# On machines where an old distro matplotlib coexists with a pip user install,
# the distro's mpl_toolkits (a regular package) shadows the user install's
# namespace portion and breaks the 3D toolkit. Point the package path at the
# matplotlib that was actually imported, then register the 3d projection.
import mpl_toolkits
_user_mt = os.path.join(os.path.dirname(os.path.dirname(matplotlib.__file__)),
                        "mpl_toolkits")
if os.path.isdir(_user_mt) and _user_mt not in mpl_toolkits.__path__:
    mpl_toolkits.__path__.insert(0, _user_mt)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.projections import register_projection
register_projection(Axes3D)

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "scripts", "gerber"))

import p2_model as M
from gerber_outline import parse

GERBER_DIR = os.path.join(REPO, "design", "gerbers", "P2_BASKET_Apr26", "Gerber")
CU_COLOR = "#b87333"


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def poly_contains(poly, x, y) -> bool:
    inside = False
    n = len(poly)
    for i in range(n):
        (x0, y0), (x1, y1) = poly[i], poly[(i + 1) % n]
        if (y0 > y) != (y1 > y):
            xi = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if x < xi:
                inside = not inside
    return inside


def ray_extent(poly, phi_deg=30.0, r_max=700.0, dr=0.25):
    """Intervals of r where the ray at phi crosses the polygon."""
    c, s = math.cos(math.radians(phi_deg)), math.sin(math.radians(phi_deg))
    rs = np.arange(0, r_max, dr)
    inside = np.array([poly_contains(poly, r * c, r * s) for r in rs])
    ivals, start = [], None
    for r, i in zip(rs, inside):
        if i and start is None:
            start = r
        elif not i and start is not None:
            ivals.append((start, r)); start = None
    if start is not None:
        ivals.append((start, rs[-1]))
    return ivals


def resample(poly, k=72):
    """Resample closed polygon boundary to k points by arclength."""
    p = np.array(poly + [poly[0]])
    seg = np.hypot(*np.diff(p, axis=0).T)
    cum = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, cum[-1], k, endpoint=False)
    x = np.interp(t, cum, p[:, 0])
    y = np.interp(t, cum, p[:, 1])
    return np.column_stack([x, y])


def prism_faces(poly, z0, z1, k=72):
    """Faces (list of Nx3 arrays) of an extruded polygon for Poly3DCollection."""
    b = resample(poly, k)
    bot = np.column_stack([b, np.full(len(b), z0)])
    top = np.column_stack([b, np.full(len(b), z1)])
    faces = [bot, top]
    for i in range(len(b)):
        j = (i + 1) % len(b)
        faces.append(np.array([bot[i], bot[j], top[j], top[i]]))
    return faces


def ring_faces(outer, inner, z0, z1, k=72):
    """Faces of an extruded ring (outer minus inner), matched by arclength."""
    o, n = resample(outer, k), resample(inner, k)
    faces = []
    for i in range(k):
        j = (i + 1) % k
        for zz in (z0, z1):  # top + bottom annulus as quads
            faces.append(np.array([[*o[i], zz], [*o[j], zz],
                                   [*n[j], zz], [*n[i], zz]]))
        faces.append(np.array([[*o[i], z0], [*o[j], z0], [*o[j], z1], [*o[i], z1]]))
        faces.append(np.array([[*n[i], z0], [*n[j], z0], [*n[j], z1], [*n[i], z1]]))
    return faces


def add_prism(ax, poly, z0, z1, color, alpha, hole=None, k=72, lw=0.15):
    faces = (ring_faces(poly, hole, z0, z1, k) if hole is not None
             else prism_faces(poly, z0, z1, k))
    ax.add_collection3d(Poly3DCollection(
        faces, facecolor=color, alpha=alpha, edgecolor="k", linewidths=lw))


def model_overlays(ax, lw=1.2):
    for params, color, ls, lab in [
            (M.BOARD, "crimson", "-", "board (Edge_Cuts)"),
            (M.FRAME, "0.35", "-", "gas frame footprint"),
            (M.OPENING, "0.35", "--", "frame opening / windows")]:
        p = M.wedge_outline(**params)
        xs = [q[0] for q in p] + [p[0][0]]
        ys = [q[1] for q in p] + [p[0][1]]
        ax.plot(xs, ys, ls, color=color, lw=lw, label=lab, zorder=5)
    th = np.linspace(math.radians(0.97), math.radians(59.0), 100)
    for r in (M.ACTIVE["r_in"], M.ACTIVE["r_out"]):
        ax.plot(r * np.cos(th), r * np.sin(th), ":", color="royalblue",
                lw=1.0, zorder=5)
    ax.plot([], [], ":", color="royalblue", label="active area (bulk mask)")
    ax.plot(0, 0, "+", color="k", ms=10, mew=1.5, zorder=6)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — board copper from the production gerbers
# ─────────────────────────────────────────────────────────────────────────────

def draw_gerber(ax, gf, lw_scale):
    # Zone pours are single serpentine outlines with 10^4-10^5 vertices that
    # snake around every pad; a plain fill turns them into solid blobs, so
    # draw them as faint boundaries instead (the pads carry the geometry).
    small = [r for r in gf.regions if 2 < len(r) <= 10000]
    zones = [r for r in gf.regions if len(r) > 10000]
    polys = [MplPolygon(np.array(r), closed=True) for r in small]
    ax.add_collection(PatchCollection(polys, facecolor=CU_COLOR,
                                      edgecolor="none", alpha=0.95, zorder=2))
    if zones:
        ax.add_collection(LineCollection(
            [np.array(z) for z in zones], colors=CU_COLOR, linewidths=0.1,
            alpha=0.35, zorder=1))
    by_w = {}
    for s in gf.segments:
        w = gf.apertures.get(s.aperture)
        w = w.size if w else 0.15
        by_w.setdefault(round(w, 3), []).append(s)
    for w, segs in by_w.items():
        lines = [s.sample(16) for s in segs]
        ax.add_collection(LineCollection(
            lines, colors=CU_COLOR, linewidths=max(w * lw_scale, 0.25),
            alpha=0.9, zorder=2, capstyle="round"))
    pats = []
    for f in gf.flashes:
        ap = gf.apertures.get(f.aperture)
        if ap is None:
            continue
        if ap.template in ("R", "O") and len(ap.params) >= 2:
            pats.append(Rectangle((f.x - ap.params[0]/2, f.y - ap.params[1]/2),
                                  ap.params[0], ap.params[1]))
        elif ap.template == "C" and ap.params:
            pats.append(Circle((f.x, f.y), max(ap.params[0], 0.05) / 2))
        elif len(ap.params) >= 3:
            # macro rotated rectangle: params = [w, h, rotation_deg]
            w, h, rot = ap.params[0], ap.params[1], ap.params[2]
            pats.append(Rectangle((f.x - w/2, f.y - h/2), w, h,
                                  angle=rot, rotation_point="center"))
        else:
            pats.append(Circle((f.x, f.y), 0.25))
    ax.add_collection(PatchCollection(pats, facecolor=CU_COLOR,
                                      edgecolor="none", alpha=0.95, zorder=3))


def fig_board(outdir):
    fcu = parse(os.path.join(GERBER_DIR, "P2_BASKET-F_Cu.gbr"))
    bcu = parse(os.path.join(GERBER_DIR, "P2_BASKET-B_Cu.gbr"))

    fig = plt.figure(figsize=(19, 7.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.25, 1.0])
    axf = fig.add_subplot(gs[0]); axb = fig.add_subplot(gs[1])
    axz = fig.add_subplot(gs[2])

    for ax, gf, title in [(axf, fcu, "F.Cu — readout pads (18 µm)"),
                          (axb, bcu, "B.Cu — ground plane (18 µm)")]:
        draw_gerber(ax, gf, lw_scale=0.55)
        model_overlays(ax)
        ax.set_xlim(-30, 680); ax.set_ylim(-60, 580)
        ax.set_aspect("equal"); ax.grid(alpha=0.2, lw=0.4)
        ax.set_title(title)
        ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    axf.legend(loc="upper right", fontsize=8)

    # zoom: pads + vias + connector fan-out near the inner radius
    draw_gerber(axz, fcu, lw_scale=3.2)
    model_overlays(axz)
    axz.set_xlim(100, 210); axz.set_ylim(-5, 90)
    axz.set_aspect("equal"); axz.grid(alpha=0.2, lw=0.4)
    axz.set_title("F.Cu zoom — inner rings, vias, fan-out")
    axz.set_xlabel("x [mm]")

    fig.suptitle("P2 readout board — production gerbers (P2_BASKET_Apr26) with model overlays",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(outdir, "p2_board_copper.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — cross-section of the model at phi = 30 deg
# ─────────────────────────────────────────────────────────────────────────────

def xsec_rects(items, phi=30.0):
    """(r0, r1, z0, z1, color, alpha, name) for every model volume on the ray."""
    out = []
    for L in items:
        hole = getattr(L, "hole", None)
        ivals = ray_extent(L.poly, phi)
        if hole is not None:
            hi = ray_extent(hole, phi)
            cut = []
            for a, b in ivals:
                pieces = [(a, b)]
                for ha, hb in hi:
                    pieces = [(p, min(q, ha)) for p, q in pieces if p < min(q, ha)] + \
                             [(max(p, hb), q) for p, q in pieces if max(p, hb) < q]
                cut += pieces
            ivals = cut
        for a, b in ivals:
            out.append((a, b, L.z0, L.z0 + L.t, L.color, L.alpha, L.name))
    return out


def fig_xsec(outdir):
    layers, frames, windows, key_z = M.build_stack()
    rects = xsec_rects(layers + frames + windows)

    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(19, 6.4), width_ratios=[2.2, 1.0, 1.0])

    for ax in (ax1, ax2):
        for r0, r1, z0, z1, c, a, _ in rects:
            ax.add_patch(Rectangle((r0, z0), r1 - r0, z1 - z0,
                                   facecolor=c, alpha=max(a, 0.35),
                                   edgecolor="k", lw=0.2))
        ax.invert_yaxis()  # beam comes from the top of the plot

    ax1.annotate("", xy=(355, -1.5), xytext=(355, -11.5),
                 arrowprops=dict(arrowstyle="-|>", color="k", lw=1.6))
    ax1.text(360, -8.5, "beam (default aim r=355 mm)", fontsize=9)
    ax1.set_xlim(0, 700); ax1.set_ylim(20.5, -13.5)
    ax1.set_xlabel("r along $\\phi=30°$ [mm]"); ax1.set_ylabel("z [mm]")
    ax1.set_title("True scale — bulged windows, frame, gas volumes")
    ax1.grid(alpha=0.25, lw=0.4)

    zs = {L.name: (L.z0, L.t) for L in layers}
    z_mesh = zs["Micromesh"][0]
    z_pcb_end = zs["PCB_Cu_B"][0] + zs["PCB_Cu_B"][1]
    ax2.set_xlim(340, 372); ax2.set_ylim(z_pcb_end + 0.15, z_mesh - 0.12)
    ax2.set_xlabel("r [mm]")
    ax2.set_title("Zoom: mesh / amp gap / PCB")
    ax2.grid(alpha=0.25, lw=0.4)
    for name, lab in [("Micromesh", "mesh"), ("AmpGas", "amp gap"),
                      ("PCB_Cu_F", "F.Cu"), ("PCB_FR4", "FR4"),
                      ("PCB_Cu_B", "B.Cu")]:
        z0, t = zs[name]
        ax2.annotate(f"{lab} ({t*1000:.0f} µm)", xy=(372, z0 + t/2),
                     xytext=(373.5, z0 + t/2),
                     fontsize=8, va="center", annotation_clip=False)

    # not-to-scale stackup with thickness labels
    ax3.set_title("Stack (not to scale)")
    rows = [L for L in layers] + [f for f in frames if f.name == "GasFrame"]
    rows = layers
    y = 0
    for L in rows:
        ax3.add_patch(Rectangle((0, y), 1, 1, facecolor=L.color,
                                alpha=max(L.alpha, 0.4), edgecolor="k", lw=0.4))
        t_lab = f"{L.t*1000:.0f} µm" if L.t < 1 else f"{L.t:.0f} mm"
        ax3.text(1.05, y + 0.5, f"{L.name}   ({t_lab})", va="center", fontsize=9)
        y += 1
    ax3.text(1.05, y + 0.7,
             f"windows: 40 µm mylar,\nsag {M.BULGE_FRONT:.0f} mm front / "
             f"{M.BULGE_BACK:.0f} mm back", fontsize=9, va="bottom")
    ax3.set_xlim(0, 4.2); ax3.set_ylim(-0.5, y + 3)
    ax3.invert_yaxis(); ax3.axis("off")

    fig.suptitle("P2 wedge Geant4 model — cross-section at $\\phi=30°$ "
                 "(z = 0 at front window plane, beam along +z)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(outdir, "p2_stack_xsec.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3/4 — 3D renders
# ─────────────────────────────────────────────────────────────────────────────

EXPLODE_GROUPS = {  # layer-name prefix -> explode group index
    "FrontWindow": -2, "FrontGas": -1,
    "DriftCathode": 0, "DriftGas": 1, "Micromesh": 2, "AmpGas": 2,
    "PCB": 3, "BackGas": 4, "BackWindow": 5,
    "GasFrame": 1, "GasFrameBack": 4}

# Label -> the explode groups it points at, for the annotated exploded views.
EXPLODE_LABELS = [  # wrapped so a big font still fits a narrow column
    ("front window\n(bulged mylar)",                   (-2,)),
    ("drift cathode\n+ drift gas + frame",   (0, 1)),
    ("mesh + amp gap",                                 (2,)),
    ("readout PCB\n(Cu/FR4/Cu)",                       (3,)),
    ("back gas\n+ carbon back frame",                  (4,)),
    ("back window\n(bulged mylar)",                    (5,)),
]
LABEL_FS = 15


def grp(name):
    for k in sorted(EXPLODE_GROUPS, key=len, reverse=True):
        if name.startswith(k):
            return EXPLODE_GROUPS[k]
    return 0


def _poly_area(poly):
    p = np.asarray(poly)
    x, y = p[:, 0], p[:, 1]
    return abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2


def group_anchors(zexag, explode):
    """(group index) -> (mid-z, footprint) of that group in the exploded
    render, plus the full z range the exploded stack occupies."""
    layers, frames, windows, _ = M.build_stack()
    zs: dict[int, list[float]] = {}
    poly: dict[int, list] = {}
    lo, hi = 1e9, -1e9
    for L in layers + frames + windows:
        g = grp(L.name)
        zs.setdefault(g, []).append(L.z0 + L.t / 2)
        if g not in poly or _poly_area(L.poly) > _poly_area(poly[g]):
            poly[g] = L.poly  # the group's widest piece = its silhouette
        lo = min(lo, L.z0 * zexag + g * explode)
        hi = max(hi, (L.z0 + L.t) * zexag + g * explode)
    anchors = {g: ((sum(v) / len(v)) * zexag + g * explode, poly[g])
               for g, v in zs.items()}
    return anchors, lo, hi


def draw_model_3d(ax, zexag=1.0, explode=0.0):
    """Draw the full model. zexag scales z; explode adds spacing between groups."""
    layers, frames, windows, key_z = M.build_stack()

    def Z(z, group):
        return (z * zexag + group * explode)

    skip_gas_3d = {"FrontGas", "BackGas", "DriftCathode_Gas"}
    for L in layers + frames + windows:
        if L.material == "gas" and L.name in skip_gas_3d:
            continue  # flat gas slabs clutter the render; frame shows the gap
        if "Window_Gas" in L.name:
            L.alpha = 0.10  # keep the bulge volume, but very faint
        g = grp(L.name)
        hole = getattr(L, "hole", None)
        k = 64 if L.t > 0.05 or hole is not None else 48
        add_prism(ax, L.poly, Z(L.z0, g), Z(L.z0 + L.t, g),
                  L.color, min(max(L.alpha, 0.35), 0.95), hole=hole, k=k)

    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_xlim(-30, 680); ax.set_ylim(-60, 580)
    return key_z


def fig_3d(outdir):
    fig = plt.figure(figsize=(19, 6.8))
    # +z is downstream, so the front (drift/window) side is at negative z:
    # view it from below (negative elevation).
    views = [("front / drift side", -34, -118, 6),
             ("edge-on — window bulge", 2, -90, 6),
             ("back / readout side", 30, 62, 6)]
    for i, (title, elev, azim, zex) in enumerate(views, 1):
        ax = fig.add_subplot(1, 3, i, projection="3d")
        draw_model_3d(ax, zexag=zex)
        ax.set_zlim(-90, 160)
        ax.set_box_aspect((710, 640, 250))
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(f"{title}  (z ×{zex})", fontsize=11)
        ax.set_zlabel("")
        ax.set_zticks([])
        if "edge-on" in title:
            ax.set_yticks([]); ax.set_ylabel("")
    fig.suptitle("P2 wedge Geant4 model — 3D (z exaggerated for visibility)",
                 fontsize=13)
    fig.tight_layout()
    out = os.path.join(outdir, "p2_3d_overview.png")
    fig.savefig(out, dpi=160); plt.close(fig)
    print("wrote", out)

    zexag, explode = 4, 28
    anchors, _, _ = group_anchors(zexag, explode)
    fig = plt.figure(figsize=(11.5, 7.0))
    ax = fig.add_subplot(projection="3d")
    draw_model_3d(ax, zexag=zexag, explode=explode)
    ax.set_zlim(-160, 260)
    ax.set_box_aspect((710, 640, 420))
    ax.view_init(elev=18, azim=-105)
    ax.set_axis_off()  # presentation view: the model floats, no x/y/z axes
    fig.subplots_adjust(left=-0.13, right=0.70, top=1.10, bottom=-0.38)
    fig.canvas.draw()

    # Leaders land on the right-hand silhouette of each piece and run out to a
    # label column; the text sits at the anchor height, pushed apart where two
    # pieces project too close together to give the type room.
    tips, bottom = [], 1.0
    for lab, groups in EXPLODE_LABELS:
        z = sum(anchors[g][0] for g in groups) / len(groups)
        poly = max((anchors[g][1] for g in groups), key=_poly_area)
        pts = [_fig_xy(fig, ax, x, y, z) for x, y in resample(poly, 96)]
        tips.append(max(pts, key=lambda p: p[0]))
        bottom = min(bottom, min(p[1] for p in pts))
    # Stack the labels by measured height so the whole column is comfortably
    # shorter than the figure, then slide it down if it would run off the top.
    figh = fig.get_size_inches()[1]
    hs = [(lab.count("\n") + 1) * LABEL_FS * 1.35 / 72 / figh
          for lab, _ in EXPLODE_LABELS]
    pad = 0.022
    ys = [t[1] for t in tips]
    ys[0] = max(ys[0], bottom + 0.03 + hs[0] / 2)
    for i in range(1, len(ys)):  # tips run bottom-to-top in stack order
        ys[i] = max(ys[i], ys[i - 1] + (hs[i - 1] + hs[i]) / 2 + pad)
    over = (ys[-1] + hs[-1] / 2) - 0.97
    if over > 0:
        ys = [y - over for y in ys]
    for (lab, _g), (tx, ty), y in zip(EXPLODE_LABELS, tips, ys):
        ax.annotate(lab, xy=(tx, ty), xycoords="figure fraction",
                    xytext=(0.735, y), textcoords="figure fraction",
                    ha="left", va="center", fontsize=LABEL_FS,
                    # relpos pins the leader to the left edge of the text block
                    # instead of letting it leave from under the words.
                    arrowprops=dict(arrowstyle="-", lw=1.1, color="0.25",
                                    shrinkA=6, shrinkB=3, relpos=(0.0, 0.5)))
    out = os.path.join(outdir, "p2_3d_exploded.png")
    fig.savefig(out, dpi=160)  # no tight bbox: it crops the outer labels
    plt.close(fig)
    print("wrote", out)

    fig_3d_exploded_lr(outdir)


def _fig_xy(fig, ax, x, y, z):
    """3D data point -> figure-fraction coordinates (roll included)."""
    xs, ys, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
    return fig.transFigure.inverted().transform(ax.transData.transform((xs, ys)))


def fig_3d_exploded_lr(outdir):
    """The exploded view rolled 90 deg, so the stack reads left to right, with
    the labels off to the side on leader arrows. No axes — it floats."""
    zexag, explode = 4, 62
    anchors, zlo, zhi = group_anchors(zexag, explode)
    pad = 0.06 * (zhi - zlo)

    fig = plt.figure(figsize=(13.5, 10))
    ax = fig.add_subplot(projection="3d")
    draw_model_3d(ax, zexag=zexag, explode=explode)
    ax.set_zlim(zlo - pad, zhi + pad)
    ax.set_box_aspect((710, 640, (zhi - zlo + 2 * pad)))
    # roll=-90 lays the stack down horizontally, front window (beam side) left.
    ax.view_init(elev=18, azim=-105, roll=-90)
    ax.set_axis_off()
    fig.subplots_adjust(left=-0.03, right=1.03, top=1.02, bottom=-0.02)

    # Anchor each label on the mid-plane of its group, at the wedge centroid.
    # The text alternates top / bottom and each row is spread evenly across the
    # frame: anchors run left to right in group order, so the leaders never
    # cross even though the text is not directly over its piece.
    fig.canvas.draw()
    rows = {True: [i for i in range(len(EXPLODE_LABELS)) if i % 2 == 0],
            False: [i for i in range(len(EXPLODE_LABELS)) if i % 2 == 1]}
    slot = {}
    for idx in rows.values():
        for j, i in enumerate(idx):
            slot[i] = 0.18 + (0.64 * j / max(len(idx) - 1, 1))
    for i, (lab, groups) in enumerate(EXPLODE_LABELS):
        z = sum(anchors[g][0] for g in groups) / len(groups)
        poly = max((anchors[g][1] for g in groups), key=_poly_area)
        top = (i % 2 == 0)
        # Land the leader on the silhouette edge of the piece facing its row,
        # so it never has to cross the stack to reach the label.
        pts = [_fig_xy(fig, ax, x, y, z) for x, y in resample(poly, 96)]
        fx, fy = (max if top else min)(pts, key=lambda p: p[1])
        ax.annotate(lab, xy=(fx, fy), xycoords="figure fraction",
                    xytext=(slot[i], 0.975 if top else 0.025),
                    textcoords="figure fraction",
                    ha="center", va="top" if top else "bottom", fontsize=LABEL_FS,
                    arrowprops=dict(arrowstyle="-", lw=1.1, color="0.25",
                                    shrinkA=6, shrinkB=3))
    out = os.path.join(outdir, "p2_3d_exploded_lr.png")
    fig.savefig(out, dpi=160)  # no tight bbox: it crops the outer labels
    plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5 — stack schematic with the open questions (for the collaboration)
# ─────────────────────────────────────────────────────────────────────────────

C_OK, C_MEAS, C_ASSUME, C_GUESS = "#2e7d32", "#1565c0", "#ef6c00", "#c62828"


def fig_questions(outdir):
    fig, ax = plt.subplots(figsize=(15.5, 11.5))

    X0, X1 = 1.30, 4.30                      # stack bars
    XL = X0 - 0.12                           # left annotation anchor
    XQ = 5.55                                # question boxes

    # (key, label, thickness text, height, facecolor, status)
    rows = [
        ("fwin",  "front gas window — mylar", "10 µm",      0.95, "#63b8d8", C_OK),
        ("fgas",  "front gas gap",            "3.88 mm",    1.20, "#dff2fa", C_MEAS),
        ("cath",  "drift cathode — Al-mylar", "120 µm + 1 µm Al",
                                                            0.60, "#a8e0a8", C_OK),
        ("drift", "DRIFT GAS",                "4 mm",       1.30, "#b3ccff", C_MEAS),
        ("mesh",  "micromesh — woven SS 45/18 µm",
                                              "",           0.55, "#9e9e9e", C_OK),
        ("amp",   "amplification gap",        "150 µm",     0.80, "#ffc4c4", C_OK),
        ("fcu",   "F.Cu readout pads",        "18 µm × 0.98", 0.45, "#cc6619", C_MEAS),
        ("fr4",   "FR4 core",                 "200 µm",     0.70, "#5aa85a", C_MEAS),
        ("bcu",   "B.Cu signal lines",        "18 µm × 0.17", 0.45, "#cc6619", C_MEAS),
        ("bgas",  "back gas gap (carbon frame)", "1 mm",   1.00, "#dff2fa", C_OK),
        ("bwin",  "back gas window — mylar",  "10 µm",      0.95, "#63b8d8", C_OK),
    ]

    gap = 0.07
    ytop = 11.9
    y = ytop
    pos = {}
    for key, lab, tlab, h, fc, st in rows:
        y -= h
        pos[key] = (y, h)
        if key in ("fwin", "bwin"):
            # curved band suggesting the overpressure bulge (outward = away
            # from the gas volume): centreline c(x) = flat edge + sag*sin
            xs = np.linspace(X0, X1, 80)
            b = np.sin(np.pi * (xs - X0) / (X1 - X0))
            if key == "fwin":
                c = (y + 0.19*h) + 0.62*h*b
            else:
                c = (y + h - 0.19*h) - 0.62*h*b
            verts = ([(x, ci + 0.14*h) for x, ci in zip(xs, c)] +
                     [(x, ci - 0.14*h) for x, ci in zip(xs[::-1], c[::-1])])
            ax.add_patch(MplPolygon(np.array(verts), closed=True, facecolor=fc,
                                    edgecolor=st, lw=2.2))
            ax.text(XL - 0.55, y + h/2, f"{lab}\n{tlab}, bulged",
                    ha="right", va="center", fontsize=10.5, color="k")
        else:
            ax.add_patch(Rectangle((X0, y), X1 - X0, h, facecolor=fc,
                                   edgecolor=st, lw=2.2))
            ax.text((X0+X1)/2, y + h/2,
                    f"{lab}   —   {tlab}" if tlab else lab,
                    ha="center", va="center",
                    fontsize=11.5 if h >= 0.6 else 10.0,
                    fontweight="bold" if key == "drift" else "normal")
        y -= gap

    # gas frame side rails (front half: window plane -> pad plane; back half)
    yf0 = pos["fcu"][0] + pos["fcu"][1] + gap          # pad plane
    yf1 = pos["fwin"][0]                                # front window plane
    yb0 = pos["bwin"][0] + pos["bwin"][1]
    yb1 = pos["bcu"][0] - gap
    for xr in (X0 - 0.42, X1 + 0.06):
        ax.add_patch(Rectangle((xr, yf0), 0.36, yf1 - yf0, facecolor="#e6e6d9",
                               edgecolor=C_ASSUME, lw=2.0))
        ax.add_patch(Rectangle((xr, yb0), 0.36, yb1 - yb0, facecolor="#5a5a5a",
                               edgecolor=C_OK, lw=2.0))
    ax.text(X0 - 0.24, (yf0+yf1)/2, "gas frame (STEP)", rotation=90,
            ha="center", va="center", fontsize=9)
    ax.text(X0 - 0.24, (yb0+yb1)/2, "carbon frame", rotation=90,
            ha="center", va="center", fontsize=9, color="w")

    # beam arrow
    ax.annotate("", xy=((X0+X1)/2, ytop + 0.55), xytext=((X0+X1)/2, ytop + 1.55),
                arrowprops=dict(arrowstyle="-|>", lw=2.5, color="k"))
    ax.text((X0+X1)/2 + 0.12, ytop + 1.05, "e⁻ (from target side)",
            fontsize=11, va="center")

    # ── left side: what is already fixed ──────────────────────────────────
    def left_note(key, txt, color):
        yy, hh = pos[key]
        ax.text(XL - 0.55, yy + hh/2, txt, ha="right", va="center",
                fontsize=10.5, color=color,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=color, lw=1.4))
    left_note("drift", "✓ frame V2 ledge height (STEP)\n(campaign will scan 1–4 mm)", C_MEAS)
    left_note("amp",   "✓ 150 µm confirmed", C_OK)
    left_note("mesh",  "✓ woven SS: 18 µm wire, 45 µm opening;\nmodelled as 36 µm effective-density slab\n(~51 % optical transparency not modelled)", C_OK)
    left_note("fr4",   "✓ PCB from Stack_Up_P2.txt; Cu layers\ndensity-scaled by gerber-measured coverage\n(B.Cu radial 0.03 → 0.26, mean 0.17)", C_MEAS)
    left_note("fgas",  "✓ frame STEP V2: 8 mm body, foil on the\nledge 4 mm above the board; one foil\n(two-foil model retired 2026-10-06)", C_MEAS)
    left_note("bgas",  "✓ carbon frame glued to PCB back,\nmylar on its rear face; 1 mm confirmed\n(a few of this production are 3 mm)", C_OK)

    # ── right side: the questions ─────────────────────────────────────────
    def qbox(n, y_center, txt, color, targets):
        bb = ax.text(XQ, y_center, f"Q{n}.  {txt}", ha="left", va="center",
                     fontsize=10.8, color="k",
                     bbox=dict(boxstyle="round,pad=0.45", fc="#fffdf5",
                               ec=color, lw=1.8))
        for key in targets:
            yy, hh = pos[key]
            ax.annotate("", xy=(X1 + 0.55, yy + hh/2), xytext=(XQ - 0.06, y_center),
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.3,
                                        shrinkA=2, shrinkB=1,
                                        connectionstyle="arc3,rad=-0.08"))

    qbox(1, 10.75,
         "Carbon back frame: inner outline (opening) and wall width?\n"
         "Depth 1 mm confirmed; sim assumes the front frame's\n"
         "footprint and window-side opening (r 107→603, edges +3).",
         C_ASSUME, ["bgas"])
    qbox(2, 8.55,
         "Window sag 10 mm on both sides (2026-10-06):\n"
         "modelled as a 6-step terraced dome of mylar,\n"
         "chamber gas filling it.",
         C_OK, ["fwin", "bwin"])
    qbox(3, 6.55,
         "What does the frame sit on — bare PCB or the bulk border?\n"
         "STEP ledge is 4.0 mm above the frame's own face: drift gas is\n"
         "3.81 mm (bare PCB) to ~3.96 mm (bulk); sim uses 4.0 mm.",
         C_ASSUME, ["drift"])

    # legend for the status colors
    for i, (c, lab) in enumerate([
            (C_OK,     "confirmed"),
            (C_MEAS,   "measured from fab files"),
            (C_ASSUME, "assumed — please check"),
            (C_GUESS,  "unknown — best guess")]):
        ax.add_patch(Rectangle((0.10, 0.62 - 0.55*i), 0.30, 0.30,
                               facecolor="white", edgecolor=c, lw=2.2))
        ax.text(0.50, 0.77 - 0.55*i, lab, fontsize=10.5, va="center")

    ax.set_xlim(-1.6, 11.6)
    ax.set_ylim(-1.35, ytop + 2.0)
    ax.axis("off")
    ax.set_title("P2 wedge gas envelope — stack schematic (not to scale) and open questions",
                 fontsize=14, pad=12)
    fig.tight_layout()
    out = os.path.join(outdir, "p2_stack_questions.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print("wrote", out)


# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.join(REPO, "docs", "figures"))
    ap.add_argument("--only", choices=["board", "xsec", "3d", "questions"],
                    default=None)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    if args.only in (None, "board"):
        fig_board(args.out_dir)
    if args.only in (None, "xsec"):
        fig_xsec(args.out_dir)
    if args.only in (None, "3d"):
        fig_3d(args.out_dir)
    if args.only in (None, "questions"):
        fig_questions(args.out_dir)


if __name__ == "__main__":
    main()
