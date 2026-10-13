#!/usr/bin/env python3
"""Figures for the P2 mesh field map and transparency (task P0.13).

Reads the outputs of scripts/garfield/mesh_transparency.cc (copied from
/eos/user/a/akallits/p2_sim/transparency/) and writes docs/results/mesh/:

  fig_mesh_xz.png      side views: |E| colour map, field lines, wires
  fig_mesh_xy.png      top views of |E| above, through and below the wires
  fig_mesh_profile.png Ez along z through an opening and over a crossing
  fig_transp_lost.png  where the lost electrons land on the wires (top view)
  fig_transp_bars.png  transparency per gas, 4 mm vs 2 mm cell
  transparency.csv     the numbers

    python3 scripts/analysis/mesh_field_maps.py <dir with the .txt files>
"""
import csv
import math
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, ListedColormap

GASES = ["ArIso", "ArIso9010", "NeIso", "NeIso8515", "ArCO2Iso9352", "ArCF4Iso"]
LABEL = {"ArIso": "Ar/iso 95/5", "ArIso9010": "Ar/iso 90/10",
         "NeIso": "Ne/iso 95/5", "NeIso8515": "Ne/iso 85/15",
         "ArCO2Iso9352": "Ar/CO2/iso 93/5/2", "ArCF4Iso": "Ar/CF4/iso 88/10/2"}
CELLS = {"h4mm": "4 mm drift", "h2mm": "2 mm drift"}
SOLID = ListedColormap(["#555555"])

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11})


def load_grid(path):
    d = np.loadtxt(path)
    a, b = np.unique(d[:, 0]), np.unique(d[:, 1])
    shape = (len(b), len(a))
    g = {k: d[:, i].reshape(shape) for i, k in
         enumerate(["a", "b", "V", "Ex", "Ey", "Ez", "gas"])}
    g["a1"], g["b1"] = a, b
    g["E"] = np.sqrt(g["Ex"] ** 2 + g["Ey"] ** 2 + g["Ez"] ** 2)
    return g


def draw_side(ax, g, lines, y0, title, zlim=(60, 330)):
    E = np.where(g["gas"] > 0, g["E"], np.nan)
    im = ax.pcolormesh(g["a1"], g["b1"], E / 1e3, cmap="viridis",
                       norm=LogNorm(0.3, 100), shading="auto")
    solid = np.where(g["gas"] > 0, np.nan, 1.0)
    ax.pcolormesh(g["a1"], g["b1"], solid, cmap=SOLID, shading="auto")
    cs = ax.contour(g["a1"], g["b1"], g["V"], levels=np.arange(-500, 0, 10),
                    colors="white", linewidths=0.4, alpha=0.6)
    if lines is not None:
        sel = lines[np.isclose(lines[:, 1], y0)]
        for lid in np.unique(sel[:, 0]):
            p = sel[sel[:, 0] == lid]
            ax.plot(p[:, 2], p[:, 4], color="#FF5A36", lw=0.8)
    ax.set_xlim(g["a1"][0], g["a1"][-1])
    ax.set_ylim(*zlim)
    ax.set_aspect("equal")
    ax.set_xlabel("x [µm]")
    ax.set_ylabel("z [µm]  (pads at 0)")
    ax.set_title(title, fontsize=11)
    return im


def fig_side(d, out, cell):
    lines_p = d / f"maps_{cell}_ArIso_lines.txt"
    lines = np.loadtxt(lines_p) if lines_p.exists() else None
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.2), constrained_layout=True)
    for ax, y0, t in [(axes[0], 0, "y = 0: through the openings\n(wires along y seen end-on at x = ±34 µm)"),
                      (axes[1], 34, "y = 34 µm: along a wire\n(it weaves over and under the crossing wires)")]:
        g = load_grid(d / f"maps_{cell}_ArIso_xz_y{y0}.txt")
        im = draw_side(ax, g, lines, y0, t)
    cb = fig.colorbar(im, ax=axes, shrink=0.8)
    cb.set_label("|E| [kV/cm]")
    fig.suptitle(f"P2 woven mesh, COMSOL unit cell ({CELLS[cell]}): field strength (colour), "
                 "equipotentials every 10 V (white), electron field lines (red), wires (grey)",
                 fontsize=11)
    fig.savefig(out / f"fig_mesh_xz_{cell}.png", dpi=160)
    plt.close(fig)


def fig_top(d, out, cell):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), constrained_layout=True)
    for ax, z, t in [(axes[0], 185, "z = 185 µm: 20 µm above the mesh"),
                     (axes[1], 150, "z = 150 µm: through the wires"),
                     (axes[2], 110, "z = 110 µm: 20 µm below the mesh")]:
        g = load_grid(d / f"maps_{cell}_ArIso_xy_z{z}.txt")
        E = np.where(g["gas"] > 0, g["E"], np.nan)
        im = ax.pcolormesh(g["a1"], g["b1"], E / 1e3, cmap="viridis",
                           norm=LogNorm(0.3, 100), shading="auto")
        ax.pcolormesh(g["a1"], g["b1"], np.where(g["gas"] > 0, np.nan, 1.0),
                      cmap=SOLID, shading="auto")
        ax.set_aspect("equal")
        ax.set_title(t, fontsize=11)
        ax.set_xlabel("x [µm]")
        ax.set_ylabel("y [µm]")
    cb = fig.colorbar(im, ax=axes, shrink=0.85)
    cb.set_label("|E| [kV/cm]")
    fig.suptitle(f"Top views of one weave period, 136 × 136 µm ({CELLS[cell]}); grey = wire",
                 fontsize=11)
    fig.savefig(out / f"fig_mesh_xy_{cell}.png", dpi=160)
    plt.close(fig)


def fig_profile(d, out):
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    for cell, ls in [("h4mm", "-"), ("h2mm", "--")]:
        p = d / f"maps_{cell}_ArIso_xz_y0.txt"
        if not p.exists():
            continue
        g = load_grid(p)
        for x0, c, lab in [(0, "#2F6DB5", "centre of an opening"),
                           (17, "#D9822B", "halfway between opening centre and wire")]:
            i = np.argmin(np.abs(g["a1"] - x0))
            ez = np.where(g["gas"][:, i] > 0, g["Ez"][:, i], np.nan)
            ax.plot(g["b1"], ez / 1e3, ls, color=c, lw=1.5,
                    label=f"{lab}, {CELLS[cell]}")
    ax.set_yscale("log")
    ax.set_xlabel("z [µm]  (pads at 0, mesh at 135–165)")
    ax.set_ylabel("Ez [kV/cm]")
    ax.axvspan(135.5, 165, color="#DDDDDD")
    ax.legend(fontsize=9)
    ax.set_title("Field along z at y = 0: 0.63 kV/cm drift, 31 kV/cm amplification")
    fig.savefig(out / "fig_mesh_profile.png", dpi=160)
    plt.close(fig)


def read_transp(p):
    d = np.loadtxt(p)
    return d


def fold(v):
    return (v + 68.0) % 136.0 - 68.0


def fig_lost(d, out):
    cases = [("h4mm", "ArIso"), ("h4mm", "ArIso9010"), ("h4mm", "NeIso")]
    cases = [c for c in cases if (d / f"transp_{c[0]}_{c[1]}_electrons.txt").exists()]
    gp = d / "maps_h4mm_ArIso_xy_z150.txt"
    if not cases or not gp.exists():
        return
    g = load_grid(gp)
    fig, axes = plt.subplots(1, len(cases), figsize=(5 * len(cases), 5.2),
                             constrained_layout=True, squeeze=False)
    for ax, (cell, gas) in zip(axes[0], cases):
        e = read_transp(d / f"transp_{cell}_{gas}_electrons.txt")
        lost = e[:, 6] < 1
        ax.pcolormesh(g["a1"], g["b1"], np.where(g["gas"] > 0, np.nan, 1.0),
                      cmap=ListedColormap(["#CCCCCC"]), shading="auto")
        sc = ax.scatter(fold(e[lost, 3]), fold(e[lost, 4]), s=5,
                        c=e[lost, 5], cmap="coolwarm", vmin=135, vmax=165)
        ax.set_aspect("equal")
        ax.set_xlim(-68, 68)
        ax.set_ylim(-68, 68)
        ax.set_xlabel("x [µm]")
        ax.set_ylabel("y [µm]")
        ax.set_title(f"{LABEL[gas]}: {lost.sum()} of {len(e)} lost "
                     f"(T = {1 - lost.mean():.1%})", fontsize=11)
    cb = fig.colorbar(sc, ax=axes, shrink=0.85)
    cb.set_label("height where it hits the wire [µm]")
    fig.suptitle("Where the lost electrons end: on the upper half of the wires (red), "
                 "spread along them (grey = wires cut at z = 150 µm, 4 mm cell)", fontsize=11)
    fig.savefig(out / "fig_transp_lost.png", dpi=160)
    plt.close(fig)


def fig_bars(d, out):
    rows = []
    for cell in CELLS:
        for gas in GASES:
            p = d / f"transp_{cell}_{gas}_electrons.txt"
            if not p.exists():
                continue
            e = read_transp(p)
            n = len(e)
            tm, tr = e[:, 6].mean(), e[:, 11].mean()
            rows.append({"cell": cell, "gas": gas, "n": n,
                         "T_micro": round(tm, 4),
                         "err": round(math.sqrt(tm * (1 - tm) / n), 4),
                         "T_fieldlines": round(tr, 4)})
    if not rows:
        return
    with open(out / "transparency.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    x = np.arange(len(GASES))
    for k, (cell, col) in enumerate([("h4mm", "#2F6DB5"), ("h2mm", "#9DBFE3")]):
        r = {q["gas"]: q for q in rows if q["cell"] == cell}
        v = [100 * r[g]["T_micro"] if g in r else np.nan for g in GASES]
        e = [100 * r[g]["err"] if g in r else 0 for g in GASES]
        ax.bar(x + (k - 0.5) * 0.38, v, 0.36, yerr=e, color=col, label=CELLS[cell],
               capsize=2)
        for xi, vi in zip(x + (k - 0.5) * 0.38, v):
            if not np.isnan(vi):
                ax.text(xi, vi + 0.6, f"{vi:.1f}", ha="center", fontsize=8)
    ax.set_xticks(x, [LABEL[g] for g in GASES], rotation=15)
    ax.set_ylim(70, 102)
    ax.set_ylabel("electrons reaching the gap [%]")
    ax.legend(loc="lower right")
    ax.set_title("Mesh transparency with diffusion (Garfield microscopic); "
                 "field lines alone give 100 %")
    fig.savefig(out / "fig_transp_bars.png", dpi=160)
    plt.close(fig)
    for q in rows:
        print(q)


def main():
    d = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    out = pathlib.Path(__file__).resolve().parents[2] / "docs" / "results" / "mesh"
    out.mkdir(parents=True, exist_ok=True)
    for cell in CELLS:
        if (d / f"maps_{cell}_ArIso_xz_y0.txt").exists():
            fig_side(d, out, cell)
            fig_top(d, out, cell)
    fig_profile(d, out)
    fig_lost(d, out)
    fig_bars(d, out)


if __name__ == "__main__":
    main()
