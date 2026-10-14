#!/usr/bin/env python3
"""Figures for the P2 mesh field map and transparency (task P0.13).

Reads the outputs of scripts/garfield/mesh_transparency.cc (copied from
/eos/user/a/akallits/p2_sim/transparency/) and writes docs/results/mesh/:

  fig_mesh_xz.png      side views: |E| colour map, field lines, wires
  fig_mesh_xy.png      top views of |E| above, through and below the wires
  fig_mesh_profile.png Ez along z through an opening and over a crossing
  fig_transp_lost.png  where the lost electrons land on the wires (top view)
  fig_transp_bars.png  transparency per gas, 4 mm vs 2 mm cell
  fig_transp_sweep.png transparency vs drift field (sweep/), with the
                       Nikolopoulos et al. measurement for comparison
  fig_transp_sweep_gap150.png  the same, corrected cell (sweep_gap150/)
  fig_transp_ek.png    wire losses vs electron characteristic energy at the
                       SPS point (why the gases differ)
  fig_transp_nikolopoulos.png  their geometry and gas rebuilt (sweep_nik/),
                       default vs null-collision stepping, against their data
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
         "ArCO2Iso9352": "Ar/CO2/iso 93/5/2", "ArCF4Iso": "Ar/CF4/iso 88/10/2",
         "ArCO2_8515": "Ar/CO2 85/15"}
E_AMP = 31.5e3  # V/cm, our cell (mesh -450 V over 150 um, from the probe)

# Nikolopoulos, Bhattacharya, Chernyatin, Veenhof, JINST 6 (2011) P06011,
# fig. 5 (left), read off the figure (+-0.02): 55Fe amplitude normalised to
# E_drift < 600 V/cm. Ar/CO2 85/15, 18 um calendered woven mesh, 63.5 um
# pitch, 128 um gap, E_amp = 41.4 kV/cm.
NIK_EAMP = 41.4e3
NIK = [(124, 1.01), (275, 0.99), (521, 0.99), (760, 0.995), (1000, 0.97),
       (1242, 0.89), (1500, 0.84), (1735, 0.78), (1946, 0.74), (2187, 0.71),
       (2460, 0.68), (2877, 0.62), (3326, 0.59), (3846, 0.55), (4385, 0.52),
       (4853, 0.47), (5450, 0.42), (6033, 0.40)]
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
            att = ((e[:, 6] < 1) & (e[:, 2] == -7)).mean()
            rows.append({"cell": cell, "gas": gas, "n": n,
                         "T_micro": round(tm, 4),
                         "err": round(math.sqrt(tm * (1 - tm) / n), 4),
                         "attached": round(att, 4),
                         "lost_on_wires": round(1 - tm - att, 4),
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
        a = [100 * r[g]["attached"] if g in r else 0 for g in GASES]
        ax.bar(x + (k - 0.5) * 0.38, v, 0.36, yerr=e, color=col, label=CELLS[cell],
               capsize=2)
        ax.bar(x + (k - 0.5) * 0.38, a, 0.36, bottom=v, color="none",
               edgecolor="#D62728", hatch="///", lw=0,
               label="lost to attachment (gas, not mesh)" if k == 0 else None)
        for xi, vi in zip(x + (k - 0.5) * 0.38, v):
            if not np.isnan(vi):
                ax.text(xi, vi + 0.6, f"{vi:.1f}", ha="center", fontsize=8)
    ax.set_xticks(x, [LABEL[g] for g in GASES], rotation=15)
    ax.set_ylim(70, 102)
    ax.set_ylabel("electrons reaching the gap [%]")
    ax.legend(loc="lower right")
    ax.set_title("Electrons reaching the gap (Garfield microscopic, mesh −450 V, drift 632 V/cm); "
                 "field lines alone give 100 %", fontsize=10)
    fig.savefig(out / "fig_transp_bars.png", dpi=160)
    plt.close(fig)
    for q in rows:
        print(q)


def fig_sweep(d, out, sub="sweep", e_amp=None, tag="", note=""):
    e_amp = e_amp or E_AMP
    sd = d / sub
    if not sd.exists():
        return
    res = {}
    for p in sorted(sd.glob("E*_*_electrons.txt")):
        e_str, gas = p.name.split("_", 1)
        gas = gas.replace("_electrons.txt", "")
        e = read_transp(p)
        n = len(e)
        res.setdefault(gas, []).append((float(e_str[1:]), e[:, 6].mean(),
                                        (e[:, 2] == -7).mean(),
                                        e[:, 11].mean() if e[:, 7].any() else np.nan, n))
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
    cols = {"ArIso": "#D62728", "ArIso9010": "#FF9F80", "NeIso": "#2CA02C",
            "ArCO2Iso9352": "#9467BD", "ArCO2_8515": "#1F77B4",
            "ArCF4Iso": "#8C564B", "NeIso8515": "#17BECF"}
    for gas, pts in res.items():
        pts.sort()
        E = np.array([q[0] for q in pts])
        T = np.array([q[1] for q in pts])
        A = np.array([q[2] for q in pts])
        R = np.array([q[3] for q in pts])
        n = np.array([q[4] for q in pts])
        err = np.sqrt(T * (1 - T) / n)
        c = cols.get(gas, "k")
        axes[0].errorbar(E / e_amp, 100 * T, 100 * err, marker="o", ms=4,
                         color=c, label=LABEL[gas])
        if not np.all(np.isnan(R)):
            axes[0].plot(E / e_amp, 100 * R, ":", color=c, lw=1)
        # literature convention: normalise to the low-field plateau
        ref = T[E <= 600].mean() if np.any(E <= 600) else T.max()
        axes[1].plot(E / e_amp, T / ref, "-o", ms=4, color=c, label=LABEL[gas])
        for q, t, a, r in zip(E, T, A, R):
            rows.append({"E_drift": q, "gas": gas, "ratio": round(q / e_amp, 4),
                         "T_micro": round(t, 4), "attached": round(a, 4),
                         "T_fieldlines": None if np.isnan(r) else round(r, 4)})
    ne = np.array(NIK)
    for ax in axes:
        ax.set_xscale("log")
        ax.axvline(625 / e_amp, color="#888888", ls="--", lw=1)
        ax.text(625 / e_amp, 0.03, " P2 at SPS", color="#555555", fontsize=9,
                transform=ax.get_xaxis_transform(), va="bottom")
        ax.set_xlabel("E_drift / E_amp")
    axes[1].plot(ne[:, 0] / NIK_EAMP, ne[:, 1], "s", color="k", ms=5,
                 label="measured, Nikolopoulos 2011\n(Ar/CO2 85/15, 63.5 µm pitch)")
    axes[0].set_ylabel("electrons reaching the gap [%] (absolute)")
    axes[0].set_title("Absolute transparency (dotted: field lines only)", fontsize=11)
    axes[1].set_ylabel("relative to E_drift ≤ 600 V/cm")
    axes[1].set_title("As the literature quotes it: normalised to the low-field plateau",
                      fontsize=11)
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    fig.suptitle(f"Mesh transparency vs drift field (2 mm cell, E_amp = {e_amp/1e3:.1f} kV/cm, "
                 f"18 µm wire / 68 µm pitch{note})", fontsize=11)
    fig.savefig(out / f"fig_transp_sweep{tag}.png", dpi=160)
    plt.close(fig)
    with open(out / f"transparency_sweep{tag}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["gas"], r["E_drift"])))


def fig_nik(d, out):
    sd = d / "sweep_nik"
    if not sd.exists():
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
    ne = np.array(NIK)
    rows = []
    for mode, c, lab in [("x", "#D62728", "default free flights"),
                         ("xnull", "#2F6DB5", "null-collision steps")]:
        pts = []
        for p in sorted((sd / mode).glob("E*_ArCO2_8515_electrons.txt")):
            e = read_transp(p)
            # nominal -> actual drift field in the cell (probe: 614 V/cm at
            # the nominal 600, the gap is counted from the wire centres)
            pts.append((1.023 * float(p.name.split("_")[0][1:]), e[:, 6].mean(), len(e),
                        (e[:, 2] == -7).mean()))
        if not pts:
            continue
        pts.sort()
        E = np.array([q[0] for q in pts])
        T = np.array([q[1] for q in pts])
        n = np.array([q[2] for q in pts])
        err = np.sqrt(T * (1 - T) / n)
        axes[0].errorbar(E, 100 * T, 100 * err, marker="o", ms=4, color=c, label=lab)
        ref = T[E <= 600].mean()
        axes[1].errorbar(E, T / ref, err / ref, marker="o", ms=4, color=c,
                         label=f"{lab} (plateau {100 * ref:.1f} %)")
        dev = T / ref - np.interp(E, ne[:, 0], ne[:, 1])
        for q, t, dv, a in zip(E, T, dev, [q[3] for q in pts]):
            rows.append({"mode": mode, "E_drift": q, "T_micro": round(t, 4),
                         "T_rel": round(t / ref, 4), "attached": round(a, 4),
                         "rel_minus_data": round(dv, 4) if q >= ne[0, 0] else None})
        print(f"{mode}: rms(sim/plateau - data) for E>=600 V/cm = "
              f"{np.sqrt(np.mean(dev[E >= 600] ** 2)):.3f}")
    axes[1].plot(ne[:, 0], ne[:, 1], "s", color="k", ms=5,
                 label="measured (55Fe), Nikolopoulos 2011")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("E_drift [V/cm]")
        ax.legend(fontsize=9)
    axes[0].set_ylabel("electrons reaching the gap [%] (absolute)")
    axes[0].set_title("Our simulation, absolute", fontsize=11)
    axes[1].set_ylabel("relative to E_drift ≤ 600 V/cm")
    axes[1].set_title("Against the measurement (same normalisation)", fontsize=11)
    fig.suptitle("Validation: Nikolopoulos et al. geometry rebuilt (18 µm wire, 63.5 µm, "
                 "calendered 30 µm, 128 µm gap, E_amp 41.4 kV/cm), Ar/CO2 85/15",
                 fontsize=11)
    fig.savefig(out / "fig_transp_nikolopoulos.png", dpi=160)
    plt.close(fig)
    if rows:
        with open(out / "transparency_nikolopoulos.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


def fig_review(d, out):
    """SPS point per gas: first pass -> null-collision steps -> mesh bottom
    on the 150 um pillars (the 2026-10-14 review)."""
    steps = [("first pass: default steps, mesh centre at 150 µm",
              lambda g: d / f"transp_h4mm_{g}_electrons.txt", "#BBBBBB"),
             ("+ null-collision steps", lambda g: d / "check" / f"xnull_h4mm_{g}_electrons.txt",
              "#9DBFE3"),
             ("+ mesh bottom on the 150 µm pillars (final)",
              lambda g: d / "sweep_gap150" / f"E625_{g}_electrons.txt", "#2F6DB5")]
    if not all(f(GASES[0]).exists() for _, f, _ in steps):
        return
    fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
    x = np.arange(len(GASES))
    rows = []
    for k, (lab, f, c) in enumerate(steps):
        T, A, E = [], [], []
        for g in GASES:
            e = read_transp(f(g))
            t = e[:, 6].mean()
            T.append(100 * t)
            A.append(100 * ((e[:, 6] < 1) & (e[:, 2] == -7)).mean())
            E.append(100 * math.sqrt(t * (1 - t) / len(e)))
            rows.append({"step": k, "gas": g, "T": round(t, 4), "attached": round(A[-1] / 100, 4),
                         "n": len(e)})
        xs = x + (k - 1) * 0.27
        ax.bar(xs, T, 0.25, yerr=E, color=c, label=lab, capsize=2)
        ax.bar(xs, A, 0.25, bottom=T, color="none", edgecolor="#D62728", hatch="///",
               lw=0, label="lost to attachment (gas)" if k == 2 else None)
        if k == 2:
            for xi, v in zip(xs, T):
                ax.text(xi, v + 0.4, f"{v:.1f}", ha="center", fontsize=8)
    ax.set_xticks(x, [LABEL[g] for g in GASES], rotation=12)
    ax.set_ylim(85, 101.5)
    ax.set_ylabel("electrons reaching the gap [%]")
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title("SPS point (mesh −450 V, drift 625 V/cm): how each correction moved the number",
                 fontsize=11)
    fig.savefig(out / "fig_transp_review.png", dpi=160)
    plt.close(fig)
    with open(out / "transparency_review.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def fig_ek(d, out, e_drift=638.0):
    """Wire losses at the SPS point vs the characteristic energy
    eps_k = e D_T / mu = D_T^2 E / 2 (D_T in sqrt(cm), E in V/cm -> eV),
    from the Magboltz tables at the drift field in the corrected cell."""
    sd = d / "sweep_gap150"
    js = pathlib.Path(__file__).resolve().parents[2] / "stage_b" / "gas_tables" / \
        "magboltz_dry_20261008.json"
    if not sd.exists() or not js.exists():
        return
    import json
    tab = json.load(open(js))
    fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    pts = []
    for g in GASES:
        p = sd / f"E625_{g}_electrons.txt"
        if not p.exists() or g not in tab:
            continue
        e = read_transp(p)
        wire = ((e[:, 6] < 1) & (e[:, 2] == -5)).mean()
        err = math.sqrt(max(wire, 1 / len(e)) * (1 - wire) / len(e))
        t = tab[g]
        E = np.array(t["E_V_per_cm"])
        dt = np.exp(np.interp(np.log(e_drift), np.log(E), np.log(t["DT_um_sqrtcm"])))
        ek = (dt * 1e-4) ** 2 * e_drift / 2
        pts.append((ek, 100 * wire, 100 * err, g, dt))
    cols = {"ArIso": "#D62728", "ArIso9010": "#FF9F80", "NeIso": "#2CA02C",
            "NeIso8515": "#17BECF", "ArCO2Iso9352": "#9467BD", "ArCF4Iso": "#8C564B"}
    for ek, w, er, g, dt in pts:
        ax.errorbar(ek, w, er, marker="o", ms=9, color=cols[g], capsize=3)
        off = {"NeIso": (10, 8), "ArCF4Iso": (-60, 22), "NeIso8515": (-30, 30)}.get(g, (10, -4))
        ax.annotate(f"{LABEL[g]}\nD_T = {dt:.0f} µm/√cm", (ek, w), textcoords="offset points",
                    xytext=off, fontsize=9)
    ax.axvline(0.0253, color="#888888", ls=":", lw=1)
    ax.text(0.03, 0.95, "thermal\n(0.025 eV)", transform=ax.get_xaxis_transform(),
            fontsize=8, color="#555555", va="top")
    ax.set_xscale("log")
    ax.set_xlim(0.02, 1.6)
    ax.set_ylim(-0.5, 8.5)
    ax.set_xlabel("electron characteristic energy ε_k = e·D_T/μ at 638 V/cm [eV]  (Magboltz)")
    ax.set_ylabel("electrons stopped on the wires [%]")
    ax.set_title("Hotter electrons diffuse further off the field lines and hit the mesh\n"
                 "(SPS point, mesh bottom on the 150 µm pillars; attachment not included)",
                 fontsize=11)
    fig.savefig(out / "fig_transp_ek.png", dpi=160)
    plt.close(fig)
    for ek, w, er, g, dt in pts:
        print(f"{g:13s} eps_k {ek:.3f} eV  wire {w:.2f} +- {er:.2f} %")


def fig_gaps(d, out):
    """Drift-gap scan at a fixed 250 V above the mesh (gaps/)."""
    sd = d / "gaps"
    if not sd.exists():
        return
    fig, ax = plt.subplots(figsize=(8.5, 5), constrained_layout=True)
    cols = {"ArIso": "#D62728", "ArIso9010": "#FF9F80", "NeIso": "#2CA02C",
            "NeIso8515": "#17BECF", "ArCO2Iso9352": "#9467BD", "ArCF4Iso": "#8C564B"}
    rows = []
    for g in GASES:
        pts = []
        for gap in (4, 5, 6, 7, 8):
            p = sd / f"h{gap}mm_{g}_electrons.txt"
            if not p.exists():
                continue
            e = read_transp(p)
            t = e[:, 6].mean()
            wire = ((e[:, 6] < 1) & (e[:, 2] == -5)).mean()
            pts.append((gap, t, wire, math.sqrt(t * (1 - t) / len(e))))
            rows.append({"gap_mm": gap, "E_drift": round(250 / (gap * 0.1)), "gas": g,
                         "T": round(t, 4), "lost_on_wires": round(wire, 4),
                         "attached": round(1 - t - wire, 4)})
        if pts:
            a = np.array(pts)
            ax.errorbar(a[:, 0], 100 * a[:, 1], 100 * a[:, 3], marker="o", color=cols[g],
                        label=LABEL[g], capsize=2)
    ax.set_xticks([4, 5, 6, 7, 8],
                  [f"{g} mm\n{250 / (g * 0.1):.0f} V/cm" for g in (4, 5, 6, 7, 8)])
    ax.set_ylabel("electrons reaching the gap [%]")
    ax.set_ylim(90, 100.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("Drift gap at a fixed 250 V above the mesh (mesh −450 V): longer gap = weaker field",
                 fontsize=10)
    fig.savefig(out / "fig_transp_gaps.png", dpi=160)
    plt.close(fig)
    with open(out / "transparency_gaps.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


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
    fig_sweep(d, out)
    # corrected cell: mesh bottom on the 150 um pillars, null-collision steps
    fig_sweep(d, out, "sweep_gap150", 28.5e3, "_gap150",
              "; mesh bottom at 150 µm; null-collision steps")
    fig_nik(d, out)
    fig_review(d, out)
    fig_ek(d, out)
    fig_gaps(d, out)


if __name__ == "__main__":
    main()
