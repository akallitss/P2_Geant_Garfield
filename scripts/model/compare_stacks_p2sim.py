#!/usr/bin/env python3
"""
Layer-by-layer comparison of the P2 Micromegas stack in Matthieu Gailliard's
background simulation (P2Sim, detector-test mode, geometry MMP2_v4 from the
p2geometry submodule, branch detTest_filter) with the stack this repo builds
from the design files (scripts/model/p2_model.py, mirror of
DetectorConstruction::ConstructP2).

His stack is READ from the cfg files (z ranges + materials in geometry.cfg),
so the comparison follows any regeneration of MMP2_v4. Ours comes from
p2_model. Areal densities use each simulation's own material densities;
"ours" folds the copper coverage in (F.Cu pads ~0.98 over the pad field,
B.Cu traces ~0.17 band average) because that is the copper a photon sees.

Usage:
    python scripts/model/compare_stacks_p2sim.py [path/to/MMP2_v4]
writes docs/figures/p2_stack_vs_p2sim.png and prints a markdown table.
"""
from __future__ import annotations

import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import p2_model as M  # noqa: E402

DEFAULT = os.path.join(os.path.dirname(REPO), "P2Sim", "detector_geometry", "MMP2_v4")
OUT = os.path.join(REPO, "docs", "figures", "p2_stack_vs_p2sim.png")

# g/cm3, as each simulation defines them
RHO_THEIRS = dict(mylar=1.40, aluminum=2.699, copper=8.96, mmFR4=1.86,
                  mmGas=0.00167008, mmFullMesh=7.93)
RHO_OURS = dict(mylar=1.40, Al=2.699, Cu=8.96, FR4=1.85, gas=0.00169846,
                steel=8.00, Dynamask=1.30, plastic=1.20, carbon=1.60)
FCU_COV, BCU_COV = 0.978, 0.174


def read_cfg(path):
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            out[p[0]] = p[1]
    return out


def read_theirs(d):
    """[(name, z0, z1, material)] for the non-mesh volumes + mesh summary."""
    vols, mesh = [], []
    lines = open(os.path.join(d, "geometry.cfg")).read().splitlines()
    for ln in lines:
        p = ln.split()
        if not p or not p[0].endswith(".cfg") or p[0] in ("air_chamber.cfg", "GasMM.cfg"):
            continue
        c = read_cfg(os.path.join(d, p[0]))
        if p[0].startswith("mesh_tube"):
            mesh.append((p[0], c))
            continue
        vols.append((p[0][:-4], float(c["zMin"]), float(c["zMax"]), p[1]))
    gas = read_cfg(os.path.join(d, "GasMM.cfg"))
    # lateral size: the mesh tubes span xMin..xMax = -4..4 mm (the box cfg's
    # rMax/dY convention -- half or full width -- is not documented)
    vols.sort(key=lambda v: v[1])
    # fill the mother-volume gas between daughters
    full, z = [], 0.0
    for v in vols:
        if v[1] > z + 1e-9:
            full.append(("gas (mother)", z, v[1], "mmGas"))
        full.append(v)
        z = v[2]
    if float(gas["zMax"]) > z + 1e-9:
        full.append(("gas (mother)", z, float(gas["zMax"]), "mmGas"))
    xs = sorted(float(c["xMin"]) for n, c in mesh if "_x_" in n)
    r = float(mesh[0][1]["rMax"])
    pitch = (xs[-1] - xs[0]) / (len(xs) - 1)
    meshinfo = dict(n_x=len(xs), n_y=sum("_y_" in n for n, _ in mesh), r=r,
                    pitch=pitch, z=float(mesh[0][1]["zMin"]))
    box = (xs[-1] - xs[0] + pitch, xs[-1] - xs[0] + pitch)
    return full, meshinfo, box


def g_cm2(t_mm, rho):
    return t_mm * 0.1 * rho


def build_rows(theirs, mesh):
    T = {v[0]: v for v in theirs}
    gases = [v for v in theirs if v[0] == "gas (mother)"]
    t = lambda n: (T[n][2] - T[n][1])
    d, p = 2 * mesh["r"], mesh["pitch"]
    steel_eq_theirs = 2 * math.pi * d * d / 4 / p           # mm of solid steel
    steel_eq_ours = math.pi * M.MESH_WIRE ** 2 / (2 * (M.MESH_WIRE + M.MESH_OPEN))
    front_t, back_t = gases[0][2] - gases[0][1], gases[-1][2] - gases[-1][1]
    mesh_gap = T["AmpGap"][1] - T["DriftGap"][2]
    R = []

    def row(name, th_t, th_desc, th_g, ou_t, ou_desc, ou_g, mat, same=None):
        if same is None:
            same = abs(th_t - ou_t) < 1e-6 and abs(th_g - ou_g) < 1e-3 * max(th_g, 1e-12)
        R.append(dict(name=name, th_t=th_t, th=th_desc, th_g=th_g,
                      ou_t=ou_t, ou=ou_desc, ou_g=ou_g, mat=mat, same=same))

    win = t("MylarSup")
    row("outer window", win, f"{win*1e3:.0f} µm mylar, flat", g_cm2(win, 1.40),
        M.T_WINDOW, f"{M.T_WINDOW*1e3:.0f} µm mylar, 10 mm sag (dome)",
        g_cm2(M.T_WINDOW, 1.40), "mylar", same=False)
    row("front gas", front_t, f"{front_t:.3f} mm", g_cm2(front_t, RHO_THEIRS["mmGas"]),
        M.FRONT_GAP + M.BULGE_FRONT,
        f"{M.FRONT_GAP:.3f} mm at frame, {M.FRONT_GAP+M.BULGE_FRONT:.2f} at dome centre",
        g_cm2(M.FRONT_GAP + M.BULGE_FRONT, RHO_OURS["gas"]), "gas")
    row("drift cathode", t("MylarDrift"), f"{t('MylarDrift')*1e3:.0f} µm mylar",
        g_cm2(t("MylarDrift"), 1.40), M.T_CATH_MYLAR, f"{M.T_CATH_MYLAR*1e3:.0f} µm mylar",
        g_cm2(M.T_CATH_MYLAR, 1.40), "mylar")
    row("cathode Al", t("AluHV"), f"{t('AluHV')*1e3:.0f} µm Al", g_cm2(t("AluHV"), 2.699),
        M.T_CATH_AL, f"{M.T_CATH_AL*1e3:.0f} µm Al", g_cm2(M.T_CATH_AL, 2.699), "Al")
    row("drift gas", t("DriftGap"), f"{t('DriftGap'):.3f} mm",
        g_cm2(t("DriftGap"), RHO_THEIRS["mmGas"]), M.T_DRIFT,
        f"{M.T_DRIFT:.3f} mm (4.0 frame − mesh)", g_cm2(M.T_DRIFT, RHO_OURS["gas"]), "gas")
    row("micromesh", mesh_gap,
        f"Ø{d*1e3:.0f} µm wires, pitch {p*1e3:.0f} (open {(p-d)*1e3:.0f}), x+y in 1 plane", g_cm2(steel_eq_theirs, 7.93),
        M.T_MESH, f"45/18 woven, pitch 63, {M.T_MESH*1e3:.0f} µm eff. slab",
        g_cm2(steel_eq_ours, 8.00), "steel")
    row("amp gap", t("AmpGap"), f"{t('AmpGap')*1e3:.0f} µm, no pillars",
        g_cm2(t("AmpGap"), RHO_THEIRS["mmGas"]), M.T_AMP,
        f"{M.T_AMP*1e3:.0f} µm + Dynamask pillars", g_cm2(M.T_AMP, RHO_OURS["gas"]), "gas")
    row("F.Cu (pads)", t("CopperUp"), f"{t('CopperUp')*1e3:.0f} µm full Cu",
        g_cm2(t("CopperUp"), 8.96), M.T_CU_F, f"{M.T_CU_F*1e3:.0f} µm, real pads (~{FCU_COV:.2f})",
        g_cm2(M.T_CU_F * FCU_COV, 8.96), "Cu")
    row("FR4", t("FR4"), f"{t('FR4')*1e3:.0f} µm", g_cm2(t("FR4"), 1.86),
        M.T_FR4, f"{M.T_FR4*1e3:.0f} µm", g_cm2(M.T_FR4, 1.85), "FR4")
    row("B.Cu", t("CopperDown"), f"{t('CopperDown')*1e3:.0f} µm full Cu",
        g_cm2(t("CopperDown"), 8.96), M.T_CU_B,
        f"{M.T_CU_B*1e3:.0f} µm traces (~{BCU_COV:.2f}, 0.03→0.25 radially)",
        g_cm2(M.T_CU_B * BCU_COV, 8.96), "Cu")
    row("back gas", back_t, f"{back_t:.3f} mm", g_cm2(back_t, RHO_THEIRS["mmGas"]),
        M.BACK_GAP + M.BULGE_BACK,
        f"{M.BACK_GAP:.1f} mm in carbon frame, {M.BACK_GAP+M.BULGE_BACK:.0f} at dome centre",
        g_cm2(M.BACK_GAP + M.BULGE_BACK, RHO_OURS["gas"]), "gas")
    win2 = t("MylarDown")
    row("back window", win2, f"{win2*1e3:.0f} µm mylar, flat", g_cm2(win2, 1.40),
        M.T_WINDOW, f"{M.T_WINDOW*1e3:.0f} µm mylar, 10 mm sag (dome)",
        g_cm2(M.T_WINDOW, 1.40), "mylar", same=False)
    return R


COL = dict(mylar="#a8e0a8", gas="#bfe6f5", Al="#b0b0b0", steel="#808080",
           Cu="#cc6619", FR4="#339933")


def plot(R, box):
    fig, ax = plt.subplots(figsize=(15, 0.62 * len(R) + 2.6))
    ax.set_xlim(0, 15)
    ax.set_ylim(-len(R) - 1.6, 1.3)
    ax.axis("off")
    hx = dict(name=0.1, th=2.3, ou=7.6, d=12.9)
    for k, lab in [("th", "P2Sim MMP2_v4  (Matthieu, detTest)"),
                   ("ou", "this repo  (design files)"), ("d", "areal density  theirs → ours")]:
        ax.text(hx[k], 0.55, lab, fontsize=11.5, weight="bold", va="center")
    ax.text(hx["name"], 0.55, "layer  (beam →)", fontsize=11.5, weight="bold", va="center")
    for i, r in enumerate(R):
        y = -i - 0.5
        diff = not r["same"]
        for k, desc, tt in (("th", r["th"], r["th_t"]), ("ou", r["ou"], r["ou_t"])):
            ax.add_patch(Rectangle((hx[k] - 0.05, y - 0.42), 5.15, 0.84,
                                   color=COL[r["mat"]], alpha=0.55, lw=0))
            ax.text(hx[k] + 0.05, y, desc, fontsize=9.6, va="center")
        ax.text(hx["name"], y, r["name"], fontsize=10.5, va="center",
                color="#b00000" if diff else "black", weight="bold" if diff else "normal")
        ratio = r["ou_g"] / r["th_g"] if r["th_g"] else float("nan")
        ax.text(hx["d"], y, f"{r['th_g']*1e3:7.2f} → {r['ou_g']*1e3:7.2f} mg/cm²"
                + ("" if r["same"] else f"   (×{ratio:.2f})"),
                fontsize=9.6, va="center", family="monospace",
                color="#b00000" if diff else "black")
    def tot(sel, key):
        return sum(r[key] for r in R if sel(r["name"]))
    front = lambda n: n in ("outer window", "front gas", "drift cathode", "cathode Al")
    back = lambda n: n in ("micromesh", "amp gap", "F.Cu (pads)", "FR4", "B.Cu",
                           "back gas", "back window")
    y = -len(R) - 0.4
    ax.text(hx["name"], y, "to the drift gas, from the front:", fontsize=10, weight="bold", va="center")
    ax.text(hx["d"], y, f"{tot(front,'th_g')*1e3:7.2f} → {tot(front,'ou_g')*1e3:7.2f} mg/cm²",
            fontsize=9.6, family="monospace", va="center")
    y -= 0.55
    ax.text(hx["name"], y, "to the drift gas, from behind:", fontsize=10, weight="bold", va="center")
    ax.text(hx["d"], y, f"{tot(back,'th_g')*1e3:7.2f} → {tot(back,'ou_g')*1e3:7.2f} mg/cm²",
            fontsize=9.6, family="monospace", va="center")
    y -= 0.65
    ax.text(hx["name"], y,
            f"Lateral: theirs = {box[0]:.0f}×{box[1]:.0f} mm box, no frame, no pillars, flat windows.  "
            "Ours = full wedge, plastic frame 4+4 mm on 150 µm Dynamask, 1 mm carbon back frame "
            "(board edge → pad zone), Saclay V1 pillars.  Red = differs.",
            fontsize=9.2, va="center", style="italic")
    fig.tight_layout()
    fig.savefig(OUT, dpi=140)
    print("wrote", OUT)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    theirs, mesh, box = read_theirs(d)
    print(f"P2Sim stack ({d}):")
    for n, z0, z1, m in theirs:
        print(f"  {n:14s} {z0:8.3f} {z1:8.3f}  {(z1-z0)*1e3:9.1f} µm  {m}")
    print(f"  mesh: {mesh['n_x']}+{mesh['n_y']} tubes, r {mesh['r']*1e3:.0f} µm, "
          f"pitch {mesh['pitch']*1e3:.1f} µm, at z {mesh['z']}")
    R = build_rows(theirs, mesh)
    print("\n| layer | P2Sim MMP2_v4 | this repo | mg/cm² theirs → ours |")
    print("|---|---|---|---|")
    for r in R:
        flag = "" if r["same"] else " ⚠"
        print(f"| {r['name']}{flag} | {r['th']} | {r['ou']} | "
              f"{r['th_g']*1e3:.2f} → {r['ou_g']*1e3:.2f} |")
    plot(R, box)


if __name__ == "__main__":
    main()
