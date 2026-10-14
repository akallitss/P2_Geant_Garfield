#!/usr/bin/env python3
"""Drift velocity vs transverse diffusion of the P2 gases (plan item 5).

From the Stage B Magboltz JSON: one curve per gas, parametric in the drift
field (100-1500 V/cm), with the fields we care about marked. Transverse
diffusion drives the mesh losses, drift velocity the timing. Reference points
read off the ILC-TPC Magboltz plots (A. Ishikawa,
http://david.attie.free.fr/ILC-TPC/gas/, B = 0, read by eye, +-5 %) at
625 V/cm are overlaid as a cross-check of our tables.

    python3 scripts/analysis/gas_vd_dt_map.py
"""
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
JS = REPO / "stage_b" / "gas_tables" / "magboltz_dry_20261008.json"
OUT = REPO / "docs" / "results" / "gas"
LABEL = {"ArIso": "Ar/iso 95/5", "ArIso9010": "Ar/iso 90/10", "NeIso": "Ne/iso 95/5",
         "NeIso8515": "Ne/iso 85/15", "ArCO2Iso9352": "Ar/CO2/iso 93/5/2",
         "ArCF4Iso": "Ar/CF4/iso 88/10/2"}
COL = {"ArIso": "#D62728", "ArIso9010": "#FF9F80", "NeIso": "#2CA02C",
       "NeIso8515": "#17BECF", "ArCO2Iso9352": "#9467BD", "ArCF4Iso": "#8C564B"}
MARK_E = {313: "s", 450: "^", 625: "o", 750: "D", 1000: "v"}
# ILC-TPC (Magboltz, A. Ishikawa) at ~625 V/cm, read off the log plots
ILC = {"ArIso": (480, 33), "ArIso9010": (350, 46), "NeIso": (410, 33)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tab = json.load(open(JS))
    fig, ax = plt.subplots(figsize=(9, 6.2), constrained_layout=True)
    rows = []
    for g, lab in LABEL.items():
        t = tab[g]
        E = np.array(t["E_V_per_cm"])
        lE = np.log(E)
        ee = np.geomspace(100, 1500, 120)
        dt = np.exp(np.interp(np.log(ee), lE, np.log(t["DT_um_sqrtcm"])))
        vd = np.exp(np.interp(np.log(ee), lE, np.log(t["v_um_per_ns"])))
        ax.plot(dt, vd, "-", color=COL[g], lw=1.6, label=lab)
        for e0, mk in MARK_E.items():
            d0 = np.exp(np.interp(np.log(e0), lE, np.log(t["DT_um_sqrtcm"])))
            v0 = np.exp(np.interp(np.log(e0), lE, np.log(t["v_um_per_ns"])))
            ax.plot(d0, v0, mk, color=COL[g], ms=7, mec="k", mew=0.5)
            rows.append((g, e0, d0, v0))
        if g in ILC:
            ax.plot(*ILC[g], "*", color=COL[g], ms=16, mec="k", mew=0.8)
    for e0, mk in MARK_E.items():
        gap = {313: "8 mm", 450: "(sweep)", 625: "4 mm", 750: "4 mm, VMM run", 1000: ""}[e0]
        ax.plot([], [], mk, color="#999999", mec="k", label=f"{e0} V/cm {gap}".strip())
    ax.plot([], [], "*", color="#999999", mec="k", ms=12, label="ILC-TPC Magboltz (read off, 625 V/cm)")
    ax.set_xlabel("transverse diffusion D_T [µm/√cm]   (→ more mesh losses, wider clusters)")
    ax.set_ylabel("drift velocity [µm/ns]   (→ faster signal)")
    ax.set_title("P2 gases, 1 atm, 20 °C: drift velocity vs transverse diffusion, 100–1500 V/cm",
                 fontsize=11)
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    ax.grid(alpha=0.3)
    fig.savefig(OUT / "fig_vd_vs_dt.png", dpi=160)
    with open(OUT / "vd_dt_points.csv", "w") as f:
        f.write("gas,E_V_per_cm,DT_um_sqrtcm,v_um_per_ns\n")
        for g, e0, d0, v0 in rows:
            f.write(f"{g},{e0},{d0:.1f},{v0:.2f}\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
