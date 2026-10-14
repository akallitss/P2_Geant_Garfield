#!/usr/bin/env python3
"""Can O2 / H2O in the SPS gas explain the measured drift scan? (2026-10-14)

Model for the measured Landau MPV vs drift field E (drift_mesh_scan_1,
Ar/CO2/iso 93/5/2, mesh 450 V, 4 mm):

    M(E) = C * A(E; O2) * T(E; H2O) * S(E)

  A  survival against attachment over the drift, averaged over a uniform
     track: (1 - exp(-eta d)) / (eta d), eta(E) from the Magboltz tables of
     the contaminated gas (gas_tables_contam/*.gas exported to JSON),
     interpolated linearly in the O2 concentration (attachment ~ [O2]).
  T  mesh transparency, wire losses only (Garfield, corrected cell, null
     steps): dry from transparency_sweep_gap150.csv, humid from contam/.
  S  DREAM shaping loss: peak of a 180 ns CR-RC^4 response to a box of
     charge arrivals of width d / v_d(E), relative to a delta.
  C  free scale.

Fits C and the O2 level (ppm) for each H2O hypothesis, for P2_OUT and P2_MID,
and prints chi2. Caveat: Magboltz treats O2 three-body attachment with O2 as
the third body; with CO2 / isobutane as third bodies the real rate per ppm can
be larger, so the fitted "ppm" is an effective O2 level.

    python3 scripts/analysis/sps_contamination_fit.py <contam json> <contam transparency dir>
"""
import csv
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

REPO = pathlib.Path(__file__).resolve().parents[2]
DATA = pathlib.Path("/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/"
                    "mpgd26_workspace/products/analysis")
OUT = REPO / "docs" / "results" / "contamination"
D_CM = 0.3964
TP_NS = 180.0


def data(det):
    p = DATA / det / "drift_mesh_scan_1" / "scan_drift" / "20_beam_spectra" / "beam_mpv_vs_hv_spark_vetoed.csv"
    E, M, S = [], [], []
    for r in csv.DictReader(open(p)):
        if r["mpv_adc"]:
            E.append((float(r["hv"]) - 450.) / 0.4)
            M.append(float(r["mpv_adc"]))
            S.append(float(r["mpv_err"]))
    return np.array(E), np.array(M), np.array(S)


def table_interp(t, key, E):
    Et = np.array(t["E_V_per_cm"])
    y = np.maximum(np.array(t[key], float), 1e-12)
    return np.exp(np.interp(np.log(E), np.log(Et), np.log(y)))


def shaping(E, vd_um_ns, n=4):
    """Peak of CR-RC^n (peaking time TP_NS) for a box arrival of width d/v_d,
    relative to an instantaneous arrival."""
    tau = TP_NS / n
    t = np.arange(0, 2000, 0.5)
    h = (t / tau) ** n * np.exp(-t / tau)
    h /= h.max()
    out = []
    for v in np.atleast_1d(vd_um_ns):
        w = D_CM * 1e4 / v
        box = (t < w).astype(float)
        box /= box.sum()
        out.append(np.convolve(h, box)[: len(t)].max())
    return np.array(out)


def transparency(gas, E, contam_dir):
    if gas == "ArCO2Iso9352":
        rows = [r for r in csv.DictReader(open(REPO / "docs/results/mesh/transparency_sweep_gap150.csv"))
                if r["gas"] == gas]
        Es = np.array([float(r["E_drift"]) for r in rows])
        T = np.array([float(r["T_micro"]) + float(r["attached"]) for r in rows])  # wire losses only
    else:
        Es, T = [], []
        for p in sorted(pathlib.Path(contam_dir).glob(f"E*_{gas}_electrons.txt")):
            e = np.loadtxt(p)
            Es.append(float(p.name.split("_")[0][1:]))
            T.append(1. - ((e[:, 6] < 1) & (e[:, 2] == -5)).mean())
        Es, T = np.array(Es), np.array(T)
    o = np.argsort(Es)
    return np.interp(E, Es[o], T[o])


def main():
    tabs = json.load(open(sys.argv[1]))
    contam_dir = sys.argv[2]
    OUT.mkdir(parents=True, exist_ok=True)
    # attachment per ppm O2 from the O2 tables (linear in [O2])
    o2_tabs = {k: v for k, v in tabs.items() if "_O2_" in k and "H2O" not in k}
    Ef = np.geomspace(150, 1500, 60)
    eta_per_ppm = np.mean([table_interp(t, "eta_per_cm", Ef) / float(k.split("_O2_")[1])
                           for k, t in o2_tabs.items()], axis=0)
    base = tabs.get("ArCO2Iso9352_H2O_500") or next(iter(tabs.values()))
    vd = lambda E: table_interp(base, "v_um_per_ns", E)  # noqa: E731
    hyps = {"dry": "ArCO2Iso9352", "H2O 500 ppm": "ArCO2Iso9352_H2O_500",
            "H2O 2000 ppm": "ArCO2Iso9352_H2O_2000"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    results = []
    for ax, det in zip(axes, ("P2_OUT", "P2_MID")):
        E, M, S = data(det)
        ax.errorbar(E, M, S, fmt="o", color="k", label=f"{det} measured MPV")
        for hlab, gas in hyps.items():
            T = transparency(gas, E, contam_dir)
            Sh = shaping(E, vd(E))

            def model(ppm, E=E, T=T, Sh=Sh):
                eta = ppm * np.interp(E, Ef, eta_per_ppm)
                A = np.where(eta > 0, (1 - np.exp(-eta * D_CM)) / np.maximum(eta * D_CM, 1e-12), 1.)
                return A * T * Sh

            def chi2(ppm):
                m = model(ppm)
                C = np.sum(M * m / S ** 2) / np.sum(m * m / S ** 2)
                return np.sum(((M - C * m) / S) ** 2)

            r = minimize_scalar(chi2, bounds=(0, 20000), method="bounded")
            m = model(r.x)
            C = np.sum(M * m / S ** 2) / np.sum(m * m / S ** 2)
            ee = np.linspace(E.min(), E.max(), 100)
            mm = model(r.x, ee, transparency(gas, ee, contam_dir), shaping(ee, vd(ee)))
            ax.plot(ee, C * mm, label=f"{hlab}: O2 {r.x:.0f} ppm-eq, χ²={r.fun:.0f}/{len(E) - 2}")
            results.append({"det": det, "h2o": hlab, "o2_ppm_eq": round(r.x, 1),
                            "chi2": round(r.fun, 1), "ndf": len(E) - 2, "scale": round(C, 2)})
            print(results[-1])
        ax.set_xlabel("drift field [V/cm]")
        ax.set_ylabel("Landau MPV [ADC]")
        ax.set_title(f"{det}, Ar/CO2/iso 93/5/2: attachment × transparency × shaping", fontsize=10)
        ax.legend(fontsize=8)
    fig.savefig(OUT / "fig_contamination_fit.png", dpi=160)
    with open(OUT / "contamination_fit.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0]))
        w.writeheader()
        w.writerows(results)


if __name__ == "__main__":
    main()
