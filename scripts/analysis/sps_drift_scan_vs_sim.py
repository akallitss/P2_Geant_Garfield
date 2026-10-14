#!/usr/bin/env python3
"""SPS July 2026 drift scan (DREAM) vs the simulated mesh transparency
(plan item 5, 2026-10-14).

GAS: drift_mesh_scan_1 ran with Ar/CO2/iso 93/5/2. Its run_config.json
(mpgd26_workspace/eos_inventory/runs/) was corrected 2026-08-16 ("DAQ config
wrongly said Ar/Iso 95/5"); the beam-test logbook text still says 95/5. Always
take the gas from the corrected run_config, and simulate that gas.

Data: drift_mesh_scan_1 (25 July), mesh 450 V, drift 500-900 V, 4 mm gap,
Landau MPV of the cluster charge per sub-run from stage 20 of
P2_basket_analysis/sps_beam_analysis (beam_mpv_vs_hv_spark_vetoed.csv).
Drift field = (V_drift - 450 V) / 0.4 cm. Both curves normalised to their
maximum (the measurement has no absolute transparency).
Simulation: docs/results/mesh/transparency_sweep_gap150.csv (corrected cell,
mesh -450 V, null-collision steps), the run's gas (SIM_GAS).

    python3 scripts/analysis/sps_drift_scan_vs_sim.py
"""
import csv
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
RUN = "drift_mesh_scan_1"
SIM_GAS = "ArCO2Iso9352"   # from the corrected run_config: Ar/CO2/Iso 93/5/2
GAS_LABEL = "Ar/CO2/iso 93/5/2"
DATA = pathlib.Path("/local/home/ak271430/Documents/PostDocSaclay/data/SPS_Beam_Test/"
                    "mpgd26_workspace/products/analysis")
OUT = REPO / "docs" / "results" / "mesh"


def data(det):
    p = DATA / det / "drift_mesh_scan_1" / "scan_drift" / "20_beam_spectra" / "beam_mpv_vs_hv_spark_vetoed.csv"
    E, M, S = [], [], []
    for r in csv.DictReader(open(p)):
        if not r["mpv_adc"]:
            continue
        E.append((float(r["hv"]) - 450.) / 0.4)
        M.append(float(r["mpv_adc"]))
        S.append(float(r["mpv_err"]))
    return np.array(E), np.array(M), np.array(S)


def sim(gas="ArIso"):
    p = OUT / "transparency_sweep_gap150.csv"
    E, T = [], []
    for r in csv.DictReader(open(p)):
        if r["gas"] == gas:
            E.append(float(r["E_drift"]))
            T.append(float(r["T_micro"]))
    o = np.argsort(E)
    return np.array(E)[o], np.array(T)[o]


def main():
    fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    for det, c in [("P2_OUT", "#2F6DB5"), ("P2_MID", "#9DBFE3")]:
        E, M, S = data(det)
        ax.errorbar(E, M / M.max(), S / M.max(), fmt="o", color=c, ms=6,
                    label=f"SPS {det}: Landau MPV (DREAM), mesh 450 V")
    Es, Ts = sim(SIM_GAS)
    ax.plot(Es, Ts / Ts.max(), "-s", color="#D62728", label=f"Garfield T, {GAS_LABEL} (corrected cell, null steps)")
    ax.axvline(625, color="#888888", ls="--", lw=1)
    ax.text(630, 0.52, "DREAM nominal\n700/450 V", fontsize=8, color="#555555")
    ax.axvline(750, color="#BBBBBB", ls=":", lw=1)
    ax.text(755, 0.45, "VMM nominal\n750/450 V", fontsize=8, color="#555555")
    ax.set_xlim(0, 2600)
    ax.set_ylim(0.3, 1.08)
    ax.set_xlabel("drift field [V/cm]")
    ax.set_ylabel("relative to the maximum")
    ax.set_title(f"{GAS_LABEL}, 4 mm ({RUN}): measured signal vs simulated mesh transparency", fontsize=11)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3)
    fig.savefig(OUT / "fig_sps_drift_scan_vs_sim.png", dpi=160)
    E, M, S = data("P2_OUT")
    print("E   data/max   sim/max")
    for e, m in zip(E, M / M.max()):
        print(f"{e:6.0f} {m:7.3f} {np.interp(e, Es, Ts) / Ts.max():7.3f}")


if __name__ == "__main__":
    main()
