#!/usr/bin/env python3
"""How much of the ionization in detected events starts close to the mesh?

Plan item 4 (2026-10-14). The Garfield transparency T is computed for
electrons that arrive from far above the mesh (diffusion has decorrelated them
from their start point). Electrons born within ~50-100 um of the mesh may see
a different T. This counts, per Stage A file, the share of drift-gas electrons
(nPrimary) born within 25 / 50 / 100 / 200 um of the mesh, in events with at
least one drift electron, split by where the event's photon converted (gas /
mesh / other solid). AmpGas electrons (already past the mesh) are reported
separately.

    python3 scripts/analysis/near_mesh_ionization.py <dir with p2_*.root> [out.csv]

Run on lxplus (LCG view: uproot, numpy). Reads stage_b.frame / clusters.
"""
import csv
import pathlib
import re
import sys

import numpy as np
import uproot

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from stage_b.clusters import read_clusters  # noqa: E402
from stage_b.frame import read_frame  # noqa: E402

CUTS_UM = [25, 50, 100, 200]
GAS = ("DriftGas", "AmpGas")


def conv_class(v):
    v = v.decode() if isinstance(v, bytes) else str(v)
    if v in GAS:
        return "gas"
    if "mesh" in v.lower():
        return "mesh"
    return "other"


def one(path):
    fr = read_frame(path)
    cl = read_clusters(path)
    with uproot.open(path) as f:
        ev = f["EventTree"].arrays(["eventID", "convVolume"], library="np")
    order = np.argsort(ev["eventID"])
    ev_sorted = ev["eventID"][order]
    ev_cls = np.array([conv_class(v) for v in ev["convVolume"]])[order]
    depth_um = 1e3 * fr.drift_sign * (fr.mesh_z_mm - cl["z"])
    drift = cl["volume"] == "DriftGas"
    amp = cl["volume"] == "AmpGas"
    rows = []
    ev_ids = cl["eventID"]
    idx = np.clip(np.searchsorted(ev_sorted, ev_ids), 0, len(ev_sorted) - 1)
    cl_cls = np.where(ev_sorted[idx] == ev_ids, ev_cls[idx], "other")
    for group in ("all", "gas", "mesh", "other"):
        if group == "all":
            sel_ev = np.ones(len(ev_ids), bool)
        else:
            sel_ev = cl_cls == group
        nd = cl["nPrimary"][drift & sel_ev].sum()
        na = cl["nPrimary"][amp & sel_ev].sum()
        n_events = len(np.unique(ev_ids[drift & sel_ev]))
        if nd == 0:
            continue
        r = {"group": group, "events": n_events, "n_drift": int(nd),
             "amp_share": round(na / (nd + na), 5)}
        for c in CUTS_UM:
            near = drift & sel_ev & (depth_um < c)
            r[f"share_lt{c}um"] = round(cl["nPrimary"][near].sum() / nd, 5)
        rows.append(r)
    return rows


def main():
    d = pathlib.Path(sys.argv[1])
    out = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else pathlib.Path("near_mesh_ionization.csv")
    allrows = []
    for p in sorted(d.glob("p2_*_gamma*keV_th0_t0.root")) + sorted(d.glob("p2_*_electron100MeV_th0_t0.root")):
        m = re.match(r"p2_(\w+?)_(\w+?mm)_(\w+?)_th0", p.name)
        try:
            rows = one(p)
        except Exception as e:  # noqa: BLE001
            print(f"skip {p.name}: {e}")
            continue
        for r in rows:
            r = {"gas": m.group(1), "gap": m.group(2), "particle": m.group(3), **r}
            allrows.append(r)
            if r["group"] in ("all", "mesh"):
                print(f"{r['gas']:10s} {r['gap']:8s} {r['particle']:14s} {r['group']:5s} "
                      f"ev={r['events']:7d} <100um={100 * r['share_lt100um']:.2f}% "
                      f"<50um={100 * r['share_lt50um']:.2f}% amp={100 * r['amp_share']:.2f}%")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(allrows[0]))
        w.writeheader()
        w.writerows(allrows)
    print("wrote", out)


if __name__ == "__main__":
    main()
