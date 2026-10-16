#!/usr/bin/env python3
"""Primary ionization: Heed vs Geant4 on a bare slab and in the full P2 stack,
option4 vs PAI in the gas (2026-10-16; answers Dylan's MX17 handoff C1-C4).

Inputs (all from the lxplus work dir /afs/cern.ch/work/a/akallits/p2_pai):
  out/<tag>_s<k>_t0.root   mm_sim -m slab | -m p2, see jobs.txt there
                           (tag prefix np_ = --no-pixe, jobs2.txt)
  logs/<tag>_s<k>.out      mm_sim stdout + "wall_s <seconds>"
  heed/heed_<case>.txt     scripts/garfield/heed_mip.cc (col 2 = n_electrons)

Per case: mean / median / MPV of nPrimDrift (events with no drift ionization
added back as zeros from RunMeta.thrown), the MPV of the drift energy deposit,
ratios to Heed, CPU per event, and the largest per-volume step count of an
event (a stuck low-energy track, C4). For the full stack, nPrimDrift is split
by where the ionizing track was born (slab too) (ClusterTree.birthVolume): the beam
particle itself, tracks born in the gas, and tracks born in each solid.

    python3 scripts/analysis/heed_vs_g4_slab.py <work dir> <outdir>
"""
import collections
import csv
import pathlib
import re
import sys

import numpy as np
import uproot

HEED = {"slab_ArIso_4mm": "heed_ArIso_4mm", "slab_ArCF4Iso_4mm": "heed_ArCF4Iso_4mm",
        "slab_NeIso_4mm": "heed_NeIso_4mm", "slab_ArIso_30mm": "heed_ArIso_30mm",
        "slab_ArIso_4mm_mu": "heed_ArIso_4mm_mu", "p2_ArIso": "heed_ArIso_4mm"}
GAS_VOLS = {"DriftGas", "AmpGas"}


def mpv(x, width):
    """Mode of x on bins of the given width (bin centre)."""
    if len(x) == 0:
        return float("nan")
    h, e = np.histogram(x, np.arange(0, np.percentile(x, 99) + width, width))
    return float(0.5 * (e[h.argmax()] + e[h.argmax() + 1]))


def load(work, tag):
    files = sorted((work / "out").glob(f"{tag}_s*_t0.root"))
    n, ed, thrown, maxsteps = [], [], 0, 0
    clus = collections.Counter()
    for fn in files:
        with uproot.open(fn) as f:
            th = int(f["RunMeta"]["thrown"].array(library="np").sum())
            ev = f["EventTree"].arrays(["nPrimDrift", "edepDrift"], library="np")
            n.append(np.concatenate([ev["nPrimDrift"], np.zeros(th - len(ev["nPrimDrift"]), int)]))
            ed.append(ev["edepDrift"])
            thrown += th
            if "VolumeTree" in f:
                ns = f["VolumeTree"]["nSteps"].array(library="np")
                maxsteps = max(maxsteps, int(ns.max()) if len(ns) else 0)
            c = f["ClusterTree"].arrays(["volume", "trackID", "nPrimary", "birthVolume"],
                                        library="np")
            drift = c["volume"] == "DriftGas"
            for tid, bv, k in zip(c["trackID"][drift], c["birthVolume"][drift],
                                  c["nPrimary"][drift]):
                key = "beam particle" if tid == 1 else (
                    "born in gas" if bv in GAS_VOLS else f"born in {bv}")
                clus[key] += int(k)
    wall = []
    for lg in sorted((work / "logs").glob(f"{tag}_s*.out")):
        m = re.search(r"wall_s (\d+)", lg.read_text())
        if m:
            wall.append(int(m.group(1)))
    return (np.concatenate(n), np.concatenate(ed), thrown, maxsteps, clus,
            sum(wall) / thrown * 1e3 if wall and thrown else float("nan"))


def heed(work, name):
    p = work / "heed" / f"{name}.txt"
    return np.loadtxt(p, usecols=1).astype(int) if p.exists() else None


def main():
    work, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    tags = sorted({re.sub(r"_s\d+_t0\.root$", "", p.name) for p in (work / "out").glob("*_s*_t0.root")})
    rows, split = [], {}
    for tag in tags:
        base = re.sub(r"_(opt4|pai)_st\d+$", "", re.sub(r"^np_", "", tag))
        n, ed, thrown, maxsteps, clus, ms_ev = load(work, tag)
        h = heed(work, HEED.get(base, ""))
        r = dict(case=tag, N=thrown, mean=round(n.mean(), 2), median=float(np.median(n)),
                 mpv=mpv(n[n > 0], 2), edep_mpv_eV=round(mpv(ed[ed > 0], 20)),
                 ms_per_event=round(ms_ev, 2), max_steps_one_volume=maxsteps)
        if h is not None:
            r.update(heed_mean=round(h.mean(), 2), heed_median=float(np.median(h)),
                     ratio_mean=round(n.mean() / h.mean(), 3),
                     ratio_median=round(np.median(n) / np.median(h), 3))
        rows.append(r)
        if clus:
            tot = sum(clus.values())
            split[tag] = {k: round(v / thrown, 2) for k, v in clus.most_common()}
            split[tag]["total per event"] = round(tot / thrown, 2)
        print(r)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out / "heed_vs_g4_slab.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    with open(out / "p2_stack_origin_split.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "origin", "primaries_per_event"])
        for tag, d in split.items():
            for k, v in d.items():
                w.writerow([tag, k, v])
                print(f"{tag:24s} {k:28s} {v}")


if __name__ == "__main__":
    main()
