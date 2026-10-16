#!/usr/bin/env python3
"""Per-chamber hit-time statistics vs drift voltage from the SPS combined hits.

Runs on lxplus (needs EOS). For drift_mesh_scan_1, each drift point and each
P2 chamber (FEU 3 = P2_IN, 4 = P2_MID, 5 = P2_OUT): per event, the largest
hit in the chamber with amplitude in [AMP_MIN, AMP_MAX], not saturated or
truncated. Writes the time histogram, quantiles and the per-event MID/OUT
time minus P2_IN time (P2_IN HV is fixed, so that difference cancels the
trigger). Times are the DREAM 30 % leading edge, ftst-corrected [ns].

    python3 sps_hit_times_extract.py <run dir> <out.json> [max hits per file]
"""
import glob
import json
import sys

import numpy as np
import uproot

FEU = {3: "P2_IN", 4: "P2_MID", 5: "P2_OUT"}
POINTS = ["drift_450", "drift_500", "drift_550", "drift_600", "drift_650",
          "drift_700", "nominal_00", "drift_750", "drift_800", "drift_850",
          "drift_900"]
AMP_MIN, AMP_MAX = 200.0, 3500.0
EDGES = np.arange(0, 961, 10.0)


def leading_hits(files, max_hits):
    """{feu: (eventId, time, amp)} for the largest good hit per event."""
    br = ["eventId", "feu", "amplitude", "time", "saturated", "trunc_left",
          "trunc_right"]
    keep = {f: [] for f in FEU}
    for fn in files:
        n = 0
        for a in uproot.open(fn)["hits"].iterate(br, step_size=2_000_000, library="np"):
            ok = ((a["amplitude"] >= AMP_MIN) & (a["amplitude"] <= AMP_MAX)
                  & ~a["saturated"] & ~a["trunc_left"] & ~a["trunc_right"])
            for f in FEU:
                s = ok & (a["feu"] == f)
                keep[f].append(np.stack([a["eventId"][s].astype(float),
                                         a["time"][s], a["amplitude"][s]], 1))
            n += len(a["eventId"])
            if max_hits and n >= max_hits:
                break
    out = {}
    for f, parts in keep.items():
        x = np.concatenate(parts) if parts else np.zeros((0, 3))
        # largest amplitude per event
        o = np.lexsort((-x[:, 2], x[:, 0]))
        x = x[o]
        first = np.r_[True, x[1:, 0] != x[:-1, 0]]
        out[f] = x[first]
    return out


def q(x):
    if len(x) < 50:
        return None
    p = np.percentile(x, [5, 16, 25, 50, 75, 84, 95])
    return dict(zip(["p05", "p16", "p25", "p50", "p75", "p84", "p95"], p.tolist()))


def main():
    run, out = sys.argv[1], sys.argv[2]
    max_hits = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    res = {}
    for pt in POINTS:
        files = sorted(glob.glob(f"{run}/{pt}/combined_hits_root/*.root"))
        if not files:
            continue
        lh = leading_hits(files, max_hits)
        r = {}
        ref = lh[3]
        for f, name in FEU.items():
            x = lh[f]
            d = dict(n=int(len(x)), time=q(x[:, 1]), amp=q(x[:, 2]),
                     hist=np.histogram(x[:, 1], EDGES)[0].tolist())
            if f != 3:
                common, ia, ib = np.intersect1d(x[:, 0], ref[:, 0], return_indices=True)
                dt = x[ia, 1] - ref[ib, 1]
                d["dt_vs_IN"] = q(dt)
                d["n_pairs"] = int(len(dt))
            r[name] = d
        if 4 in lh and 5 in lh:
            common, ia, ib = np.intersect1d(lh[4][:, 0], lh[5][:, 0], return_indices=True)
            r["MID_minus_OUT"] = q(lh[4][ia, 1] - lh[5][ib, 1])
        res[pt] = r
        print(pt, {n: (r[n]["n"], r[n]["time"] and round(r[n]["time"]["p50"], 1))
                   for n in FEU.values()}, flush=True)
    res["_edges"] = EDGES.tolist()
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
