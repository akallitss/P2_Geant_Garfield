#!/usr/bin/env python3
"""Does the drift-scan transparency edge move with the local gain? (2026-10-14)

Test of the geometry hypothesis on the cosmic bench, det1 (P2_1), Ar/iso 95/5,
4 mm drift gap. If the gain non-uniformity of det1 (stage 23: smooth saddle,
87 % peak-to-peak in pad median amplitude) comes from the amplification gap,
the high-gain pads have a smaller gap / higher E_amp and their mesh
transparency must stay on the plateau up to a HIGHER drift field than the
low-gain pads (the edge scales with E_drift / E_amp). No shift between the
groups argues against a local gap variation as the cause of the gain map.

Gain map: per-pad median amplitude from the 24 h long run at the same mesh
voltage (det1_long5, mesh 415 / drift 615 V), i.e. independent of the scan.
Pads are split into gain terciles. Scan: det1_driftscan2 (mesh 415 V, drift
415-965 V, 30 min per point). Per point and group: median and 75th percentile
of the pad-hit amplitude above a signal threshold, each normalised to that
group's maximum over the scan. The threshold truncates the low-gain group
more, so the result is shown for two thresholds.

Uses the cosmic-bench analysis modules (P2_basket_analysis/cosmic_bench_analysis)
for mapping, dead connectors and hot pads, as stages 16 and 23 do.

    python3 scripts/analysis/cosmic_drift_scan_regions.py
"""
import os
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COSMIC = "/local/home/ak271430/Documents/PostDocSaclay/P2_basket_analysis/cosmic_bench_analysis"
sys.path.insert(0, COSMIC)
import p2_qa_config as qa  # noqa: E402
qa.setup_paths()
import p2_io as p2io  # noqa: E402
import p2_mapping as pmap  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = REPO / "docs" / "results" / "cosmic"
AMP_MAX, AMP_BIN = 5000.0, 10.0
THRESHOLDS = (150.0, 300.0)
MESH_V = 415.0
PLATEAU = (250, 750)     # V/cm: reference window, flat in every group
HIGH = (1125, 1375)      # V/cm: the falling edge


def channel_table(cfg):
    return pmap.build_channel_table(
        cfg.run_config_path, cfg.MAP_CSV_PATH, det_type=cfg.DET_TYPE,
        det_name=cfg.DET_NAME, strategy="reverse",
        drop_connectors=cfg.DEAD_CONNECTORS,
        strategy_overrides=cfg.STRATEGY_OVERRIDES)


def pad_hist(hits_dir, ct, pads, min_amp, t_max_h=None):
    """Per-pad amplitude histogram (n_pads x n_bins), streamed."""
    m = ct[ct["mapped"] & ct["channel_id"].isin(pads)]
    lut = pd.Series(m["channel_id"].to_numpy(),
                    index=pd.MultiIndex.from_arrays([m["feu"].astype(int), m["channel"].astype(int)]))
    row_of = {int(p): i for i, p in enumerate(pads)}
    nb = int(round(AMP_MAX / AMP_BIN))
    hist = np.zeros((len(pads), nb), dtype=np.int64)
    for h in p2io.iter_hits(hits_dir, ["channel", "amplitude", "feu"], feus=ct.attrs["feus"],
                            min_amp=min_amp, t_max_h=t_max_h, progress=False):
        h = h[h["amplitude"] > THRESHOLDS[0]]
        if not len(h):
            continue
        cid = lut.reindex(pd.MultiIndex.from_arrays(
            [h["feu"].astype(int), h["channel"].astype(int)])).to_numpy()
        ok = pd.notna(cid)
        if not ok.any():
            continue
        rows = np.array([row_of[int(c)] for c in cid[ok]])
        cols = np.clip((h["amplitude"].to_numpy()[ok] / AMP_BIN).astype(int), 0, nb - 1)
        np.add.at(hist, (rows, cols), 1)
    return hist


def boot_err(hist_1d, q, thr, nb=300, seed=1):
    """Bootstrap error of the quantile: multinomial resampling of the histogram."""
    c = (np.arange(len(hist_1d)) + 0.5) * AMP_BIN
    h = np.where(c > thr, hist_1d, 0)
    n = int(h.sum())
    if n < 30:
        return np.nan
    rng = np.random.default_rng(seed)
    qs = [c[(np.cumsum(rng.multinomial(n, h / n)) >= q * n).argmax()] for _ in range(nb)]
    return float(np.std(qs))


def quantile(hist_1d, q, thr):
    c = (np.arange(len(hist_1d)) + 0.5) * AMP_BIN
    h = np.where(c > thr, hist_1d, 0)
    n = h.sum()
    if n < 30:
        return np.nan, n
    return c[(np.cumsum(h) >= q * n).argmax()], n


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    long_cfg, scan_cfg = qa.RUNS["det1_long5"], qa.RUNS["det1_driftscan2"]
    ct = channel_table(scan_cfg)
    long_hits = os.path.join(long_cfg.subrun_dir(long_cfg.SUB_RUN), "combined_hits_root")
    drop = set(p2io.drop_pads_for(long_cfg, channel_table(long_cfg), hits_dir=long_hits))
    pads = np.sort([p for p in ct[ct["mapped"]]["channel_id"].unique() if p not in drop])

    # gain map from the long run (first 6 h are plenty: ~10^5 hits)
    hl = pad_hist(long_hits, ct, pads, long_cfg.MIN_AMP, t_max_h=6.0)
    med = np.array([quantile(r, 0.5, 300.0)[0] for r in hl])
    nl = hl.sum(axis=1)
    good = np.isfinite(med) & (nl >= 50)
    pads, med = pads[good], med[good]
    t1, t2 = np.nanpercentile(med, [33.3, 66.7])
    group = np.where(med < t1, 0, np.where(med < t2, 1, 2))
    names = [f"low gain (median < {t1:.0f})", f"mid gain", f"high gain (median > {t2:.0f})"]
    print(f"{len(pads)} pads; gain terciles at {t1:.0f} / {t2:.0f} ADC; "
          f"group medians {[np.median(med[group == g]) for g in range(3)]}")
    geo = ct[ct["mapped"]].drop_duplicates("channel_id").set_index("channel_id")
    pd.DataFrame({"channel_id": pads, "long_median_adc": med, "group": group,
                  "pad_cx": geo.loc[pads, "pad_cx"].to_numpy(),
                  "pad_cy": geo.loc[pads, "pad_cy"].to_numpy()}).to_csv(OUT / "det1_gain_groups.csv", index=False)

    rows = []
    for sub in sorted(os.listdir(scan_cfg.run_dir)):
        if not sub.startswith("drift_scan_det1_415_"):
            continue
        vd = float(sub.split("_")[-1])
        E = (vd - MESH_V) / (scan_cfg.DRIFT_GAP_MM / 10.0)
        hits = os.path.join(scan_cfg.subrun_dir(sub), "combined_hits_root")
        hs = pad_hist(hits, ct, pads, scan_cfg.MIN_AMP)
        for g in range(3):
            hg = hs[group == g].sum(axis=0)
            for thr in THRESHOLDS:
                q50, n = quantile(hg, 0.5, thr)
                q75, _ = quantile(hg, 0.75, thr)
                rows.append(dict(E=E, group=g, thr=thr, n=n, q50=q50, q75=q75,
                                 q50_err=boot_err(hg, 0.5, thr), q75_err=boot_err(hg, 0.75, thr)))
        print(f"  {sub}: E {E:.0f} V/cm, hits/group {[int(hs[group == g].sum()) for g in range(3)]}")
    d = pd.DataFrame(rows).sort_values(["thr", "group", "E"])
    d.to_csv(OUT / "det1_drift_scan_by_gain.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    col = ["#2F6DB5", "#888888", "#D62728"]
    for ax, (thr, q) in zip(axes, [(THRESHOLDS[0], "q50"), (THRESHOLDS[1], "q50"), (THRESHOLDS[0], "q75")]):
        for g in range(3):
            s = d[(d.thr == thr) & (d.group == g) & (d.E > 0)]
            ref = s[(s.E >= PLATEAU[0]) & (s.E <= PLATEAU[1])][q].mean()
            y = s[q] / ref
            err = s[q + "_err"] / ref
            ax.errorbar(s["E"], y, err, fmt="o-", color=col[g], ms=4, label=names[g])
        ax.set_xlabel("drift field [V/cm]")
        ax.set_ylabel(f"{'median' if q == 'q50' else '75th percentile'} pad amplitude / plateau "
                      f"({PLATEAU[0]}-{PLATEAU[1]} V/cm)")
        ax.set_title(f"hits > {thr:.0f} ADC, {q}", fontsize=10)
        ax.set_ylim(0.4, 1.1)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    # edge ratio R = <amplitude>_HIGH / <amplitude>_PLATEAU per group, with error
    summ = []
    for thr in THRESHOLDS:
        for q in ("q50", "q75"):
            for g in range(3):
                s = d[(d.thr == thr) & (d.group == g)]
                a = s[(s.E >= PLATEAU[0]) & (s.E <= PLATEAU[1])]
                b = s[(s.E >= HIGH[0]) & (s.E <= HIGH[1])]
                ma, mb = a[q].mean(), b[q].mean()
                ea = np.sqrt((a[q + "_err"] ** 2).sum()) / len(a)
                eb = np.sqrt((b[q + "_err"] ** 2).sum()) / len(b)
                R = mb / ma
                summ.append(dict(thr=thr, q=q, group=g, R=R, R_err=R * np.hypot(ea / ma, eb / mb)))
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "det1_edge_ratio_by_gain.csv", index=False)
    print(summ.round(4).to_string())
    fig.suptitle("Cosmic det1, Ar/iso 95/5, mesh 415 V: drift scan split by local gain "
                 "(gain map from the 24 h run)", fontsize=11)
    fig.savefig(OUT / "fig_det1_drift_scan_by_gain.png", dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
