# Measured speed and campaign cost

**Written:** 2026-08-05, on `lxplus925` (16 cores, AlmaLinux 9, Geant4
11.2.p01 + LCG_106 ROOT from CVMFS). This is the **first time the P2
simulation has run at all**, so these are also the first numbers on whether
the campaign is affordable. Short answer: **it is, comfortably — compute is
not the binding constraint on this campaign.**

Reproduce with **`bash scripts/benchmark_lxplus.sh`** (or `quick` for just
the cost-model fit).

> **Methodology note.** lxplus login nodes are shared and were running at
> load 9–15 throughout, so **wall-clock timings there are meaningless.** All
> per-event costs below are **user CPU time** (`/usr/bin/time %U`), which is
> only mildly affected by contention. An early run of this benchmark was
> discarded because two instances were accidentally running at once.

---

## 1. The cost model

Every job pays a **fixed startup cost** (geometry construction, physics
table building) plus a marginal per-event cost. Fitted from event-count
scans at 119 MeV e⁻ (200 / 1 000 / 5 000 / 20 000) and 60 keV γ
(5 000 / 50 000), single-threaded:

| | startup | marginal | throughput |
|---|---|---|---|
| **γ, 60 keV** | 6.6 s | 0.50 ms/event | **≈ 2 000 γ/s/core** |
| **e⁻, 119 MeV** | 6.4 s | 1.82 ms/event | **≈ 550 e⁻/s/core** |

The fit predicts 56.5 s of CPU for 10⁵ γ (measured 56.0) and 10.0 s for
2 000 e⁻ (measured 10.4), so it is good to a few percent.

**Independently reproduced** by a second, later run of
`scripts/benchmark_lxplus.sh quick` on the same node at a different load:
0.502 ms/γ (1 990 γ/s) and 1.77 ms/e⁻ (565 e⁻/s), startup 6.5–6.9 s. That
pass added a third photon point (2×10⁵) which confirms strict linearity —
predicted 107.3 s CPU, measured 107.26 s. Take the spread between the two
passes (≈3 %) as the reproducibility on a shared node.

**The 6.5 s startup is the reason not to submit tiny Condor jobs.** At
10⁴ events a photon job is ~57 % startup overhead; at 10⁶ it is 1 %. Size
jobs at ≥10⁶ photons or ≥10⁵ electrons.

## 2. Raw measurements (single thread unless noted)

| run | N | CPU user [s] | rate [ev/s] | output | written |
|---|---|---|---|---|---|
| γ 8 keV, ArIso, skip-empty | 10⁵ | 45.6 | 2190 | 7.93 MB | 5551 |
| γ 60 keV, ArIso, skip-empty | 10⁵ | 56.0 | 1790 | 575 kB | 347 |
| γ 100 keV, ArIso, skip-empty | 10⁵ | 54.6 | 1830 | 322 kB | 172 |
| γ 150 keV, ArIso, skip-empty | 10⁵ | 54.9 | 1820 | 224 kB | 109 |
| γ 60 keV, **NeIso**, skip-empty | 10⁵ | 56.7 | 1770 | 426 kB | **301** |
| γ 60 keV, ArIso, **no** skip-empty | 10⁵ | 57.2 | 1750 | **1.89 MB** | 100000 |
| e⁻ 30 MeV, ArIso | 2000 | 11.4 | — | 2.17 MB | 2000 |
| e⁻ 119 MeV, ArIso | 2000 | 10.4 | — | 2.01 MB | 2000 |
| e⁻ 155 MeV, ArIso | 2000 | 10.0 | — | 2.00 MB | 2000 |
| e⁻ 119 MeV, NeIso | 2000 | 9.7 | — | 1.73 MB | 2000 |
| μ 200 GeV, ArIso | 2000 | 8.8 | — | 1.72 MB | 2000 |
| γ 60 keV, 8 threads | 4×10⁵ | 208.5 | 11 019 (wall) | 2.52 MB | 1384 |
| e⁻ 119 MeV, 8 threads | 8000 | 21.4 | — | 8.01 MB | 8000 |

(The short 2000-event runs are startup-dominated — do not read a throughput
off them; use §1.)

Observations worth carrying:

- **Photon cost is nearly independent of energy** across 8–150 keV. It is
  dominated by geometry navigation of the un-interacting photon, not by
  physics. So the expensive photon points are the ones with high statistics,
  not the ones at high energy.
- **8 keV is the cheapest but bulkiest** — 5.5 % of photons interact, giving
  14× the output of the 60 keV run.
- **Gas choice barely changes runtime** (Ne vs Ar within 2 %).
- **200 GeV muons are the cheapest particle we run.** The test-beam
  validation points (plan §4.3) are free.

## 3. Multithreading

8 threads on 400 k photons: 36.3 s wall against an extrapolated 243 s
single-threaded → **6.7× speedup, 84 % efficiency**, and total CPU is the
same as the single-thread run per event (no MT overhead). On a shared login
node at load ~12 that is a floor, not a ceiling.

For HTCondor the useful mode is still **1 thread per job, many jobs** —
scaling is embarrassingly parallel across run points and avoids the shared
ROOT-file merge. MT is for interactive work.

> Multithreading only works as of 2026-08-05: every `-t > 1` run used to
> core-dump in `TTree` construction because ROOT's global class registry is
> not thread-safe. Fixed with `ROOT::EnableThreadSafety()` in `mm_sim.cc`.

## 4. Output volume

| | per event | note |
|---|---|---|
| γ 60 keV with `--skip-empty` | **5.75 B/thrown** | 3.3× smaller than without |
| γ 60 keV without | 18.9 B/thrown | |
| γ 8 keV with `--skip-empty` | 79 B/thrown | high interaction probability |
| e⁻ 119 MeV | **1.0 kB/event** | dominated by ClusterTree |

The electron figure roughly doubled with the provenance rework: each cluster
now carries four `char[32]` process/volume strings. If storage ever binds,
replacing them with enum codes is the obvious fix — but at the campaign
scale below it does not bind.

## 5. What the campaign costs

Using the §1 model against the run matrices in `SIM_CAMPAIGN_PLAN.md`:

| phase | runs | events | core-hours |
|---|---|---|---|
| Phase 1 electron scan (13 gases × ~40 points × 10⁴) | 520 | 5.2×10⁶ | 3 |
| Phase 1 photon scan (13 gases × 17 E × 3 angles × 10⁶) | 663 | 6.6×10⁸ | 92 |
| Phase 1 photon, 10⁷ at the 4 critical energies | 156 | 1.6×10⁹ | 216 |
| Phase 3 matrix, electrons (7 gaps × 8 gases × 5 pts × 10⁵) | 280 | 2.8×10⁷ | 14 |
| Phase 3 matrix, photons (7 gaps × 8 gases × 4 pts × 10⁷) | 224 | 2.2×10⁹ | 310 |
| **total** | **~1 850** | **~4.5×10⁹** | **≈ 630** |

Wall-clock on HTCondor: **~13 h at 50 slots, ~6 h at 100, ~3 h at 200.**
Storage: **≈ 59 GB** (revising the plan's <50 GB estimate slightly upward).

### The conclusion that matters

**The Geant4 stage is not the bottleneck and should not be economised on.**
630 core-hours is a few hours of a modest Condor allocation. Three
consequences for planning:

1. **Statistics are cheap — take more than the minimum.** Going from 10⁶ to
   10⁷ photons at *every* Phase-1 energy (not just the four critical ones)
   costs ~800 extra core-hours and removes all statistical caveats from the
   conversion-layer budget. That is a good trade.
2. **The real costs are elsewhere**: Magboltz gas tables (~1 core-hour per
   gas per E-point — plan P0.10, and the one place to cache aggressively),
   Stage B/C development, and analysis time. Budget attention accordingly.
3. **Rerunning after a bug fix is affordable.** Given that three bugs
   surfaced on the first run alone, and the P0.1 gas-mixture fix is still
   pending, plan on re-running the whole campaign at least once and do not
   treat any single production pass as final.

## 6. Caveats

- Single machine, single node type. HTCondor workers vary; treat these as
  ±30 % for planning.
- Measured with the *current* physics configuration: EM option 4 plus
  deexcitation with `deexcitationIgnoreCut`, and a 1 µm-cut `P2FineCut`
  region over gas/mesh/pad-Cu. **The PAI model (plan P0.6) is not yet in**
  and is expected to roughly double the EM cost in the gas regions — a
  small fraction of the total for photon runs, larger for electron runs.
- The Phase-1/Phase-3 run counts are read off the plan's matrices; if the
  gas list or angle grid changes, rescale linearly.
