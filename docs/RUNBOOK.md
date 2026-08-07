# Runbook — running this repo on lxplus

**Written:** 2026-08-07, for the handover to Alexandra (CEA Saclay).
**Read `HANDOFF.md` first** for what the simulation *is*; this document is
only about running it.

---

## 0. Read this before planning any work

The stated plan is "fork the repo, plug in the pillar file, run the full
chain". Two parts of that sentence do not yet correspond to anything in the
repo, and it is better to know now:

### ⛔ "The full chain" does not exist yet — only Stage A does

The campaign is designed as three stages (`SIM_CAMPAIGN_PLAN.md` §2):

| stage | what it does | state |
|---|---|---|
| **A** — Geant4 | beam/background → ionization clusters in the gas | ✅ **works, this is what you can run today** |
| **B** — drift/avalanche/induction → VMM | clusters → per-pad charge and time | ⛔ **not written.** P0.9. The pad-map loader and a standalone VMM emulator exist (`vmm/`), but nothing joins Stage A output to them |
| **C** — electronics/analysis | PDO/TDO, thresholds, clustering, observables | 🟡 partial — the VMM emulator (`vmm/`) and the timing toy (`vmm/time_resolution.py`) run standalone, on synthetic input, not on Stage A output |

So today you can produce ionization/edep/conversion-budget physics — which is
what the gas decision actually rests on (`PHOTON_DISCRIMINATION_NOTES.md`) —
but **not** pad charge, pad multiplicity, efficiency-vs-threshold, or a
time-resolution number derived from real tracks.

Before writing Stage B from scratch, read `HANDOFF_MX17_RESPONSE.md` §2.2:
MX17 has a working digitizer (`MX17_Geant/response/digitizer/`, ~1700 lines
with selftests, as of 2026-08-07) whose decomposition is identical to ours up
to two plug-ins. **The plan of record is to lift it, not to rewrite it.**

### ⛔ There is no pillar code to plug a file into

`grep -ri pillar src/ include/` returns nothing. The amplification gap is
pure gas. P0.16 is not "obtain a file and drop it in" — it is **obtain the
file *and* write the geometry code that consumes it**: place pillars as
daughters of `AmpGas` (as `MX17_Geant/shared/MX17ModuleGeometry.hh` does),
apply the ~4.8 % effective-gas correction, and teach Stage B that electrons
landing on a pillar are lost.

What we already have: the design pattern from the CERN bulk mask
(`design/gerbers/bulk_masks_CERN/P2_Mask2.gbr` — Ø 0.5 mm, pitch 2.000 mm,
41 366 pillars, 4.8 % of the gap) and the material (Dynamask dry film,
ρ ≈ 1.2–1.4 g/cm³). What is missing is the map the *analysis* uses, so sim
and data share one definition.

### ⛔ HTCondor submission does not support `p2` mode

`scripts/submit_condor*.py` still speak the inherited MX17 modes. They pass
`-g/-p/-e/-n/-o/-s` and `-a` (Al thickness, a vacuum-mode flag) and **no
`-m`, no `--drift-gap`, no `--gun-*`, no `--skip-empty`**. Because `p2` is
the *default* mode, they will happily run — producing p2-geometry jobs at
default gaps, with no angle control, tagged with MX17-style names that
`collect_results.py` then parses. That is a silent-wrong-answer path, not a
crash. **P0.7 must be done before any campaign submission.**

Everything below therefore describes **interactive / single-node Stage A
running**, which does work.

---

## 1. Environment

```bash
ssh lxplus
git clone <your fork> p2_geant && cd p2_geant
source scripts/setup_lxplus.sh     # Geant4 11.2 + ROOT via LCG, from CVMFS
bash scripts/build.sh              # -> build/mm_sim
```

`setup_lxplus.sh` probes for the newest available Geant4 and LCG view and
falls back gracefully; it prints what it picked. Verified locally against
Geant4 11.2.0 / ROOT 6.36.

**If you will also run Garfield++ (Stage B, gas tables, field solves): do
not use the CVMFS Garfield.** See `research/TOOLCHAIN_NOTES.md` §1 — the LCG
builds are 281–664 commits behind and are missing, among other things, the
neBEM OpenMP race fix and the interface-crossing check that makes mesh
transparency correct. Magboltz gas tables are unaffected either way.

## 2. Smoke test

```bash
./build/mm_sim -m p2 -g ArIso -p electron -e 119 -n 500 -t 4 -o smoke.root
```

Expect: overlap checks all `OK!`, the geometry banner, then per-thread files
`smoke.root_t{0..3}.root`. Roughly 500 electron events per thread-minute.

```bash
./build/mm_sim --list-gases          # the mixtures, campaign ones starred
./build/mm_sim -h                    # all flags
```

## 3. Things that will bite you

**Output is per worker thread, not one file.** `-t 4` gives four files.
`hadd` them, or use `scripts/collect_results.py`, which does it for you.

**Never divide by `GetEntries()`.** With `--skip-empty` (which you want for
photon runs — they are ≥99.9 % empty) the trees hold only events that made a
cluster. The denominator is `RunMeta.thrown`, summed over worker rows.
`OUTPUT_FORMAT.md` §0. At 60 keV the two differ by ~300×.

**Check `RunMeta.gitDirty`.** Every output file records the commit it came
from and whether the tree was dirty. A production run from a dirty tree is
not reproducible; `collect_results.py` warns, but only if you read the warning.

**Do not pool files whose `geometryHash` differs.** It hashes every
geometry-affecting parameter. `collect_results.py` warns on a mixed merge.

**A pencil beam on a pad boundary silently biases pad-level observables.**
The default aim point is now a pad *ring centre* (r = 349.286 mm, φ = 30°),
and `mm_sim` warns if you move it near a boundary. But a single point is
still a single point: **for any pad-multiplicity, charge-sharing or
positional observable, pass `--beam-spread 11.43`** (one ring pitch) so the
impact point averages over the pad cell. MX17 lost a first result to exactly
this trap — a spurious 3× left/right asymmetry that was entirely the beam
position (`MX17_Geant/design/RESPONSE_SIM_PLAN.md` §7).

**`--gun-theta` pivots about the drift mid-plane**, not the window, so an
angle scan keeps illuminating the same pads. The standoff auto-raises to
clear the window bulge; you will see a line about it.

**Physics gotcha, already handled but worth knowing:** the γ production cut
used to suppress Cu-K fluorescence at 8.05 keV — the dominant gas-sensitive
wall channel. Fixed with `SetDeexcitationIgnoreCut(true)` plus a 1 µm
fine-cut region. If you change the physics list, re-check that the Cu-K,
Fe-K and Cr-K lines are still present.

## 4. What is worth running today

These are the Stage-A-only studies that produce real results now:

- **Photon conversion-layer budget per gas** (`SIM_CAMPAIGN_PLAN.md` §4.2) —
  the decision variable for Ar vs Ne. Needs `--skip-empty` and 10⁶–10⁷
  photons per point. The first run of this already overturned the premise of
  the gas argument: walls convert 90 % of fakes in argon, 97 % in neon.
- **Electron/muon sensitivity scans** (§4.1, §4.3), now including the angle
  axis that P0.3 unblocked.
- **The 200 GeV muon test-beam anchor** (§4.3) — `-p muon -e 200000`.
  Note this is blocked on the as-run conditions in
  `testbeam/TB_CONDITIONS.md`, which is the one concrete thing being asked
  of you (`TESTBEAM_PLAN.md`).

Each of these wants P0.7 (condor) first if run at campaign statistics.

## 5. Python tooling

`scripts/collect_results.py` needs `pandas` + `uproot`; the LCG view has
both. (It will not run on a bare local Python without pandas.)
`scripts/plot_track_stats.py` additionally needs `scipy`.

`vmm/nl_map.py` and `vmm/time_resolution.py` are standalone (numpy only) and
run anywhere.

## 6. Current blocking list

Before a *campaign* — as opposed to interactive Stage A runs — these are
open. Full detail in `SIM_CAMPAIGN_PLAN.md` §3.

| | |
|---|---|
| P0.7 | condor + `collect_results.py` do not speak `p2` |
| P0.9 | Stage B does not exist (lift MX17's — check §2.2 first) |
| P0.10 | no Magboltz tables; generate wet variants alongside dry |
| P0.16 | pillars: file **and** geometry code both missing |
| P0.18 | `p2_fcu_coverage = 0.983` is inflated — the gerber script mis-reads `RotRect` apertures. Do not quote it |
| P0.4b/P0.4c | segment + exit-particle exports |
| P0.6 | PAI model in the gas region |
| — | Alexandra's geometry answers (`HANDOFF.md` §3) still gate production runs |
| — | the phase-space file (P0.15) — the highest-value external ask we have |
