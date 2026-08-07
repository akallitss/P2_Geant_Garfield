# Runbook — running this repo on lxplus

**Written:** 2026-08-07, for the handover to Alexandra (CEA Saclay).
**Read `HANDOFF.md` first** for what the simulation *is*; this document is
only about running it.

---

## 0. Read this before planning any work

The stated plan is "fork the repo, plug in the pillar file, run the full
chain". Two parts of that sentence do not yet correspond to anything in the
repo, and it is better to know now:

### ⚠ "The full chain" is Stage A complete, Stage B skeletal, Stage C partial

The campaign is designed as three stages (`SIM_CAMPAIGN_PLAN.md` §2):

| stage | what it does | state |
|---|---|---|
| **A** — Geant4 | beam/background → ionization clusters in the gas | ✅ **works, this is what you can run today** |
| **B** — drift/avalanche/induction → VMM | clusters → per-pad charge and time | 🟡 **skeleton works** (2026-08-07). `python3 -m stage_b.run stageA.root` produces per-pad charge and time end to end. Four physics inputs are still placeholders — read `stage_b/README.md` §3 before quoting anything from it |
| **C** — electronics/analysis | PDO/TDO, thresholds, clustering, observables | 🟡 partial — the VMM emulator (`vmm/`) and the timing toy (`vmm/time_resolution.py`) run standalone, on synthetic input, not on Stage A output |

So today you can produce ionization/edep/conversion-budget physics — which is
what the gas decision actually rests on (`PHOTON_DISCRIMINATION_NOTES.md`) —
plus pad charge and pad multiplicity through Stage B. What you cannot yet get
is anything that needs the electronics joined on (efficiency vs threshold, a
time resolution from real tracks): Stage C still runs standalone on synthetic
input, not on Stage B output.

Stage B's placeholders matter for how you read it: gas transport is not
Magboltz (P0.10), the induced current is a delta rather than an electron
spike plus ion tail (P0.17), there are no pillars (P0.16, ~4.8 % efficiency
bias), and mesh transparency is a constant (P0.13). All four are recorded in
every Stage B output file.

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

### ✅ HTCondor submission now speaks `p2` — use the right script

**Use `scripts/submit_condor_p2.py`** for campaign runs. It knows `-m p2` and
the P2 geometry flags, the angle scan, `--skip-empty`, the §9 naming scheme,
and `--beam-spread`; it refuses to submit an unknown gas (it asks the binary
via `--list-gases`), and it writes a manifest of what was intended so you can
diff it against the `RunMeta` of what came back.

```bash
python3 scripts/submit_condor_p2.py --dry-run --scan photon
python3 scripts/submit_condor_p2.py --scan all --outdir /eos/user/<u>/<you>/p2
```

*Historical note, because it explains a class of result you may inherit:*
until 2026-08-07 `submit_condor.py` passed no `-m`, and since `p2` is the
default mode it silently ran P2 wedge geometry under vacuum-mode tags. It now
passes `-m vacuum` explicitly. Separately, `collect_results.py` had the
campaign gas names missing from a hardcoded regex and **silently dropped**
every file it could not parse; it now groups on `RunMeta` and reports
anything it cannot place. Any pre-08-07 output is suspect on both counts.

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

### Garfield++, if you need it

```bash
bash scripts/setup_garfield.sh --check     # what would be used, and why
source scripts/setup_garfield.sh           # set it up
bash scripts/setup_garfield.sh --build     # build the pin (~15 min, one time)
```

**Which one you need depends entirely on what you are doing, and for the
thing you probably want first, the CVMFS one is fine.**

*Magboltz gas tables (P0.10) — CVMFS is fine.* Magboltz is vendored inside
Garfield at 11.19 (Jan 2024) and has not moved: between the LCG Garfield and
the current pin, `magboltz.f` differs only by a Fortran continuation marker
(354 lines) and a missing comma in a FORMAT *print* statement. No cross
section, no transport change; the Penning table is byte-identical and was
re-probed to confirm it. So `source scripts/setup_lxplus.sh` and use the LCG
view's Garfield — that unblocks Stage B's largest placeholder with no build.

*Field solves and induction (P0.13, P0.17) — you need the pin.* LCG_108 is
664 commits behind and LCG_109 is 281. Every API still exists, which is why
this fails silently rather than loudly. Missing: the neBEM OpenMP race in the
SVD inversion (wrong field solves on a multi-core box, no error), the
interface-crossing check that stops electrons tunnelling through mesh wires —
which *is* the mesh-transparency observable — and
`AvalancheMicroscopic::GetIons()` for the ion tail.

*Can you source someone else's install?* There is one on lxplus at
`/afs/cern.ch/user/d/dneff/work/garfield_install/lcg109-927e5c21`, and
`P2_SHARED_INSTALL=<path> source scripts/setup_garfield.sh` will use it. But
it is in another user's AFS area, so it needs an ACL grant from them
(`fs setacl -dir <path> -acl <your-username> read`, per directory), it dies
with their account, and CERN is deprecating AFS. `--build` gives you your own
in about 15 minutes and depends on nobody. Prefer that.

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
The pad plane is polar, so there are **two** independent ways to sit on a
boundary, and this took two attempts to get right: the original default sat
35 µm from a *radial* boundary, and moving it to a ring centre put it 61 µm
from the centre of an *azimuthal* inter-pad gap — where only 25.7 % of events
deposited in pad copper at all, against 100 % at a true pad centre. The
default is now a pad centre in both coordinates, and `mm_sim` checks it
against the real gerber-derived artwork, warning separately for a radial gap,
an azimuthal gap, or an aim point off the pad field entirely.

But a single point is still a single point: **for any pad-multiplicity,
charge-sharing or positional observable, pass `--beam-spread 11.43`** (one
ring pitch) so the impact point averages over the pad cell. MX17 lost a first
result to exactly this trap — a spurious 3× left/right asymmetry that was
entirely the beam position (`MX17_Geant/design/RESPONSE_SIM_PLAN.md` §7).

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

Submit them with `scripts/submit_condor_p2.py` (P0.7, done 2026-08-07).

And through Stage B, on any of the above:

```bash
python3 -m stage_b.selftest                        # 16 checks, no Geant4
python3 -m stage_b.run <stageA.root> -o padhits.root
```

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
| P0.9 | Stage B skeleton done; Stage C not yet joined to it |
| P0.10 | no Magboltz tables; generate wet variants alongside dry |
| P0.16 | pillars: file **and** geometry code both missing |
| P0.18 | `p2_fcu_coverage = 0.983` is inflated — the gerber script mis-reads `RotRect` apertures. Do not quote it |
| P0.4b/P0.4c | segment + exit-particle exports |
| P0.6 | PAI model in the gas region |
| — | Alexandra's geometry answers (`HANDOFF.md` §3) still gate production runs |
| — | the phase-space file (P0.15) — the highest-value external ask we have |
