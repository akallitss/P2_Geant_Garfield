# Handoff — P2 wedge Geant4 simulation

**Written:** 2026-08-05
**Status:** repo seeded, design data gathered and characterized.

> **Status 2026-08-05, end of day — START HERE for a new session:**
>
> 1. **Geometry**: implemented as the default `p2` mode → `P2_MODEL.md`
>    (what's measured vs guessed). **Not yet compiled** (no local Geant4;
>    first build on lxplus). **Alexandra (CEA Saclay) is reviewing the
>    stack numbers — production runs are gated on her answers.**
> 2. **Campaign**: fully planned → **`SIM_CAMPAIGN_PLAN.md`** (the working
>    document: physics framing, staged Geant4→Garfield/Magboltz→VMM chain,
>    phases P0–P4, energy/angle/gap/gas run matrices, gas candidate table
>    with flammability + procurement constraints, bookkeeping). Two sourced
>    research notes back it: `research/TOOLCHAIN_NOTES.md` and
>    `research/GAS_NOTES.md`.
> 3. **VMM3a emulator**: implemented + unit-tested → **`../vmm/`**
>    (`vmm/README.md`: model, verified chip facts, Stage-B interface,
>    NEEDS-DATA questions for BASKET). ATLAS Athena reference sources kept
>    verbatim in `../vmm/reference/`.
> 4. **Test beam (added later on 2026-08-05)**: an SPS test with 200 GeV
>    muons was just taken (data analysis in progress, external). The
>    comparison plan is **`TESTBEAM_PLAN.md`**: collect as-run conditions,
>    agree data deliverables with the analyzers, replicate in sim
>    (`-p muon -e 200000`, normal incidence = default gun; campaign plan
>    §4.3), tune {gain, θ, ε_mesh, ENC}, feed back into the campaign's
>    Phase-2 numbers. **Follow-up handed to Alexandra** (Dylan has no
>    logbook/run data): her single concrete ask is to fill
>    **`testbeam/TB_CONDITIONS.md`** (⬜ NEEDED template) +
>    `testbeam/tb_run_manifest.csv`; everything downstream blocks on that.
> 5. **Photon premise revised (2026-08-05, late) — read this before doing
>    any gas work.** The background band is **50–100 keV + tail**, above the
>    Ar K-edge. Analytic pass: **`research/PHOTON_DISCRIMINATION_NOTES.md`**
>    (reproduce with `python3 scripts/photon_budget.py`). Findings: Ar/Ne
>    interaction probability differs by 5.4× at 50 keV falling to 2.1× at
>    150 keV; **mesh + copper pad plane convert 15–65× more photons than
>    the gas does**, a gas-independent floor; net Ar-vs-Ne difference in
>    total fake sources ≈ 4 %, against neon's 30× worse zero-cluster floor.
>    ⇒ **argon promoted to co-baseline**, and the *conversion-layer budget*
>    replaces the sensitivity curve as the gas decision variable. The
>    **energy-deposition discrimination** idea (MIP Landau vs conversion
>    deposits) is real at 3–30× MIP — but the dominant wall-photoelectron
>    population sits at ~3–5× MIP, inside the Landau tail, and is track-like
>    in time too, so it defeats both the charge and the timing handles. Plan
>    §0, §1.1, §3, §4.2, §5 step 3a, §7.4, §8, §10 all amended.
> 6. **Implemented and verified on lxplus, 2026-08-05 (late).** The code now
>    builds, runs, and passes overlap checks — **first working P2 run.**
>    Three bugs were found and fixed on the way, two of them fatal:
>    - **`P2::WedgeOutline` produced self-intersecting polygons** for the gas
>      frame and the frame opening (the `y ≤ topCut` chord was forced in even
>      when the 60° radial edge reached the outer arc first, folding the
>      profile over itself). `G4ExtrudedSolid` core-dumped in
>      triangularisation with no indication of which solid was at fault.
>      Only the *board* profile was ever validated against the gerber, which
>      is why this survived. Fixed in `src/P2Wedge.cc` + the Python mirror,
>      plus a self-intersection guard that now throws a named error instead.
>    - **ROOT thread-safety**: worker threads build TTrees concurrently and
>      raced in `TROOT::GetListOfTypes`; `-t 4`/`-t 8` core-dumped while
>      `-t 1` was fine. Fixed with `ROOT::EnableThreadSafety()` in `mm_sim.cc`.
>    - **γ production cut suppressing fluorescence** (the one predicted from
>      analysis) — fixed via `SetDeexcitationIgnoreCut(true)` + a 1 µm
>      `P2FineCut` region over gas/mesh/pad-Cu.
>    - **P0.1 gas composition — two bugs, both fixed** (`research/GAS_FIX_NOTES.md`):
>      volume fractions were passed to `AddMaterial`, which takes *mass*
>      fractions (Ne/Iso 95/5 by volume was being built as 86.8/13.2 by mass,
>      i.e. 2.6× the intended isobutane); and every component density was the
>      0 °C STP value inside a material declared at 20 °C, so all gases were
>      7 % too dense. Now derived from molar mass and volume fraction by
>      construction, and verified against Geant4's own `printMaterial`.
>
>    New output: `ClusterTree` carries full provenance (creator process,
>    origin process/volume/position/energy, ancestor ID, converting parent's
>    birth energy and volume, time); `EventTree` carries `interactionClass`
>    plus per-layer sums including the **previously unscored** `FrontGas`,
>    `DriftCathode_Gas`, `BackGas`; a new **`VolumeTree`** records edep in
>    *every* volume. `--skip-empty` added. New files: `include/TrackOrigin.hh`,
>    `include/TrackingAction.hh`, `src/TrackingAction.cc`.
>
>    **First physics result** (4×10⁵ γ at 60 keV per gas, **after** the
>    P0.1 gas fix): 3.53×10⁻³ (Ar/Iso) vs 3.11×10⁻³ (Ne/Iso) of photons make
>    gas ionization — argon → neon removes only **(12 ± 3) %** of
>    photon-induced hits, against 77 % if only the gas mattered. Conversion
>    budget **Micromesh 2.10×10⁻³, PCB_Cu_F 1.03×10⁻³, DriftGas 3.10×10⁻⁴** —
>    walls are 90 % of fakes in argon, 97 % in neon, and the mesh rate agrees
>    between gases to 1.4 %. Fluorescence lines all present at the right
>    energies (Cu-K 7.98, Fe-K 6.36, Cr-K 5.38 keV), proving the deexcitation
>    fix. **Corrections to the analysis**: wall photoelectrons deposit a
>    *median 567 eV* against a measured MIP median of 738 eV (MIP-sized, not
>    3–5× MIP), and **argon separates slightly better than neon** — wall
>    fakes sit at 0.77× the Ar MIP median but 1.05× the Ne one.
> 7. **⬜ Spatial distributions are the biggest remaining unknown**
>    (Dylan, 2026-08-05) — we do not know how signal electrons or background
>    photons are distributed over the wedge, in position or angle.
>    **`NEEDED_INPUTS.md`** is the new register: the ask is a phase-space
>    file at the wheel plane for both species (P0.15 consumes it). Matters
>    because accidental-coincidence rates go as the square/cube of the
>    *local* hit probability, so radial non-uniformity can outweigh the
>    entire gas choice.
> 8. **Next actions** (safe before Alexandra's answers): the Phase 0 task
>    list in `SIM_CAMPAIGN_PLAN.md` §3 — **P0.1 first: the gas-mixture
>    builder passes volume fractions to `G4Material::AddMaterial`, which
>    takes MASS fractions** (`src/DetectorConstruction.cc` `makeMix2/3`);
>    harmless-ish for ArIso, badly wrong for Ne mixes. Then gun-angle CLI,
>    cluster provenance, `--skip-empty`, first lxplus build.
>    **Newly added and high-value (2026-08-05, late):**
>    - **P0.12** — per-event interaction classification + edep by region.
>      This is what makes "plot event-by-event edep separated by interaction
>      type" a one-liner. Needs **a `TrackingAction` + `G4VUserTrackInformation`,
>      neither of which exists in this repo** (`ActionInitialization.cc`
>      registers only primary/run/event/stepping), because one-level
>      `parentID` cannot attribute a δ-ray back to its photoelectric
>      ancestor. Do it together with P0.4. Also: `FrontGas`,
>      `DriftCathode_Gas` and `BackGas` are currently **unscored**, so
>      "converted somewhere useless" is indistinguishable from "no
>      conversion".
>    - **P0.13** — woven/cylindrical mesh geometry, promoted out of Phase 4
>      to blocking: the mesh is the largest single fake source and the
>      effective-density slab errs in both directions at once.
>    - **A likely live bug**: `PhysicsList::SetCuts()` sets a 0.1 mm γ
>      production cut ≈ 5–10 keV in copper, so **Cu-K fluorescence at
>      8.05 keV is probably being suppressed** — that is the one
>      gas-sensitive wall channel. Fix with
>      `G4EmParameters::SetDeexcitationIgnoreCut(true)` (folded into P0.6).
> 9. The code now builds and runs on lxplus (item 6); nothing committed to
>    git yet (working tree is the deliverable so far). Benchmarks and
>    campaign cost estimates: **`BENCHMARKS.md`**.
>
> The remainder of this file is the original repo handoff — geometry/plan
> parts superseded as noted below, §3 open questions still current.

> **Superseded 2026-08-05 (later the same day):** steps 1–2 below are done —
> the P2 geometry is implemented as a new default `p2` mode
> (`DetectorConstruction::ConstructP2()`), including the gas envelope
> (frame + 3 mylar sheets, bulged windows). See **`docs/P2_MODEL.md`** for the
> as-built model and the assumptions that still need collaboration input
> (largely the open questions in §3 below). Step 3 (pad scoring in analysis)
> and step 4 (gas choice) unchanged. **Not yet compiled on lxplus.**
> The MX17 modes still work as described.

The goal is a single-detector simulation of the P2 wedge Micromegas, in the same
style as MX17_Geant's `vacuum` mode.

---

## 1. What is in this repo right now

### Copied verbatim from MX17_Geant

Source: `/home/dylan/CLionProjects/MX17_Geant` @ `8aaa0d2` (2026-05-30,
"Try again fix for sr90 liquid scint run"). Copied with `rsync`, excluding
`.git/`, `build/`, `cmake-build-debug/`, `.idea/`.

```
mm_sim.cc                 CLI entry point, mode dispatch
CMakeLists.txt            Geant4 + optional ROOT
include/  src/            the 10 G4 user classes (1929 lines total)
macros/run_default.mac
scripts/                  build.sh, setup_lxplus.sh, submit_condor*.py, analysis
ls_calibration/  sr90_calibration/
README.md                 <- still the MX17 README, superseded by this doc
```

Nothing in there knows about P2 yet. `README.md` at the repo root has been
replaced; the original MX17 text is preserved in git history of the source repo.

### New — P2 design data (`design/`, 54 MB)

Pulled from
`/media/dylan/LaCie/Extras/Physics/Post-Doc-Saclay/P2-detector-design/`.
Full inventory and provenance in [`design/README.md`](../design/README.md).

### New — gerber tooling (`scripts/gerber/`)

| File | Purpose |
|---|---|
| `gerber_outline.py` | dependency-free RS-274X reader (lines, arcs, flashes, regions, apertures) |
| `analyze_p2_geometry.py` | dumps outline / copper / drill geometry for any gerber set |
| `p2_wedge_model.py` | the idealised wedge profile + a figure comparing it to the gerber |
| `analyze_p2_readout.py` | pad map from the channel-mapping files |

### New — geometry reference (`docs/P2_GEOMETRY.md`)

**Read this before writing any geometry code.** Every P2 dimension recovered
from the fab files, with its source. Summary of the load-bearing numbers:

| | P2 wedge | MX17 (what the code builds today) |
|---|---|---|
| Active shape | 60° annular sector | 40 × 40 cm square |
| Radii | r = 95 → 650 mm (board), 120 → 590 mm (active) | — |
| Apex | on the beam axis, at gerber (0,0) | — |
| Readout | 1280 pads, 42 rings × (9…51), ~11.4 × 11.9 mm | unsegmented |
| Pillars | Ø0.5 mm, 2.000 mm pitch, 4.8 % of gap area | not modelled |
| Readout PCB | 2-layer, **0.236 mm** total (18 µm Cu / 200 µm FR4 / 18 µm Cu) | Kapton + 4×(Cu+FR4) + 5 mm Rohacell + Al |
| Drift gap | **3 mm or 4 mm** (two frame variants) — unconfirmed | 30 mm |

---

## 2. Recommended plan

MX17's `kVacuum` mode is the right template — it is already
"single detector, gun upstream, score drift + amplification gas". Work through
it in this order.

### Step 1 — wedge profile

`src/DetectorConstruction.cc:275-325` builds the vacuum-mode stack as a series
of `G4Box` slabs via the `MakeSlab`/`PlaceSlab` lambdas, all `detXY = 40 cm`
square. Every layer shares one profile, so replacing that single profile
converts the whole stack.

Use `G4ExtrudedSolid` built from `p2_wedge_model.outline()`. Dump the polygon
once to a generated header rather than parsing gerbers at runtime:

```bash
python scripts/gerber/p2_wedge_model.py --out /tmp/check.png   # prints the parameters
```

`docs/P2_GEOMETRY.md` §1 lists two cheaper approximations (`G4Tubs`, or a
boolean intersection) if the extruded solid proves awkward. A plain `G4Tubs`
with `rmin=95, rmax=650, sphi=-0.882°, dphi=61.763°` is a perfectly good first
pass and gets you running today.

Keep the apex at the origin so simulation coordinates match the gerber and the
mapping files directly — no offset bookkeeping.

### Step 2 — layer thicknesses

Rewrite the `t*` block at `src/DetectorConstruction.cc:210-228`. The PCB stack
in particular collapses from ~5.2 mm to 0.236 mm. Several MX17 layers (Rohacell,
Al foil, the 4× Cu/FR4 repetition) have no P2 counterpart and should go rather
than be given invented thicknesses.

Layers that are **not** in the gerbers and must come from the collaboration
before this step is trustworthy — see open questions below.

### Step 3 — pad-level scoring

`SteppingAction` already records `IonizationCluster` with `x, y, z` per step in
the gas (`include/EventData.hh:8-18`), and `SensitiveDetector` is deliberately a
no-op. So the cheapest correct route is:

- keep scoring exactly as-is,
- assign clusters to pads **in analysis**, using
  `scripts/gerber/analyze_p2_readout.py::load()` for the pad centres.

Only move pad assignment into C++ if per-pad output turns out to be needed at
run time. If you do, note that the pad map is not a regular grid — it is 42
rings with a different pad count each, so it needs a two-step (ring, then
azimuth) lookup, not a divide.

### Step 4 — gas and the rest

`GetGasMixture()` already offers ArCF4, ArIso, ArCO2, NeCF4, HeEth and more —
whichever P2 uses is very likely already defined. `PhysicsList`, `RunAction`
(CSV + optional ROOT), the Sr-90 spectrum sampler and the HTCondor submit
scripts should all carry over untouched.

### Step 5 — housekeeping

The executable, project and output names are still `mm_sim` /
`MicromegasSim` / `mm_output`. Rename when convenient; it is cosmetic and not
worth doing before the geometry works.

---

## 3. Open questions — need input before the geometry is trustworthy

These are **not** derivable from anything in `design/`. The gerbers describe
only the readout PCB and the bulk masks.

1. **Drift gap.** 3 mm or 4 mm? Two frame variants exist
   (`design/mechanical/P2_Frame_V{1_3mm,2_4mm}.stp`) and nothing says which is
   the baseline. This is the single most important unknown — it sets the
   primary ionization.
2. **Amplification gap.** Not in any file here. A bulk Micromegas with 0.5 mm
   pillars is typically 128 µm; MX17 uses 150 µm. Confirm.
3. **Gas mixture and pressure.** Assumed 1 atm.
4. **Drift cathode / entrance window.** Material and thickness. MX17 uses
   40 µm mylar + 0.1 µm Al + 50 µm Kapton + 9 µm Cu — almost certainly not what
   P2 has.
5. **Mesh.** MX17 models 30 µm stainless as a solid slab. For P2 the woven-mesh
   optical transparency matters more; decide whether to keep the solid-slab
   approximation.
6. **Which pad mapping is current.** `design/mapping/connector_*.txt` and
   `design/mapping/Mapping/connector_*.txt` agree on the ring structure and on
   1201 of 1280 pad positions (93.8 %), but differ on the rest and use opposite
   radial ordering. Ask which one the DAQ uses before trusting per-channel
   results. `analyze_p2_readout.py` reports this comparison on every run.
7. **Beam / source.** What is being simulated — MESA beam electrons at some
   energy, cosmics, a calibration source? This determines whether the P2
   equivalent of `kVacuum` even wants a gun on the axis or off it.
8. **Pillar treatment.** Pillars occupy 4.8 % of the amplification gap. Either
   ignore them, or fold them into an effective gas density — modelling 41 366
   individual pillars is not worth it.

---

## 4. Verification already done

- The idealised wedge model traces the Edge_Cuts gerber over the whole profile
  (see `docs/figures/p2_wedge_outline.png`) — the only unmodelled feature is a
  small chamfer on the bottom-right corner.
- The board outline derived from Edge_Cuts is independently confirmed by the
  bounding box of `P2_Mask1.gbr`, including the y ≤ 540.01 mm top cut.
- `Phi` in the mapping files is radians, confirmed against the `X`/`Y` columns.
  It is unlabelled in the file headers and sits next to a radius in mm, so it is
  an easy thing to get wrong.
- Pillar pitch measured as exactly 2.0000 mm in three independent windows of the
  mask, and consistent with the pillar count over the active area.

## 5. Not done

- No Geant4 source file has been touched.
- Nothing has been committed — the working tree is the deliverable. `git log` is
  empty; the GitHub repo `Dyn0402/p2_geant` is still empty upstream.
- The ~280 MB `.stl` frame meshes were left on the LaCie drive (STEP copied
  instead).
- `design/gerbers/K59V_Fx2Mec8/` has gerbers and drill only; the KiCad project,
  schematic, DXF exports, 12.7 MB STEP and six backup zips stayed on LaCie.
