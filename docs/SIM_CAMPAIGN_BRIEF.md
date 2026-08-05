# Simulation campaign brief — P2 wedge

**Written:** 2026-08-05. Audience: whoever (human or model) designs the
simulation campaign next. **Superseded for planning purposes (same day):
the campaign plan now exists — [`SIM_CAMPAIGN_PLAN.md`](SIM_CAMPAIGN_PLAN.md).**
This brief remains the capabilities/assumptions snapshot it was; its §5
"suggested axes" were absorbed and extended by the plan.

> **Stale as of 2026-08-05 night in two respects.** (1) The gas premise
> flipped — see [`research/PHOTON_DISCRIMINATION_NOTES.md`](research/PHOTON_DISCRIMINATION_NOTES.md);
> argon is now a co-baseline. (2) The output format below predates the
> provenance rework — current schema is
> [`OUTPUT_FORMAT.md`](OUTPUT_FORMAT.md). Measured speeds:
> [`BENCHMARKS.md`](BENCHMARKS.md). Open external asks:
> [`NEEDED_INPUTS.md`](NEEDED_INPUTS.md).

Companion docs:

| Doc | Contents |
|---|---|
| [`P2_MODEL.md`](P2_MODEL.md) | as-built Geant4 geometry, measured-vs-guessed table, questions for Saclay (+ `figures/p2_stack_questions.png` one-pager) |
| [`P2_EXPERIMENT.md`](P2_EXPERIMENT.md) | P2@MESA / BASKET physics context, kinematics, rates — all claims sourced |
| [`P2_GEOMETRY.md`](P2_GEOMETRY.md) | every dimension extracted from the fab files |
| [`MX17_README.md`](MX17_README.md) | inherited CLI / output / condor documentation |
| [`HANDOFF.md`](HANDOFF.md) | original repo handoff (partially superseded) |

---

## 1. Status snapshot

- Geant4 `p2` mode (default) builds the full single wedge: bulged front mylar
  window → 5 mm front gas → aluminised drift mylar → **3 mm drift gas
  (sensitive)** → 30 µm mesh → **150 µm amp gas (sensitive)** → PCB
  (18 µm Cu / 200 µm FR4 / 18 µm Cu) → 5 mm back gas → bulged back window,
  plus frame rings. Gerber coordinates (apex = beam axis at x=y=0).
- **Not yet compiled** — no Geant4 on the dev machine; first build must happen
  on lxplus (`source scripts/setup_lxplus.sh && bash scripts/build.sh`) and
  should include an overlap-check smoke run (`./build/mm_sim -n 100`).
- Geometry figures + a Python mirror of the geometry exist
  (`scripts/model/`, `docs/figures/`).
- Alexandra (CEA Saclay) is reviewing the stack numbers; answers to the seven
  questions in `figures/p2_stack_questions.png` may change defaults —
  **the campaign plan should treat the flagged parameters as provisional.**

## 2. Working assumptions (agreed with Dylan, 2026-08-05)

1. **Magnetic field: none at the detector.** There is a solenoid field
   (~0.6 T) around the target region, but the backward wheels sit **outside**
   it. Consequence: electrons arrive at the wedge as **straight tracks**; the
   field's only role is to set their incidence angle and position spread
   upstream. No B-field needs to be implemented in this simulation.
   (Provenance: Dylan, near-certain; consistent with fringe-field collimation
   analysis in `P2_EXPERIMENT.md` §3 — ~10° incidence implies B ≈ 45–70 mT at
   the wheels, effectively outside the bore.)
2. **Incidence angle: central value ~10° from the wedge normal, but scan
   0°–40°.** 10° is Dylan's fringe-field number (not publicly documented);
   30–40° would apply if the wheels were inside the bore. The scan brackets
   both. Angle is w.r.t. the beam axis = wedge normal (z).
3. **Electron energy: baseline ≈ 119 MeV** (elastic recoil of the 155 MeV
   beam at θ = 140°–150°; see `P2_EXPERIMENT.md` §3), **144 MeV** for the
   200 MeV parallel-running option. **An energy scan is wanted** alongside
   the angle scan — radiative losses in the 60 cm LH₂ target smear the
   spectrum downward, so scanning from a few tens of MeV up to ~160 MeV is
   the sensible bracket (exact grid: campaign designer's choice).
4. **Drift gap: 3 mm baseline** (confirmed), **scan 4 mm → 1 mm**. Mechanics
   prefers large; physics may prefer small for signal-to-noise (to be
   quantified by this campaign — that is one of its main purposes). The
   geometry rebuilds consistently for any `--drift-gap` value.
5. **Amplification gap: 150 µm** (confirmed).
6. **Gas: `ArIso` (Ar/iC₄H₁₀ 95/5) for now**; a multi-gas campaign comes
   later. The gas library already has ArCF4, ArCO2, NeCF4, ArCF4Iso, etc.,
   with W-values in `SteppingAction::kWValues`.
7. **Back side**: gas volume behind the PCB, mylar glued to the rear face of
   a back frame; only its 5 mm depth is invented. (STEP holds the front frame
   half only.)
8. **Window bulge**: ~3 mbar → 10 mm front / 5 mm back sag (Hencky estimate;
   1–10 mbar range only spans 7–15 mm, so this is not a sensitive parameter).
9. **Mylar foils**: outer windows 40 µm plain, drift 12 µm + 0.1 µm Al —
   provisional pending Alexandra (Q1 of the questions figure).

## 3. What the tool can do today

### Geometry / run knobs (CLI)

```
./build/mm_sim [-m p2] [-g ArIso] [-p electron] [-e 119] [-n N] [-o out]
               [-s seed] [-t threads]
               [--drift-gap mm] [--amp-gap um]
               [--bulge-front mm] [--bulge-back mm]
               [--gun-x mm] [--gun-y mm]
```

Defaults: p2 mode, ArIso, 155 MeV electrons, gun aimed at (307.4, 177.5) mm
(mid-active-area, r = 355 mm, φ = 30°), fired along +z from 5 cm upstream of
the front window. Every `p2_*` default lives in `include/SimConfig.hh`.

### Scoring & output

- Sensitive volumes: `DriftGas` and `AmpGas`. Every energy-depositing step in
  them becomes an ionization cluster: position (x, y, z in mm, gerber frame),
  edep (eV), nPrimary (edep/W with Poisson-ish rounding), track/parent ID,
  particle name, kinetic energy. Step limit 100 µm in both gas volumes.
- **ROOT output** (when built with ROOT): `EventTree` (per event: edepDrift,
  edepAmp, nPrimDrift, nPrimAmp, nClusDrift, nClusAmp, primInDrift,
  primInAmp) and `ClusterTree` (per cluster: eventID, trackID, parentID,
  x, y, z, edep, nPrimary, ke, volume, particle). CSV fallback with the same
  content (`*_events.csv`, `*_clusters.csv`), one file per thread.
- In p2 mode the per-layer accumulators inherited from MX17 also fire for
  `Micromesh` and `PCB_*` (harmless extras; only written for `full`-type
  modes' branches).
- **Pad assignment is done offline**: map clusters to the 1280 pads with
  `scripts/gerber/analyze_p2_readout.py::load()` (returns conn, chan, pad,
  x, y, r, phi per channel). Two mapping revisions disagree on 79/1280
  channels — which one the DAQ uses is still an open question; pad-level
  results should carry that caveat. **`Phi` in mapping files is radians.**
- HTCondor submission scripts exist (`scripts/submit_condor*.py`) but were
  written for MX17 modes — they need a once-over before campaign use (mode
  names, new flags).

## 4. Capability gaps the campaign plan must schedule

1. **First lxplus build + geometry validation.** Nothing has been compiled
   since the P2 geometry was written. Includes checking the overlap-checker
   output for the extruded/boolean solids and one visual sanity check.
2. **Gun direction control — required for the angle scan.**
   `PrimaryGeneratorAction` currently fires along +z only. Needed: polar
   angle w.r.t. the wedge normal (suggest `--gun-theta` [deg], plus
   `--gun-phi` for the azimuth of the tilt), with the gun position adjusted
   so the aim point stays on the wedge. Small, contained change.
   Consider also a beam-spot/angular-spread option for realism later.
3. **Photon-background mode.** `-p gamma -e 0.05…0.15` (MeV) works already;
   what's missing is only campaign design: the dominant real background is
   50–150 keV photons at ~10³× the electron rate (`P2_EXPERIMENT.md` §4), so
   photon-response runs (conversion probability per layer, signal size
   spectrum) belong in the plan.
4. **Analysis scripts for the new observables** (see §5) — nothing
   P2-specific exists yet beyond the pad map loader.
5. **Pending Alexandra**: mylar thicknesses/aluminisation, back-frame depth,
   frame opening/material, overpressure, mesh spec, pad-mapping revision.
   Re-run affected campaign points if answers move the defaults.
6. Deliberate approximations to keep in mind when interpreting results:
   solid-slab mesh (no optical transparency), pillars not folded into the amp
   gap (4.8% volume), terraced window domes (exact for normal incidence,
   approximate for oblique tracks), air world.

## 5. Suggested campaign axes (suggestions, not decisions)

Primary physics question: **signal vs drift gap** — how does primary
ionization statistics (and, downstream, S/N and position resolution) trend
from 4 mm down to 1 mm, across the realistic incidence range?

- **Axis 1 — drift gap:** {1, 1.5, 2, 2.5, 3, 3.5, 4} mm.
- **Axis 2 — incidence angle:** {0, 10, 20, 30, 40}° (10° = expected).
  Note path length through the drift gap grows as 1/cosθ (3 mm → 3.9 mm at
  40°), and oblique tracks spread clusters across pads.
- **Axis 3 — energy:** e.g. {30, 60, 90, 119, 144, 155} MeV electrons
  (119 = elastic baseline; low end probes the radiative tail — grid TBD).
- **Axis 4 — impact position:** inner / mid / outer radius (e.g. r = 150,
  355, 550 mm at φ = 30°, plus one φ near a wedge edge) — pad cells are
  area-equalized so effects should be small; verifies that.
- **Photon runs:** {50, 80, 100, 150} keV γ at the same positions/angles.
- Later: gas scan, overpressure/sag variants (cheap sensitivity check),
  mesh-transparency model comparison.

Suggested per-run observables: distribution of nPrimDrift (mean, Fano-ish
width, zero-cluster fraction = inefficiency proxy), cluster z-profile within
the drift gap, transverse cluster spread and pad multiplicity after mapping,
edep spectra, fraction of events with delta rays escaping into pads,
photon-conversion probability per layer (photon runs). Statistics: 10⁴–10⁵
events/point is cheap (thin detector); the zero-cluster tail at 1 mm drift
sets the floor.

## 6. Conventions cheat-sheet

- x/y = gerber frame, apex (beam axis) at (0,0); all mapping files share it.
- z = 0 at the **front window plane**, +z downstream (beam direction);
  pad plane at z = +8.19 mm (3 mm drift); cluster z tells you where in the
  drift gap ionization happened.
- Energies in output: edep in eV, ke in MeV. CLI `-e` is MeV.
- W-values per gas: `src/SteppingAction.cc::kWValues` (ArIso = 26.0 eV).
- Keep `include/P2Wedge.hh` + `include/SimConfig.hh` (C++),
  `scripts/model/p2_model.py` (Python mirror), and `docs/P2_MODEL.md` in
  sync when any geometry number changes.
