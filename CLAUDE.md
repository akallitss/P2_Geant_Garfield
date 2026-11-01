# P2 Micromegas simulation: working rules

Geant4 + Garfield++/Heed + Magboltz + COMSOL simulation of the P2 wedge
Micromegas. Its purpose is to recover the VMM threshold efficiency by raising the gain while
keeping timing ≲ 20 ns. **Read `docs/HANDOVER.md` for the current state,
open problems and to-do list.** `docs/HANDOFF.md` is from August 2026 and only of historical interest.

These rules come from real mistakes made during this project. They apply to
people and to AI assistants working in this repo.

## Geometry

- **Every geometry number must cite a design file** (`design/`,
  `../Detector_Drawings/`) or an explicit statement from the detector team, with its date.
  `docs/P2_GEOMETRY.md` is the register. If a number is not in a drawing,
  **ask** rather than guess, and label any assumption as an assumption. Don't use values
  from guesses or from colleagues' simulations (e.g. P2Sim's mesh pitch).
- Measured mesh: 18 µm wire, 50 µm opening, 68 µm pitch. Amplification gap: 150 µm. All mylar
  ~50 µm. The HV is negative on the mesh and drift, and the pads are grounded.

## Run conditions

- For a comparison with data, take the gas and HV from the **authoritative record**.
  For SPS July 2026 that is the corrected `run_config.json`: Ar/CO2/iso 93/5/2 for runs ≤ 58 and Ar/CF4/iso
  88/10/2 for runs 59–72. **There was no Ar/iso 95/5 at SPS**; that gas is the cosmic bench.
  The beam logbook and the DAQ config are wrong about the gas.

## Before presenting a number

- **Validate it.** Compare it with published measurements in the same geometry and gas where
  possible, or with a like-for-like validation run. Then test the method's
  systematics (step size, boundary handling, how losses are classified) and run at least one
  convergence check. State assumptions as assumptions, and don't over-claim on slides.
- **Compare full distributions, not means.** Never call two results "in
  agreement" from the mean or median alone. Show the spectrum shape, per-range fractions,
  the low tail (P(n < threshold)) and a KS distance. Lead with the observable
  that matters for the physics question; here that is usually the low tail, because it sets the
  VMM inefficiency.
- Split losses by cause (wire vs attachment vs pillars vs shaping).
- Don't present a hypothesis (e.g. "bent PCB") as an explanation until data
  support it.

## Simulation settings that must not regress

- Production Geant4: `submit_condor_p2.py --physics heed` (the default). This
  means `--heed --no-pixe --pai config`. Don't use `nPrimDrift` from the legacy batches
  (`rnd2026`, `rnd2026_mylar50`). Check RunMeta `heed=1`, `pixe=0`.
- Garfield transparency: `EnableNullCollisionSteps(true)`, and
  `EnableInterfaceCrossings()` with the patch in `scripts/garfield/patches/`
  applied to the pinned Garfield 927e5c21. Classify losses by electron status.
  COMSOL `meshHeight` is the wire-centre plane (164.85 µm for the 150 µm gap).
- Magboltz tables have no points between 3 and 15 kV/cm, so don't interpolate α across that range
  (`mesh_gain.cc` integrates α itself).

## Practicalities

- Setup (Feynman or laptop, no CERN account): `docs/SETUP.md`. `$P2_DATA` stands in for
  `/eos/user/a/akallits/p2_sim` in every path. Use `submit_condor_p2.py --slurm` (Feynman) or
  `--local N` (laptop) instead of HTCondor. Before trusting a new machine, run the checks in
  its §4.
- Production runs come from a clean git checkout, so that RunMeta records the commit.
- Compiled Garfield programs only; ROOT 6.38 cling cannot JIT them.
- Don't edit or pull sources that queued HTCondor jobs recompile from.
- `$P2_DATA` is shared, read-only reference data. Your outputs go under `$P2_OUT` (on lxplus, EOS),
  never in git and never into `$P2_DATA`. Results that back a claim get copied into `docs/results/<topic>/`
  with a README saying what produced them (script, commit, cluster id).
- Add a dated entry to the logbook (`docs/logbook/project/`) for each completed step.
- Slides: white background, black text, Arial, title + bullets, at most one figure.
