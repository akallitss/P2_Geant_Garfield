# Reply to the MX17 handoff "primary-ionization check, Heed vs Geant4" (P2 side, 2026-10-16)

Reproduced independently in `P2_Geant_Garfield` with our own code, not `slab.cc`:
- `mm_sim -m slab`: a bare gas box, 10×10 cm, in a G4_Galactic world.
- `--pai config|params`: PAI in DriftGas+AmpGas only, in a region with its own cuts object.
- `--gas-step`: our 100 µm user step limit in the gas, on or off.
- `--no-pixe`.

Setup: Geant4 **11.2.1** (LCG_106), option4, 100 MeV e⁻ at normal incidence, 100k events per slab point and 10k per full-stack point.
Analysis: `scripts/analysis/heed_vs_g4_slab.py`. Tables: `heed_vs_g4_slab.csv` and `p2_stack_origin_split.csv` in this directory.

## Summary

- C1, C4 and C6 are confirmed. C2 is confirmed on the bare slab.
- C3 is confirmed only for 11.3: on 11.2.1 the AddPAIModel route works.
- C5 does not apply to us.
- **The C2 implication is contested.**
  - In the P2 stack the material adds only **×1.09–1.12 on the mean and ×1.00–1.03 on the median**, not ×1.30 / ×1.15.
  - Most of our 1.35 came from something your slab did not have: **PIXE**. Our PhysicsList has `SetPixe(true)`.

## C1. Heed: confirmed

Same `heed_mip.cc`, seed 1, our pinned Garfield build on lxplus:

| case | ours | yours |
|---|---|---|
| Ar/iso 4 mm | 39.88 | 39.88 |
| Ar/CF4/iso 4 mm | 41.52 | 41.52 |
| Ne/iso 4 mm | 18.70 | 18.69 |
| Ar/iso 30 mm | 315.0 | 314.0 |
| Ar/iso, 150 GeV µ | 42.47 | 42.62 |

The 0.3 % differences in the last two rows are presumably the Garfield version.

## C2. Bare slab: confirmed, once PIXE is off

Values are the mean / median of nPrimDrift. Heed: Ar/iso 39.88 / 29.

| slab, 3.964 mm, no step limit | ours, PIXE on | ours, PIXE off | yours |
|---|---|---|---|
| Ar/iso opt4 | 44.64 / 34 | **41.16 / 32** | 41.27 / 32 |
| Ar/iso PAI | 43.11 / 32 | **39.63 / 29** | 39.78 / 29 |
| Ar/CF4/iso opt4 | 46.63 / 36 | 43.32 / 34 | 44.28 / 34 |
| Ar/CF4/iso PAI | 45.63 / 34 | 42.19 / 32 | 42.67 / 32 |
| Ne/iso opt4 | 19.80 / 14 | 19.00 / 14 | 18.24 / 13 |
| Ne/iso PAI | 19.68 / 14 | 18.95 / 14 | 18.14 / 13 |
| Ar/iso 30 mm opt4 | 335.4 / 288 | 327.2 / 280 | 324.2 / 281 |
| Ar/iso 30 mm PAI | 327.1 / 278 | 318.2 / 270 | 316.0 / 269 |
| Ar/iso µ 150 GeV opt4 | 43.43 / 34 | 43.23 / 34 | 43.51 / 34 |
| Ar/iso µ 150 GeV PAI | 42.68 / 32 | — | 42.51 / 32 |

**PIXE adds +8 % for electrons and nothing for muons.** Geant4 has e± and hadron PIXE cross sections but none for muons.
- With PIXE on, the 100 MeV e⁻ makes 0.81 extra low-energy electrons per event in the gas (40 eV – 2 keV, Ar L/K Auger), against 0.05 for µ.
- Option4 uses Penelope e⁻ ionisation below 100 keV, which already makes its own shell vacancies. Geant4 itself warns at start-up that PIXE plus Penelope double-counts de-excitation.
- PIXE also makes the result **depend on the step length**. Our 100 µm gas step limit adds another +6–7 % (Ar/iso slab PAI 43.1 → 46.1). With PIXE off the step limit does nothing: 41.16 vs 40.97.
- **Please check whether MX17 has PIXE on.** If it does, your full-stack numbers carry this bias too.

**Remaining small disagreements:**
- **Ne/iso:** our median equals Heed's (14, mean ratio 1.01). You get 13 (0.93). Our Ne W is 35.88 eV (Ne 36.4, iso 26.0, linear). Could you check the Ne W and density in `slab.cc`? Your open question 3 does not reproduce here.
- **Ar/CF4/iso opt4:** ours is 2 % below yours, probably the CF4 W (ours 34.0, Garfield 34.3).

## C2 implication: stack excess, contested

Full P2 stack, Ar/iso, 100 MeV e⁻ at 0°, 50 µm mylar geometry, **PIXE off**, no step limit:

| | full stack | bare slab | stack / slab, mean | stack / slab, median |
|---|---|---|---|---|
| opt4 | 44.74 / 32 | 41.16 / 32 | 1.087 ± 0.017 | 1.00 |
| PAI | 44.37 / 30 | 39.63 / 29 | **1.119 ± 0.019** | 1.03 |

`ClusterTree.birthVolume` is new, and gives where each ionising track was born. Per event, with PAI and PIXE off:

| born in | primaries per event |
|---|---|
| cathode mylar | 1.06 |
| mesh | 1.04 |
| front window | 0.46 |
| cathode Al | 0.32 |
| pad copper | ~0.3 |
| **total from solids** | **~3.2 of the 4.7** stack excess |

P(n ≤ 10) is the same in stack and slab (5.8 %). The material raises the mean through a tail and does not move the bulk.

**Our old 1.35 decomposes as follows** (production settings: opt4, PIXE on, 100 µm step; stack 52.95 / 37 vs Heed 39.88 / 29, i.e. 1.328 ± 0.020):

| factor | mean |
|---|---|
| PIXE (×1.15 in the stack with the step limit) | 1.148 |
| step limit without PIXE | 1.03 ± 0.02 |
| Urban vs PAI in the stack | 1.01 |
| material | 1.12 |
| PAI vs Heed | 0.99 |
| **product** | **1.33** |

So we agree that the ionisation model is a small part. We disagree that the rest is material. For us the largest piece was PIXE, and the material is about ×1.1 on the mean.

## C3. PAI routes: version-dependent

On **11.2.1** `G4EmParameters::AddPAIModel(p, "P2GasPAI", "pai")` **works**. It is event-by-event identical to the EmConfigurator route at the same seed, and different from opt4:
- 3000 events: opt4 44.49 / 35, params 41.39 / 31, config 41.26 / 31; the first eight events are identical in params and config.

The source agrees with your reading:
- In 11.2.1 `G4VEnergyLossProcess::isIonisation` defaults to `true` (header l.512), so the activator sees it.
- In 11.3.0 it defaults to `false` (l.514), so your no-op diagnosis is right for 11.3.

We use the EmConfigurator route anyway, since it works in both versions.

## C4. Region cuts: confirmed in the source, not hit here

- 11.2.1 `G4EmModelManager::Initialise` also matches regions by `G4ProductionCuts` pointer (l.488–490).
- Our PAI region `P2GasPAI` holds DriftGas and AmpGas only (the solids stay in `P2FineCut`) and has its own 1 µm cuts object.
- Across all runs, the largest step count of any volume in any event is ≤ 630, with no 1e5-step tracks.
- Time per event in the full stack with PAI is the same as opt4 within the ~20 % job-to-job noise (15–25 ms per event).

## C5. W from Heed: does not apply

`heed_mip.cc` writes `cl.electrons.size()`, and our analysis reads only that column. We never form a W from Heed's cluster energy.

## C6. W table: confirmed

Garfield `MediumGas::GetGasInfo` gives Ar 26.4, Ne 35.4, iC4H10 23.4, CF4 34.3, CO2 33.0. Fano factors: Ar 0.17, iC4H10 0.26, CO2 0.32, CF4 0.26.
- Our `GasMixtures.cc` has iC4H10 26.0, Ne 36.4 and CF4 34.0. Not changed yet; it is on our list together with the PIXE default.
- Note: our Ar/iso W prints as 26.36 eV, not 26.38.

## Not settled

- Penning / Jesse effect on W (your question 2): a shared systematic, still open.
- The material contribution against data: our test is the SPS MIP spectrum through the chain (our validation step 3).
