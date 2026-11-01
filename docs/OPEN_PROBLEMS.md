# Open problems and to-do

State as of 2026-10-21 (last physics commit 7c24ff0). Read `HANDOVER.md` first.
This file has the detail once you start working.

## 1. The transparency edge (main open problem)

Simulated mesh transparency starts to fall **~150–200 V/cm too early** compared with data:
- SPS `drift_mesh_scan_1` (Ar/CO2/iso 93/5/2, mesh 450 V): at 1125 V/cm, data 0.85 vs sim 0.77.
- Cosmic det1 (Ar/iso 95/5, mesh 415 V): still 0.77 at 1375 V/cm.
- The SPS data also lose signal **below 500 V/cm**, where T = 1. SPS relative yield vs drift
  field (V/cm): 1000: 1.00, 750: 0.93, 625: 0.82, 500: 0.63, 375: 0.36, 250: 0.15.

What is known:
- **Contamination fit** (`scripts/analysis/sps_contamination_fit.py`,
  `docs/results/contamination/`): an O2-like attachment × T × 180 ns shaping fits the
  *whole* SPS scan once a low-field loss is allowed. Most of the "edge too early" then
  goes away, because max-normalising was biased by the low-field loss. But the fit needs
  ~2000–4600 ppm O2-equivalent, which is implausible and different per chamber.
  Not proven. Real three-body attachment with CO2/iso may be stronger than Magboltz's O2-only model.
- Drift-foil charging is excluded: the Al faces the gas, and there is no time or
  scan-direction dependence (`../P2_basket_analysis/RESULTS_gap_timing.md`).
- Not the gas tables: Magboltz agrees with ILC-TPC to within 5 %.
- A bent or non-planar PCB is a candidate, but **it must be shown in data, not assumed**.

**Gain** (`docs/results/gain/`): DREAM calibration 18 ADC/fC (200 fC range).
Penning r = 0.36 (0.33–0.39) for Ar/CO2/iso and 0.35 for Ar/iso. Open: the simulated
e-folding is 31–35 V at every gap, but the data give 22–28 V. Absolute gain cannot tell you the gap.

**Geometry runs, finished but NOT compared to data**:
- What: COMSOL + Garfield, HTCondor DAG 12548862.
- Variations:
  - gaps 120/135/165 µm
  - wire Ø16/20 µm
  - spacings 66/73 µm
  - a 415 V cosmic variant
  - 430/470 V × 4 gaps
- Output: 121 transparency results and 33 gain results, as the `TRANSPARENCY`/`GAIN` lines in
  `$P2_DATA/condor_transp/geom_<tag>_E<E>_<gas>.out`.

## 2. To-do, in priority order

1. **Geometry vs data.** Plot T(E) for every variation above on the SPS scan (450 V,
   ArCO2Iso9352) and the cosmic det1 scan (`_vm415`, ArIso). Do it with a free scale *and* with
   the low-field loss modelled as in the contamination fit, not max-normalised. Question:
   can any plausible geometry move the edge by 150–200 V/cm? Start from
   `scripts/analysis/sps_drift_scan_vs_sim.py` and `sps_contamination_fit.py`.
2. **Attachment fit to the SPS yield curve** above (data-session request c).
3. The same fit on the **cosmic det1 scan**. It has little low-field loss but still falls at
   high field, so it is the cross-check of the normalisation explanation.
4. Data-session requests at cosmic conditions (Ar/iso 95/5, 4 mm):
   - (a) the frac30/peak timing shift for an amplification-gap change: do −21 / −11 ns per e-fold come out?
   - (b) a gain-vs-gap table at 420 V, to turn std ln G = 0.269 into σ_gap.
   - (d) cosmic 125 V/cm (+120 ns, σ 82 ns) vs Magboltz v_d and diffusion.
5. **Efficiency chain, the deliverable**:
   1. Put T (wire losses only) and pillar losses into Stage B.
   2. Run Stage B → VMM on `rnd2026_heed`.
   3. Get VMM efficiency vs gain for the rnd2026 options (iso 5 → 10 %, drift 4 → 5 mm, Ne/iso).
6. Gain-slope mismatch: try a field-dependent Penning r, and check the SPS MPV near threshold.
7. Small items:
   - The Geant4 mesh default is still 45/63 µm (`include/SimConfig.hh`, `p2_mesh_open_um`),
     but the measured mesh is 50/68 µm.
   - W values: ours vs Garfield's (iso 26.0 vs 23.4 eV).
   - Report the interface-crossing bug upstream to Garfield++ with our patch.
   - **Transparency statistics.** Ar/iso 95/5 at 625 V/cm reruns at **0.934 ± 0.002**
     (20 seeds × 1000 electrons; the lxplus build and the Feynman build agree to 0.0004,
     2026-11-01). The published 0.943 (`results/mesh/`) was a single 3000-electron run, 2σ high.
     The other entries in the table are single runs too. Rerun the final table with ≥ 20 000
     electrons per point before quoting it to 0.1 %.

## 3. Validation plan (agreed 2026-10-15, one stage at a time)

1. v_d from SPS hit times. **Done:** agrees at ≥ 625 V/cm; data are slower at 250–500 V/cm,
   cause open (`docs/results/drift_velocity/`).
2. Fe55 spectrum and gain slope on the cosmic bench.
3. SPS MIP spectrum and pads per hit through the chain.
4. SPS timing through Stage B → VMM.
5. Drift scans plus a gas-flow scan.

If simulation and data disagree in shape, upgrade the model in this order rather than
tuning parameters (`archive/TESTBEAM_PLAN.md` §6):
1. the field map / mesh geometry
2. Heed re-ionization in Stage B
3. the full VMM waveform with time-walk
4. rate and pile-up effects, which are out of scope, so flag them

## 4. Inputs we don't have

- **Detector:** mesh calendering and thickness; the SPS HV-scan run list; overpressure
  (mbar); the 5 mm frame drawing.
- **Matthieu (P2Sim):** the beam current behind the MHz rates in `BkgHistograms.root`; the
  energy cut of "low E".
- **Signal electrons at the wheel:** we use 100 MeV, a pencil beam with 11.4 mm spread, and
  an angle scan. The real position/angle phase space would replace the scan
  (`archive/NEEDED_INPUTS.md` §1).
- **VMM3a noise at the real pad capacitance** (75–170 pF): not published, but measurable on the
  bench. It sets how low the threshold can go.
- **Gain curves for Ne mixtures** in a 128–150 µm bulk Micromegas: none published, and no Penning
  coefficient for Ne/iC4H10. A ⁵⁵Fe scan on a prototype would settle it before trusting the Ne/iso
  option's gain.

## 5. Warnings for future work

- **Induced signals / weighting field:** don't treat the gaps between pads as grounded
  copper. MX17 measured a 27 % error in prompt charge from exactly that
  (`archive/HANDOFF_MX17_RESPONSE.md` §2.6).
- **Neighbour logic is off in the baseline.** On this pad plane VMM neighbour logic reads
  chip-channel neighbours, which mostly aren't spatial neighbours (`vmm/README.md` §3.1).
  Always state the setting with any pad multiplicity.
- **Background photons:** the walls and mesh produce most fake hits, not the gas
  (`docs/research/PHOTON_DISCRIMINATION_NOTES.md`). So changing the gas changes the fake rate
  less than you'd expect.
