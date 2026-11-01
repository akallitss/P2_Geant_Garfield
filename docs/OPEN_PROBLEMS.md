# Open problems and to-do

State as of 2026-10-21 (last physics commit 7c24ff0). Read `HANDOVER.md` first.
This file has the detail once you start working.

## 0. The questions behind the work (Alexandra's notes, 2026-11-01)

Everything below serves one end question: **how do a larger drift gap and/or a different gas
change the efficiency and the timing resolution of the detector, and which of the two matters
more for track reconstruction and physics, in numbers?** For each question: what exists already,
and whether simulation or data can answer it.

1. **Efficiency and timing vs drift gap and gas.** *Simulation can answer this, without Stage B.*
   The timing toy `vmm/time_resolution.py` already runs the whole per-pad chain (primaries →
   drift → Polya gain → VMM shaper → threshold crossing). It outputs efficiency vs threshold
   *and* σ_t, with the gap as a parameter, but its gas inputs are August placeholders. Updating
   it with the validated inputs is the route to the answer (§2, item 5a). Stage B only adds
   pad-level effects (charge sharing, pillars, beam position, photon events), which hardly
   depend on the gap. It is the later, complete check (item 5b).
2. **Timing resolution vs drift gap, quantitatively.** *Simulation can answer this.*
   - With a wider gap, the first cluster is spread over a longer drift, so σ_t grows. The size
     depends on cluster statistics, v_d and diffusion.
   - Starting points:
     - first prediction: 10–13 ns (Ar mixes), 15–25 ns (Ne mixes), from `vmm/time_resolution.py`
       and `research/TIME_RESOLUTION_NOTES.md`
     - Magboltz v_d and diffusion are validated against SPS hit times at ≥ 625 V/cm
       (`results/drift_velocity/`)
   - Open: the data are slower than Magboltz at 250–500 V/cm, so stay in the validated field range
     or understand that first.
   - This is the same toy as question 1 (§2, item 5a).
3. **Timing or efficiency: which matters more for tracking and physics, in numbers?** *Neither the
   detector simulation nor detector data can answer this alone.* It needs the requirement from the
   tracking/physics side: the time window needed against pile-up at the P2 rates, and what an
   efficiency loss per wheel costs in reconstructed tracks. Get those numbers from the collaboration
   (Matthieu's rates are one input), then express the gap/gas options in them.
4. **Stability before discharges take over, in terms of applied voltages.** *Mostly data.*
   - Spark probability is not simulated.
   - Simulation provides the gain at each voltage (`mesh_gain.cc`, Penning r ≈ 0.36), and hence the
     charge per avalanche to compare with the Raether limit.
   - Spark rates vs HV are measured: the SPS stage `26_hv_spark_qa.py` (and the cosmic bench's
     `07_hv_spark_qa.py`) in `P2_basket_analysis` turn the HV-monitor current spikes into spark rate
     vs mesh HV per detector ("the HV-headroom curve"). These outputs are **not** in the data
     bundle (`sps_products` only has spark-*vetoed* MPVs), so ask for the `spark_qa_<det>.csv` files.
   - Result: the usable voltage window per gas. The gain gained by any option must fit in it.
5. **Drift voltage.** It hardly changes the avalanche itself, but it changes:
   - what reaches it: the mesh transparency, §1 here
   - timing, through v_d and diffusion

   So drift scans are the handle for understanding timing.
6. **Drift gap → more primaries, but (hypothesis) not a higher spark probability.** Sparks are set
   by the amplification field and the charge in a single avalanche, not by the gap length. But a
   larger gap does put more total charge per track into the amplification gap. **Check this
   against data** (spark rate vs drift setting at fixed mesh voltage: `26_hv_spark_qa` on the SPS
   drift scans) before relying on it.
7. **Gas in the drift gap.** The composition sets:
   - the number of primaries
   - v_d and diffusion, hence timing
   - the electron temperature, hence transparency (`results/mesh/`)
   - through the quencher fraction, the gain and stability (more isobutane: higher stable gain,
     cooler electrons)

   The rnd2026 options (iso 5 → 10 %, Ne/iso) change all of these at once, so report each effect
   separately.
8. **IMPORTANT: simulation answers the drift-gap question, but is much weaker for the
   amplification gap.** The gain depends exponentially on the gap and field, the Penning
   transfer r is only known from a fit to data, and the simulated gain slope does not match
   the data (e-folding 31–35 V simulated vs 22–28 V measured, §1). The chamber's real gap is also
   non-uniform (gain spread ×4). So treat amplification-gap and gain predictions as
   **relative, calibrated to data**, never as absolute numbers. Drift-side predictions
   (primaries, transparency vs field, timing at validated fields) are on much firmer ground.

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

## 2. To-do

Start with item 1. The main line after that is **item 5a**, which answers the drift-gap and gas
question. Items 2–4 and 6 are side threads on understanding the data; they can run in parallel.

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
5. **The drift-gap and gas answer, the deliverable:**
   - **5a. Update the timing toy** (`vmm/time_resolution.py`) with the validated inputs, then use
     it to compare the options. It runs today
     (`python3 vmm/time_resolution.py --gas ArIso_95_5 --gap 3.964 --tp 200`), but:
     0. **Housekeeping:**
        - It uses its own gas names (`ArIso_95_5`, …), not the `mm_sim` ones.
        - It lacks three of the four rnd2026 gases (ArIso9010, NeIso, NeIso8515).
        - Its default gap is 3 mm instead of 3.964.
        - The `T_drift` column it prints is the 3 mm placeholder, even with `--gap`.

        Align the names with `src/GasMixtures.cc` and add the missing gases.
     1. **Drift velocity and diffusion:** replace the placeholder `t_drift_3mm_ns` / `sigma_L` with
        `stage_b.gas.MagboltzTable().for_gas(gas, E_drift)`.
        `scripts/analysis/sps_drift_velocity_check.py` already does this for the DREAM electronics.
     2. **Primaries:** replace the PDG averages with Heed. Either sample the drift-gas clusters
        of `rnd2026_heed` events (cluster z and electrons per cluster, both gaps, real angles), or
        run `heed_mip` for other gaps. The low tail is what matters, and PDG averages miss it.
     3. **Transparency:** survival probability per electron from the Garfield table at that
        field and gas (`results/mesh/`).
     4. **Gain:** keep it a free mean value. Compare options **at equal gain** (a relative
        statement, safe from the amplification-gap problem in §0, item 8). Then compare them at
        the gain each can reach stably (§0, item 4).
     5. **Validate before predicting:** reproduce the SPS hit-time widths at ≥ 625 V/cm
        (`results/drift_velocity/`) and the VMM efficiency at the SPS settings. The as-run gain,
        peaking time, per-chip thresholds and DAC → fC conversions are in `vmm/README.md` §6b. The
        data-side threshold analysis and its efficiency model are in
        `P2_basket_analysis/sps_beam_analysis/vmm_threshold_tradeoff/`.
     6. **Output:** efficiency vs threshold and σ_t for drift gaps {3.964, 5, maybe 6} mm ×
        gases {ArIso, ArIso9010, NeIso, NeIso8515}, at equal gain, with the low tails shown.
   - **5b. Stage B → VMM, the complete check** (later): put T (wire losses only) and pillar
     losses into Stage B, export per-pad electron arrival lists (the VMM shaper tier needs them;
     Stage B now writes only the summed charge and first time), and run it on `rnd2026_heed`. It
     cross-checks 5a and adds charge sharing, pillars and the photon (e/γ) events.
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
