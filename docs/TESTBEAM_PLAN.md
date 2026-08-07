# SPS muon test beam ↔ simulation comparison — plan

**Written:** 2026-08-05. Status: **planning — test-beam data analysis in
progress (external), sim runs gated with the campaign** (geometry review,
`SIM_CAMPAIGN_PLAN.md` §10.1). Companion docs: `SIM_CAMPAIGN_PLAN.md`
(campaign this plan feeds into; muon validation run point is its §4.3),
`P2_MODEL.md` (as-built geometry), `vmm/README.md` (VMM3a emulator).

> **Follow-up handoff (2026-08-05):** Dylan does not have the logbook or
> run data — **Alexandra owns the follow-up**. The concrete ask is one
> document: fill in **`testbeam/TB_CONDITIONS.md`** (a fill-in template of
> the §1 checklist, every field marked ⬜ NEEDED, annotated with the sim
> setting each answer feeds) plus the run manifest
> `testbeam/tb_run_manifest.csv`, and agree the §2 data-deliverable format
> with whoever analyzes. Everything in §3–§7 is blocked on that; nothing
> else here needs action until the answers land.

---

## 0. Purpose

A test with **200 GeV muons at SPS** just finished (2026-08); the data is
being analyzed. This is the first real-world measurement of the detector the
campaign simulates, and it lands exactly on the campaign's weakest spots:
the **absolute gain scale, Polya θ, noise/threshold model, and
diffusion/pad-sharing** are currently literature values and assumptions
(campaign risk §10.2). The plan here: replicate the as-run test-beam
conditions in the full Stage A→B→C chain, compare a fixed set of
observables, tune the small set of free parameters, and propagate the tuned
model back into the campaign's headline numbers.

What the muon data anchors — and what it can't:

- 200 GeV muons are clean MIPs on the relativistic-rise plateau, same
  specific-ionization regime as the 30–160 MeV MESA electrons in mm of gas
  (both β≈1; muons just don't shower and barely scatter). So the test beam
  anchors the **electron-signal side** of the campaign's ROC almost
  directly: charge spectrum, cluster size, efficiency vs threshold.
- It says **nothing about the photon side** (50–150 keV rejection) — that
  stays simulation + lab sources (⁵⁵Fe).
- If the test ran an Ar-based gas, the anchor transfers to Ne mixtures only
  through the assumed Penning r — the prototype ⁵⁵Fe/⁹⁰Sr gain scan in Ne
  (campaign §7.2) is still needed.

## 1. Step 0 — pin down the as-run conditions (blocking everything else)

Nothing can be simulated faithfully until these are collected from the
test-beam team / logbook / DAQ configs. **The fill-in template exists:
`testbeam/TB_CONDITIONS.md`** (every item below expanded into ⬜ NEEDED
fields, each annotated with the sim setting it feeds), with the
machine-readable **run manifest** `testbeam/tb_run_manifest.csv` beside it
(one row per run: runID, HV_mesh, HV_drift, gas, threshold, angle,
position, N_triggers, good-flag, remarks) — same no-un-manifested-runs
discipline as campaign §9. **Assigned to Alexandra** (holds logbook + run
data).

| # | item | why it matters |
|---|---|---|
| T0.1 | Detector(s) under test: which wedge/prototype, serial, any geometry deltas vs `P2_MODEL.md` (drift gap, mesh, window stack, no-bulge if run near atmospheric) | Stage A geometry |
| T0.2 | Gas: mixture, purity, flow, P/T logs over the run period | Stage A gas + Magboltz table; P/T drift is a few-% gain effect |
| T0.3 | HV points: mesh + drift voltages, and whether an **HV scan** was taken (the single most valuable dataset for gain anchoring) | gain calibration §4 |
| T0.4 | Readout: VMM3a hybrids? peaking time, gain [mV/fC], threshold settings (DAC + measured fC), neighbor logic on/off, which mapping revision | Stage C emulator settings; the mapping revision decides the **readout order** (only 11/1280 channels agree between revisions) and therefore what NL reads and how channel-level comparisons are made |
| T0.4b | **Timing configuration** *(added 2026-08-06)*: TDO mode (threshold-crossing vs time-at-peak), TAC ramp (60/100/350/650 ns), BC clock, per-channel time calibration applied? | the σ_t comparison (§2.6) is undefined without these — `research/TIME_RESOLUTION_NOTES.md` §1 |
| T0.5 | Trigger + tracking: telescope type/resolution, scintillator coincidence, geometry of reference detectors, **and the resolution of whatever defines t = 0** | which efficiency/residual analyses are possible; σ²_measured = σ²_pad + σ²_reference, so σ_reference is mandatory for the timing comparison |
| T0.6 | Beam: μ⁺ or μ⁻, momentum spread, spot size/divergence, incidence angle(s), which pad regions illuminated (positions scanned?) | gun config; `--gun-x/y` aim points |
| T0.7 | Known pathologies: dead/noisy channels, HV trips, saturation | comparison masks |
| T0.8 | Data format + where it lives; who is analyzing what; agreement on the deliverable format of §2 | avoids duplicated/mismatched analysis |

## 2. Data-side deliverables (to request from / agree with the analyzers)

Per run (or per HV point), as histograms (ROOT/npz) **plus one summary CSV
row** so the sim comparison is scriptable:

1. **Pedestal/noise per channel** → ENC in fC at the real pad capacitance
   (75–170 pF) — direct input to the Stage C threshold model, replacing the
   estimated several-thousand-e⁻ ENC.
   **Ask for dedicated pedestal runs, and for the covariance — not just the
   per-channel σ** *(MX17, 2026-08-06 — `HANDOFF_MX17_RESPONSE.md` §2.4)*.
   MX17 derives, from dedicated pedestal acquisitions on DREAM, both the
   per-channel σ **and the common-mode covariance** between channels; they
   then inject correlated noise in simulation and let the analysis do common-
   mode subtraction exactly as it does on data. The approach transfers
   directly to a VMM bench. Why it matters here: a per-channel σ alone
   implies *independent* noise, and independent noise averages down across a
   cluster while common-mode does not — so a threshold model built on σ alone
   will be optimistic about the fake-hit rate at low threshold, which is
   precisely the regime the efficiency/ROC argument lives in. Requesting a
   pedestal run is cheap and can be done without beam.
2. **Cluster charge spectrum** (Landau): MPV + width, vs HV if scanned.
3. **Cluster size / pad multiplicity** distributions, **with the
   neighbor-logic state noted** — NL-on and NL-off multiplicities are not
   comparable numbers, and on this pad plane NL adds *channel* neighbors
   that are only sometimes physical neighbors (`vmm/README.md` §3.1). We
   expect ~1.1 pads/track at normal-ish incidence with NL off; a
   significantly larger measured value with NL off would mean the transverse
   spread (diffusion + track projection) is bigger than modeled, which
   propagates straight into the ROC.
4. **Efficiency vs threshold** (offline threshold scan on the data) and
   **vs HV** — needs the telescope or a layer-tag; the direct analog of the
   campaign's Phase-2 ε_e.
5. Spatial residuals vs telescope (if available) — secondary for us, but
   validates diffusion + pad pitch.
6. **Timing — promoted to a headline comparison (2026-08-06).** We now have
   a prediction to test: per-pad σ_t ≈ **10–13 ns in argon mixtures**,
   15–25 ns in neon mixtures, dominated by primary-ionization statistics
   (`research/TIME_RESOLUTION_NOTES.md`). Requested, per run:
   - the **TDO residual distribution vs the reference time**, delivered
     **both raw and time-walk-corrected** (we predict both; a σ quoted
     without saying which is not comparable), with σ_reference stated;
   - the **(PDO, TDO) two-dimensional distribution**, not just projections —
     that plot *is* the time-walk curve and is the sharpest single test of
     the Stage B+C chain (its slope tests drift velocity and shaper model
     together);
   - **σ_t vs peaking time**, if more than one t_p was run (free model check,
     see §2.9);
   - the drift-time *span* (as before) → Magboltz drift-velocity check;
   - per-channel time offsets: calibrated out or not? An uncalibrated
     channel-to-channel spread inflates the pooled σ and reads as poor
     detector resolution.
   Feeds the tightened coincidence window (`research/TIMING_PSD_NOTES.md`
   §4, campaign §5 step 5), whose floor is a few σ_t.
7. Stability: MPV vs time and vs P/T → validates (or motivates) a P/T gain
   correction term.
8. Rate of anomalously large signals (δ-rays, showers from the occasional
   radiative muon event) — tail of the charge spectrum, not just MPV.
9. *(added 2026-08-05)* If any runs were taken at **more than one peaking
   time** (T0.4): PDO spectra per t_p — the MIP MPV vs t_p is the
   ballistic-deficit curve that the pulse-shape-discrimination idea rides
   on (`research/TIMING_PSD_NOTES.md` §2); flag which t_p settings exist
   even if no scan was intended.

## 3. Sim-side — replicate the test beam

Works today in `p2` mode with minor additions; run **after** P0.1 (gas
mass-fraction bug) and the Phase-0 output upgrades, on the frozen geometry.

- **Beam**: `-p muon -e 200000` (CLI takes MeV; `mu-` vs `mu+` per T0.6 —
  add `muon+` if that's what ran). The default gun is already **normal
  incidence** (+z pencil beam); aim with `--gun-x/y` at the illuminated
  region. Add a small `--beam-spread-mm` / divergence option (extends P0.3)
  once T0.6 numbers exist — a pencil beam vs a cm-scale spot changes
  pad-sharing comparisons.
- **Geometry deltas** from T0.1 (e.g. `--bulge-front 0 --bulge-back 0` if
  the chamber ran at ~atmospheric overpressure with negligible sag).
- **Physics**: standard EM + PAI in the gas region (P0.6). Keep muon
  radiative processes on — rare hard δ/brems events are in the data too and
  populate the charge-spectrum tail (§2.8).
- **Statistics**: 10⁵ general, 10⁶ for efficiency-vs-threshold points
  (the plateau inefficiency is a ≲1 % effect; need the zero/low-cluster
  tail populated).
- **Stage B**: Magboltz table for the as-run gas at the as-run drift field;
  gain/θ/transparency extracted at the as-run mesh voltage(s).
- **Stage C**: VMM emulator T1 tier with as-run peaking time, gain,
  threshold, neighbor-logic state, measured ENC (§2.1), and the as-run pad
  mapping. If NL was on, model it in **channel** space
  (`vmm.nl_map.PadMap.channel_neighbors`, with the as-run mapping revision) —
  not as geometric neighbors, which would give the wrong occupancy and the
  wrong cluster size.
- **Timing**: run `vmm/time_resolution.py` at the as-run gas, gap, peaking
  time and threshold *before* the data arrives, so the σ_t comparison stays
  blind-ish; then repeat through the full chain once Stage B exists. State
  whether the prediction is walk-corrected (we default to quoting both).
- Also simulate the **campaign reference gases** at the same 200 GeV μ
  point (campaign §4.3) so the tuned-vs-untuned shift can be quoted per gas.

## 4. Comparison — observable ↔ parameter matrix

Compare data (§2) to sim (§3) per observable; tune the **minimal parameter
set** {gain anchor (↔ Penning r), Polya θ, ε_mesh, ENC} globally — not
per-plot fudge factors. Document tuned values ± uncertainties in
`docs/testbeam/tb_tuned_params.json` (consumed by Stage B/C configs).

| observable | dominant model parameter(s) constrained |
|---|---|
| Landau MPV vs HV | **absolute gain scale + Penning r** — campaign risk §10.2, the big one |
| Landau width / resolution | Polya θ × primary-cluster statistics (PAI/Heed fidelity — the P0.6 validation gets a data point) |
| cluster size (pads) | transverse diffusion, pad sharing (NL state must match on both sides) |
| efficiency vs threshold | end-to-end chain; direct ε_e anchor |
| efficiency plateau vs HV | mesh transparency + gain model |
| pedestal/ENC | Stage C noise input (replaces estimate) |
| charge-spectrum tail | δ-ray + amp-gap direct-ionization modeling |
| drift-time span | Magboltz v_drift at operating field |
| **σ_t (walk-corrected)** | **primary-cluster density n_p × drift velocity v_d — σ_t ≈ (1.0–1.5)/(n_p·v_d); it constrains the *same* pair as the drift-time span, from the other end, so the two together over-determine v_d and test the ionization model** |
| **(PDO, TDO) walk curve** | **shaper model + threshold calibration; the shape is a stronger test than either projection** |

Method: χ²/pulls per observable at the tuned point; quote which
observables the tune *cannot* reconcile — those are model-shape failures,
not calibration, and route to §6.

## 5. Feedback into model + campaign

1. Freeze `tb_tuned_params.json`; make Stage B/C read it as the default
   parameter source.
2. **Re-run the campaign Phase-2 headline numbers** (ε_e at fixed γ
   rejection, threshold in fC) tuned vs untuned; the shift is the quoted
   modeling systematic — this is the concrete "feedback to our model".
3. Cross-check the tuned chain against the lab data we already have
   (`sr90_calibration/`, ⁵⁵Fe when taken): one parameter set must
   reproduce both beam and source data or the model is over-tuned.
4. Update `P2_MODEL.md` / campaign risk table: strike anchored items,
   keep what remains assumption (photon side, Ne transfer).

## 6. If it disagrees — escalation paths (don't tune past these)

Persistent shape disagreements map to known model simplifications; upgrade
in this order (all are campaign Phase-4 items, pulled forward on demand):

- cluster-size / residual shapes → measured mesh transparency curve or
  field-map, woven-mesh geometry instead of the effective slab;
- charge spectrum shape at low charge → Heed re-ionization (Stage B
  high-fidelity mode) instead of PAI clusters;
- efficiency-vs-threshold shape → full VMM waveform (Athena h(t)) with
  time-walk instead of T1 charge model;
- anything vs time/rate → pile-up / baseline effects, out of current scope,
  flag explicitly.

## 7. Sequencing

| step | depends on | status |
|---|---|---|
| T0 conditions + run manifest (§1) | **Alexandra** (logbook/run data) | **handed off 2026-08-05** — template at `testbeam/TB_CONDITIONS.md`; collectable while data analysis runs |
| data deliverables agreed + produced (§2) | analyzers (via Alexandra) | external; agree format early |
| P0.1 gas fix, P0.4/P0.4b/P0.6 output+physics, first lxplus build | us | campaign Phase 0, already planned |
| beam-spread gun option (§3) | T0.6 | small; extends P0.3 |
| muon sim runs (§3, campaign §4.3) | frozen geometry + T0 | after Alexandra's review |
| comparison + tune (§4) | both streams | |
| campaign re-derivation + systematics memo (§5) | tune | the deliverable |
