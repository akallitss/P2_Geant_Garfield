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
| T0.4 | Readout: VMM3a hybrids? peaking time, gain [mV/fC], threshold settings (DAC + measured fC), neighbor logic on/off, which mapping revision | Stage C emulator settings; resolves the 79/1280 mapping discrepancy |
| T0.5 | Trigger + tracking: telescope type/resolution, scintillator coincidence, geometry of reference detectors | which efficiency/residual analyses are possible |
| T0.6 | Beam: μ⁺ or μ⁻, momentum spread, spot size/divergence, incidence angle(s), which pad regions illuminated (positions scanned?) | gun config; `--gun-x/y` aim points |
| T0.7 | Known pathologies: dead/noisy channels, HV trips, saturation | comparison masks |
| T0.8 | Data format + where it lives; who is analyzing what; agreement on the deliverable format of §2 | avoids duplicated/mismatched analysis |

## 2. Data-side deliverables (to request from / agree with the analyzers)

Per run (or per HV point), as histograms (ROOT/npz) **plus one summary CSV
row** so the sim comparison is scriptable:

1. **Pedestal/noise per channel** → ENC in fC at the real pad capacitance
   (75–170 pF) — direct input to the Stage C threshold model, replacing the
   estimated several-thousand-e⁻ ENC.
2. **Cluster charge spectrum** (Landau): MPV + width, vs HV if scanned.
3. **Cluster size / pad multiplicity** distributions (with neighbor-logic
   state noted).
4. **Efficiency vs threshold** (offline threshold scan on the data) and
   **vs HV** — needs the telescope or a layer-tag; the direct analog of the
   campaign's Phase-2 ε_e.
5. Spatial residuals vs telescope (if available) — secondary for us, but
   validates diffusion + pad pitch.
6. Timing distributions (if VMM timestamps usable) → drift-time span check
   against the Magboltz drift velocity; also per-layer/per-hit TDO
   residuals vs the track reference time — input for the tightened
   coincidence-window idea (`research/TIMING_PSD_NOTES.md` §4).
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
  mapping.
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
| cluster size (pads) | transverse diffusion, pad sharing, neighbor logic |
| efficiency vs threshold | end-to-end chain; direct ε_e anchor |
| efficiency plateau vs HV | mesh transparency + gain model |
| pedestal/ENC | Stage C noise input (replaces estimate) |
| charge-spectrum tail | δ-ray + amp-gap direct-ionization modeling |
| drift-time span | Magboltz v_drift at operating field |

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
