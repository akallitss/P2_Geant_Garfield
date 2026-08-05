# SPS muon test beam — as-run conditions (FILL-IN)

**Status: TEMPLATE — every ⬜ NEEDED below is unanswered.**
**Owner of follow-up: Alexandra** (has the logbook and run data; Dylan
2026-08-05: "we'll leave this for Alexandra to follow up on").
Created 2026-08-05 by Dylan+Claude from `../TESTBEAM_PLAN.md` §1 (the
T0.x checklist); that doc explains *why* each item is needed and what the
answers feed. Replace each ⬜ NEEDED with the value + source (logbook page,
DAQ config file, elog link). Partial answers are useful — fill what you
know, mark the rest ⬜ UNKNOWN vs ⬜ NOT-YET-LOOKED-UP.

When complete: (1) strike the TEMPLATE line above, (2) fill
`tb_run_manifest.csv` (template in this directory, one row per run),
(3) ping the simulation side — every answer maps to a concrete sim setting
(right column), and sim replication runs are blocked on this file.

---

## T0.1 Detector under test

| item | value | consumed by |
|---|---|---|
| Which detector? (P2 wedge / BASKET / other prototype; serial no.) | ⬜ NEEDED | whole geometry (Stage A) |
| Geometry deltas vs `docs/P2_MODEL.md` (drift gap, mesh, windows) | ⬜ NEEDED | `--drift-gap`, `--amp-gap`, … |
| Overpressure during the run (mbar) / visible window sag | ⬜ NEEDED | `--bulge-front/back` (0 if ~atmospheric) |
| Orientation: which face upstream; wedge rotation vs beam | ⬜ NEEDED | gun geometry |

## T0.2 Gas

| item | value | consumed by |
|---|---|---|
| Mixture + nominal fractions (vol %) | ⬜ NEEDED | `-g` flag + Magboltz table |
| Premix bottle or mixer? certificate / MFC settings | ⬜ NEEDED | fraction uncertainty on gain |
| Flow rate | ⬜ NEEDED | (context) |
| P/T logged? where? | ⬜ NEEDED | P/T gain correction (§2.7 of plan) |

## T0.3 High voltage — *the most valuable item*

| item | value | consumed by |
|---|---|---|
| Mesh (amplification) voltage(s) [V] | ⬜ NEEDED | Stage B gain point |
| Drift voltage(s) / field [V, V/cm] | ⬜ NEEDED | Magboltz drift point |
| **Was an HV scan taken? list the points** | ⬜ NEEDED | gain-vs-HV anchor — the single best dataset for calibrating the absolute gain scale (campaign risk §10.2) |
| Trips / unstable periods (runs to mask) | ⬜ NEEDED | comparison masks |

## T0.4 Readout / DAQ

| item | value | consumed by |
|---|---|---|
| Front-end: VMM3a hybrids? which DAQ (SRS?) | ⬜ NEEDED | Stage C emulator |
| Peaking time [ns] | ⬜ NEEDED | `vmm/` emulator setting |
| VMM gain [mV/fC] | ⬜ NEEDED | `vmm/` emulator setting |
| Threshold: DAC value and/or measured fC equivalent | ⬜ NEEDED | `vmm/` emulator setting |
| Neighbor logic on/off | ⬜ NEEDED | cluster-size comparison validity |
| **Which pad-mapping revision** (two revisions disagree on 79/1280 pads) | ⬜ NEEDED | resolves campaign risk §10.4 |
| DAQ config files archived? where? | ⬜ NEEDED | provenance |

## T0.5 Trigger and reference tracking

| item | value | consumed by |
|---|---|---|
| Trigger: scintillators / telescope / self-triggered streaming | ⬜ NEEDED | which efficiency definition is possible |
| Reference tracker: type, resolution, geometry vs DUT | ⬜ NEEDED | efficiency + residual analyses |
| Other detectors in the stack (multiple wedges? order, spacing) | ⬜ NEEDED | layer-tag efficiency option |

## T0.6 Beam

| item | value | consumed by |
|---|---|---|
| μ⁺ or μ⁻; beamline (H2/H4/…) | ⬜ NEEDED | `-p muon` vs `muon+` |
| Momentum + spread (nominal 200 GeV) | ⬜ NEEDED | `-e` (MeV) |
| Spot size / divergence at the DUT | ⬜ NEEDED | beam-spread gun option (to be added) |
| Incidence angle(s) — normal only, or tilts too? | ⬜ NEEDED | if tilts: needs campaign P0.3 gun-angle work |
| Position(s) on the wedge illuminated / scanned | ⬜ NEEDED | `--gun-x/y` aim points |

## T0.7 Known pathologies

| item | value | consumed by |
|---|---|---|
| Dead / noisy channels list | ⬜ NEEDED | comparison masks |
| Saturation observed? | ⬜ NEEDED | ADC-saturation modeling |
| Other logbook anomalies worth masking | ⬜ NEEDED | run selection |

## T0.8 Data and analysis ownership

| item | value | consumed by |
|---|---|---|
| Raw-data location + format | ⬜ NEEDED | access |
| Who is analyzing what; timeline | ⬜ NEEDED | avoid duplication |
| Agreement on deliverable format (`../TESTBEAM_PLAN.md` §2: histograms + one summary CSV row per run/HV point) | ⬜ NEEDED | scriptable comparison |

---

## Run manifest

Fill `tb_run_manifest.csv` (template beside this file): one row per run —
`runID, date, HV_mesh_V, HV_drift_V, gas, threshold, angle_deg, pos_x_mm,
pos_y_mm, n_triggers, good(y/n), remarks`. Same no-un-manifested-runs rule
as the simulation campaign (`../SIM_CAMPAIGN_PLAN.md` §9).
