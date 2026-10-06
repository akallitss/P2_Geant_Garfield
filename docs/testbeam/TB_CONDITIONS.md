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

## Known so far (Alexandra, 2026-10-06)

| item | value | sim setting |
|---|---|---|
| Beam | **150 GeV muons**, perpendicular to the chambers | `-p muon -e 150000`, `--gun-theta 0` |
| Drift gap | **4 mm, same on every chamber** (frame V2; 3.964 mm of gas, `P2_MODEL.md`) | default |
| Pressure | "1 atm overpressure" as stated. Read as ~atmospheric absolute pressure with a few-mbar overpressure; ⬜ confirm the mbar value (it sets the window sag) | `--bulge-front/back` |
| Gases | **Ar/CO₂/iC₄H₁₀ 93/5/2** and **Ar/CF₄/iC₄H₁₀ 88/10/2** | `-g ArCO2Iso9352`, `-g ArCF4Iso` |
| Front-end settings | peaking time **200 ns** (mostly) and **100 ns**; gain **3.0** and **4.5 mV/fC** | `vmm/` emulator |
| Noise | use the **measured DREAM pedestal noise**; placeholder kept for a modelled noise if needed. VMM noise from the SNR code of the previous SPS campaign | Stage C |
| Run list | ⬜ to pick, see "Runs wanted" below | `tb_run_manifest.csv` |

### Runs wanted for the comparison

For **each gas**, at 200 ns / the gain mostly used:

1. **Mesh-HV scan at fixed drift field**: one run per point. Gives the gain
   curve (cluster-charge MPV vs HV), the anchor for Stage B's absolute gain.
2. **Drift-field scan at fixed mesh HV**, if taken: electron transparency
   and the drift-velocity/time-spread check.
3. **One high-statistics run at the working point**: cluster charge
   spectrum, cluster size (pads/cluster), efficiency, residuals and timing.
4. **The same working point at 100 ns and/or the other gain**: tests the
   front-end emulator separately from the detector physics.
5. **Pedestal/noise runs** taken alongside 1-3 (and the SNR-analysis output
   from the previous campaign) for the noise model.

Per run: run number, chamber(s) and which front-end (VMM or DREAM) each was
read with, mesh/drift HV, gas, peaking time, gain, threshold, and whether
neighbour logic was on.

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
| Neighbor logic on/off | ⬜ NEEDED | cluster-size comparison validity — NL-on and NL-off pad multiplicities are not comparable numbers |
| **Which pad-mapping revision** — `connector_*.txt` or `Mapping/connector_*.txt`? The two describe the same pad plane (to 2.5 µm) but assign channels to pads in almost entirely different orders (11/1280 agree) | ⬜ NEEDED | resolves campaign risk §10.6; decides what neighbor logic actually reads (a 12 mm neighbor vs a 126 mm one) and every per-channel comparison |
| DAQ config files archived? where? | ⬜ NEEDED | provenance |

## T0.4b Timing configuration *(added 2026-08-06 — needed for the σ_t comparison)*

We predict a per-pad time resolution of ~10–13 ns in argon mixtures and
15–25 ns in neon mixtures (`../research/TIME_RESOLUTION_NOTES.md`). None of
that can be compared with the data without these:

| item | value | consumed by |
|---|---|---|
| TDO mode: threshold-crossing time or time-at-peak? (per-chip setting) | ⬜ NEEDED | which quantity the emulator must reproduce — the two differ in walk behaviour and in model risk |
| TAC ramp used (60 / 100 / 350 / 650 ns) | ⬜ NEEDED | TDO quantization (LSB = ramp/2⁸); only the 650 ns ramp is non-negligible |
| BC clock frequency and TDO→ns conversion actually applied | ⬜ NEEDED | absolute time scale |
| Per-channel time offsets / TAC slope calibration: done? | ⬜ NEEDED | an uncalibrated channel spread inflates the pooled σ_t and looks like poor detector resolution |
| Time-walk correction applied in the analysis? functional form? | ⬜ NEEDED | we predict both raw and walk-corrected σ_t; a σ quoted without this label cannot be compared |
| Runs at more than one peaking time? | ⬜ NEEDED | σ_t vs t_p and the ballistic-deficit curve, both for free |

## T0.5 Trigger and reference tracking

| item | value | consumed by |
|---|---|---|
| Trigger: scintillators / telescope / self-triggered streaming | ⬜ NEEDED | which efficiency definition is possible |
| **What defines t = 0, and what is its time resolution?** (scintillator σ_t, telescope, or just the trigger BC — a 25 ns BC used as reference contributes 7.2 ns on its own) | ⬜ NEEDED | σ²_measured = σ²_pad + σ²_ref; without σ_ref the detector time resolution cannot be extracted at all |
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
