# VMM3a emulation — Stage C of the simulation chain

**Written:** 2026-08-05. Status: **model implemented and unit-tested; not
yet connected to Stage B (which doesn't exist yet) — no physics run.**
Campaign context: `../docs/SIM_CAMPAIGN_PLAN.md` §2 (chain architecture),
§5 (ADC discrimination study). Toolchain citations:
`../docs/research/TOOLCHAIN_NOTES.md` §6.

## 1. Purpose

Emulate the VMM3a front-end response per pad so the campaign can apply a
*realistic* hit threshold ("final ADC cut") when computing electron
efficiency vs photon rejection. What matters for that question, in order:
(1) the charge threshold in physical units with a defensible noise floor,
(2) the ballistic deficit if a short peaking time is used, (3) neighbor
logic (changes recorded cluster size/charge → topology cuts), (4) ADC
saturation (interacts with the photon charge spectrum at high gas gain).
Waveform-level timing (µTPC, time-walk, pile-up) is out of scope until
Phase 4.

## 2. What is implemented here

| file | content |
|---|---|
| `vmm_shaper.py` | Python port of the ATLAS Athena closed-form VMM shaper (`reference/VMM_Shaper.cxx`): h(t), first-peak finding, threshold-crossing time, ballistic-deficit scaling. Self-test: `python3 vmm_shaper.py`. |
| `vmm_emulator.py` | Channel emulator: `VMMConfig` (operating point), `PadDigitizer` with tiers **T0** (charge-sum + threshold) and **T1** (shaper peak + NL + time window + ENC noise + ADC quantization/saturation). Smoke test: `python3 vmm_emulator.py`. |
| `reference/` | verbatim Athena sources (provenance + license in `reference/PROVENANCE.md`) — the ground truth for the port. |

Validated numerically: for a single electron of charge Q, the shaped peak
amplitude equals Q (t_peak ≥ 150 ns) or Q·t_peak/150 ns below — i.e. the
Athena normalization constant 1/0.411819 makes "amplitude" directly a
charge, so **thresholds are in electrons/fC, no mV conversion needed**
until the ADC model is switched on.

## 3. The shaper model (from Athena, constants by G. Iakovidis)

Per electron of charge Q arriving at t₀ (t_peak in ns, a = t_peak/1.5):

```
h(t) = Q · (1/0.411819) · s_bd · a³p₀|p₁|² ·
       [ K₀ e^(−p₀Δt) + 2|K₁| e^(−Re p₁ Δt) cos(−Im p₁ Δt + arg K₁) ]
K₀ = 1.584   K₁ = −0.792 − 0.115i   p₀ = 1.263/a   p₁ = (1.149 − 0.786i)/a
s_bd = min(1, t_peak/150 ns)        (Micromegas ion-tail ballistic deficit)
```

Channel logic (Athena `MM_ElectronicsResponseSimulation`):
- a channel **fires** if the shaped signal crosses its threshold inside the
  search window (acceptance window extended by `vmmDeadtime` = 200 ns below
  and `vmmUpperGrazeWindow` = 150 ns above);
- with **neighbor logic** on, channels adjacent to a fired one are read out
  too, re-evaluated with threshold ≈ 1 e⁻ (their recorded time still comes
  from their own — tiny — threshold crossing);
- recorded per digit: first-peak amplitude (PDO proxy) and
  threshold-crossing time (TDO proxy); digits outside the time window are
  dropped.

ATLAS defaults worth knowing (from `MM_DigitizationTool.h`, for orientation
— **not** P2 settings): global threshold 15 000 e⁻, or strip-noise-scaled
threshold ×7 (`useThresholdScaling`, their default); avalanche gain 6×10³
(Ar/CO₂ 93:7); NL **off** in ATLAS production sim; strip/ART dead time
200 ns; capacitive cross-talk to 1st/2nd strip neighbor 0.3/0.09 (strip
chambers — pad cross-talk for P2 unknown, see §6).

## 4. VMM3a chip facts (verified against the papers, `docs/research/TOOLCHAIN_NOTES.md` §6)

| parameter | value |
|---|---|
| channels / package | 64, BGA 21×21 mm² |
| shaper | 3rd-order DDF (1 real + complex-conjugate pole pair) |
| peaking time | 25 / 50 / 100 / 200 ns |
| gain | 0.5 / 1 / 3 / 4.5 / 6 / 9 / 12 / 16 mV/fC |
| threshold | global 10-bit DAC + per-channel 5-bit trim |
| outputs | 10-bit PDO (peak), 8-bit TDO (TAC: 60/100/350/650 ns ramps), 6-bit fast ADC |
| conversion / rate | ~250 ns; dead-timeless to ~4 MHz/ch (continuous mode) |
| neighbor logic | forces readout of neighbors of a triggered channel; crosses chip boundaries via bidirectional IOs |
| VMM3a extras | two-stage (mild/strong) ion-tail suppression |
| BASKET usage (P2_EXPERIMENT.md) | expected signal ~21 fC; range 62 fC–2 pC quoted; pad capacitance 75–170 pF |

Noise: the VMM3a paper (JPCS 1498 (2020) 012051, Fig. 5) shows measured ENC
of a few ×10³ e⁻ on a small MM prototype (~10 pF) as a function of gain,
and notes the DDF design targets input capacitance "much smaller than
200 pF". **No published ENC curve covers the P2 pad capacitance range
(75–170 pF) at our gain/peaking settings** — see §6.

## 5. How Stage B connects (interface contract)

Input to `PadDigitizer.digitize()` — one event:

```python
{pad_id: (arrival_times_ns, charges_electrons)}   # per-pad electron lists
```

produced by the Stage B sampler as: for each Geant4/Heed primary electron —
drift time (z/v_d + longitudinal diffusion) → mesh transparency survival →
Polya gain draw → pad assignment (transverse diffusion + pad map). The
avalanche charge arrives as one entry per primary electron (the ~ns
avalanche development is far below shaper time scales).

Output: `PadHit(pad, charge_electrons, time_ns, above_threshold, saturated)`
— feed directly to the ROC/efficiency analysis. `above_threshold=False`
marks NL-only reads (keep them separable: the topology cut may or may not
use them).

Recommended campaign scan knobs: `threshold_electrons` (the headline axis,
scan ~5–60 k e⁻ ≈ 1–10 fC), `peak_time_ns` ∈ {100, 200}, NL on/off,
`enc_electrons` ∈ {0, measured} — everything else frozen per run.

## 6. Open items (NEEDS-DATA — ask BASKET/Saclay; none block coding)

1. **Operating point**: which gain (mV/fC), peaking time, and threshold-DAC
   working point does the BASKET prototype run? (They quote 21 fC expected
   signal and 62 fC threshold-ish floor — clarify what "62 fC–2 pC" is:
   dynamic range at which gain?)
2. **Measured noise** (pedestal σ in DAC counts / equivalent fC) on real
   pads at 75–170 pF — replaces the `enc_electrons` placeholder. Published
   VMM data doesn't cover this capacitance.
3. **PDO calibration** (counts ↔ mV ↔ fC) and effective ADC range for
   saturation modeling (`adc_full_scale_mv` is a placeholder; ~8-bit
   effective resolution per the VMM papers).
4. **Neighbor-logic channel mapping**: NL acts on *chip channel* neighbors;
   after the connector mapping, chip-adjacent channels are not necessarily
   geometrically adjacent pads. Need the pad→(VMM chip, channel) map (ties
   into the two-revision mapping discrepancy, HANDOFF.md Q6) before NL
   emulation is meaningful at pad level. Until then: approximate NL with
   geometric neighbors from `analyze_p2_readout.py` and flag it.
5. **Streaming mode semantics**: which VMM readout mode will P2 stream
   (continuous mode presumably)? Affects dead-time modeling only at rate
   studies (Phase 4), not the per-event ADC cut.
6. **Pad-to-pad capacitive cross-talk**: ATLAS strip values (0.3/0.09)
   don't transfer to pads; measure or neglect with a stated bound.
7. **Ion-tail suppression** (VMM3a mild/strong): changes the effective
   ballistic deficit; not modeled — relevant only if t_peak < 150 ns is
   considered.

## 7. Planned next steps (in order, all cheap)

1. Wire a **pad adjacency helper** from `scripts/gerber/analyze_p2_readout.py`
   (ring/φ neighbors) → `NeighborMap` for the emulator.
2. **Threshold-scan driver**: given a set of Stage-B event files, produce
   ε_e(threshold) and P_γ(threshold) curves + ROC per (gas, gap) — this is
   campaign plan §5 step 4; the emulator API is already shaped for it.
3. Add **TDO quantization** (TAC ramp + BC clock) only if the coincidence-
   window study (plan §5 step 5) needs digitized times rather than ns.
4. Validate T1 against T0 on Stage-B-like toy events (expected: near-equal
   at t_peak = 200 ns for our ≤150 ns drift spans — the T0 tier is then
   preferred for the big scans).
5. *(added 2026-08-05)* **Ballistic-deficit toy study** — δ-arrival
   (photon point deposit) vs uniform-arrival-over-T_drift (track) events
   through T1 across t_peak ∈ {25, 50, 100, 200} ns × an ENC grid:
   quantify the PDO separation factor and the ε_e cost. No Stage B or
   Geant4 needed; first concrete deliverable of
   `../docs/research/TIMING_PSD_NOTES.md` (§2, §8.1).
