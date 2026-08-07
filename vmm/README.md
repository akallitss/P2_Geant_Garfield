# VMM3a emulation — Stage C of the simulation chain

**Written:** 2026-08-05. **Updated 2026-08-06** (neighbor logic re-scoped for
a *pad* detector; time resolution promoted to a first-class deliverable).
Status: **model implemented and unit-tested; not yet connected to Stage B
(which doesn't exist yet) — no physics run.**
Campaign context: `../docs/SIM_CAMPAIGN_PLAN.md` §2 (chain architecture),
§5 (ADC discrimination study). Toolchain citations:
`../docs/research/TOOLCHAIN_NOTES.md` §6. Timing study:
`../docs/research/TIME_RESOLUTION_NOTES.md`.

## 1. Purpose

Emulate the VMM3a front-end response per pad so the campaign can apply a
*realistic* hit threshold ("final ADC cut") when computing electron
efficiency vs photon rejection, and so we can predict the **time resolution**
we should see at the SPS. What matters for that, in order:
(1) the charge threshold in physical units with a defensible noise floor,
(2) **per-pad timing** — the TDO width, its contribution budget, and the
coincidence window it implies (`time_resolution.py`);
(3) the ballistic deficit if a short peaking time is used;
(4) ADC saturation (interacts with the photon charge spectrum at high gas
gain);
(5) neighbor logic — **demoted 2026-08-06** (see §3.1: the charge cloud fires
~1.1 pads, and on this pad plane NL does not read the physical neighborhood
anyway; baseline is NL **off**).

> **Detector fact that drives most of this (Dylan, 2026-08-06):** the charge
> cloud fires **~1.1 pads on average**. Inter-pad clustering is a ~10 %
> minority effect — real, worth using optimally, but not a design driver.
> Everything below that assumes a strip-like multi-channel cluster (neighbor
> logic, µTPC, charge interpolation, cluster-size cuts) has to be re-derived
> for that regime rather than inherited from ATLAS NSW.

## 2. What is implemented here

| file | content |
|---|---|
| `vmm_shaper.py` | Python port of the ATLAS Athena closed-form VMM shaper (`reference/VMM_Shaper.cxx`): h(t), first-peak finding, threshold-crossing time, ballistic-deficit scaling. Self-test: `python3 vmm_shaper.py`. |
| `vmm_emulator.py` | Channel emulator: `VMMConfig` (operating point), `PadDigitizer` with tiers **T0** (charge-sum + threshold) and **T1** (shaper peak + NL + time window + ENC noise + ADC quantization/saturation). Smoke test: `python3 vmm_emulator.py`. |
| `nl_map.py` *(new 2026-08-06)* | Pad map with **both** adjacency relations: `channel_neighbors()` (what VMM NL actually reads: chan ±1) and `geometric_neighbors()` (physically touching pads, for offline clustering). Prints the mismatch report of §3.1: `python3 nl_map.py`. |
| `time_resolution.py` *(new 2026-08-06)* | Standalone toy MC for the per-pad TDO width — clusters → drift → Polya → shaper → threshold crossing → ENC jitter → TAC quantization, with a contribution breakdown. Needs no Geant4 and no Stage B: `python3 time_resolution.py --breakdown`. Results and interpretation: `../docs/research/TIME_RESOLUTION_NOTES.md`. |
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

### 3.1 Neighbor logic on a *pad* detector — re-scoped 2026-08-06

Neighbor logic is a **strip-detector feature**. On a strip chamber, chip
channel *n±1* is also the strip physically next door, so NL recovers the
sub-threshold shoulders of a charge cluster that genuinely spans several
strips. Two things break that here.

**(a) The cluster is one pad.** The charge cloud fires ~1.1 pads on average,
so ~90 % of hits have nothing to recover; the recoverable population is the
~10 % of tracks landing near a pad edge, and on those the partner pad carries
by construction the smaller share of the charge.

**(b) Channel neighbors are not spatial neighbors.** NL acts on chip-channel
number. The P2 pad → connector-pin assignment snakes through the pad plane,
so chan ±1 is a physical neighbor only sometimes — and *how* often depends on
which mapping revision the DAQ uses. Measured with `nl_map.py` over all
1280 pads:

| | rev A (`connector_*.txt`) | rev B (`Mapping/connector_*.txt`) |
|---|---|---|
| NL partner physically touching | **7.7 %** | **73.8 %** |
| median NL-partner distance | **126 mm** | **11.9 mm** |
| 90th pct / max distance | 297 / 457 mm | 38 / 54 mm |
| fraction of the 5.78 real neighbors NL reads | **3 %** | **25 %** |
| pads whose NL set contains no real neighbor | 88.8 % | 0 % |
| non-neighbor channels forced per fired pad | 1.83 | 0.52 |
| **P(the charge-sharing partner is in the NL set)** | **2 %** | **36 %** |

Read the last row as the useful yield: even in the favourable revision, NL
catches the sharing partner about a third of the time, because chan ±1 can
cover at most 2 of the 4 pad edges (per-pad range 8–51 %). In rev A, NL is
essentially a random-pad readout — chan ±1 sits a median 126 mm away, i.e.
NL would triple the data volume while reading pads on the other side of the
wedge. **This is not a problem, it is a fact about the mapping** — but it
means NL cannot be modeled as "geometric neighbors", and `nl_map.py`
therefore exposes the two relations separately.

Consequences, adopted as the campaign baseline:

1. **Baseline: NL off.** P2 does not currently plan to read out neighbors
   (Dylan, 2026-08-06: probably not, unless we find an advantage). The
   emulator default is now `neighbor_logic=False`; NL-on is a scan variant.
   The upside it could buy is bounded above by ~10 % of tracks × 36 %
   coverage ≈ **3–4 % of tracks gaining a partial second pad**, against a
   ~3× channel-occupancy/bandwidth cost — that is the trade to quantify in
   Phase 2, and the number NL has to beat.
2. **Clustering is an offline, geometric step, and it is independent of NL.**
   Pads that each pass threshold get merged using `geometric_neighbors()`.
   That is where the small inter-pad clustering that does exist is used, and
   it works with NL off. Never pass geometric neighbors to the NL model:
   that simulates a detector we do not have.
3. **If NL is ever enabled, it must be simulated in channel space** via
   `PadMap.channel_neighbors()`, with the revision stated — otherwise
   cluster-size and occupancy predictions are wrong in both directions
   (missing the far-away forced channels, inventing the near ones).
4. **Whether NL crosses the 64-channel chip boundary** changes only 20 of
   2540 pairs (0.8 %) — irrelevant at our multiplicity, so this NEEDS-DATA
   item can be closed as "does not matter".

ATLAS defaults worth knowing (from `MM_DigitizationTool.h`, for orientation
— **not** P2 settings): global threshold 15 000 e⁻, or strip-noise-scaled
threshold ×7 (`useThresholdScaling`, their default); avalanche gain 6×10³
(Ar/CO₂ 93:7); NL **off** in ATLAS production sim; strip/ART dead time
200 ns; capacitive cross-talk to 1st/2nd strip neighbor 0.3/0.09 (strip
chambers — pad cross-talk for P2 unknown, see §7).

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
| neighbor logic | forces readout of channel-number neighbors of a triggered channel; crosses chip boundaries via bidirectional IOs (**≠ spatial neighbors here** — §3.1) |
| TDO timestamp | threshold-crossing time **or** time-at-peak, per-chip config — matters for the timing comparison (§5) |
| VMM3a extras | two-stage (mild/strong) ion-tail suppression |
| BASKET usage (P2_EXPERIMENT.md) | expected signal ~21 fC; range 62 fC–2 pC quoted; pad capacitance 75–170 pF |

Noise: the VMM3a paper (JPCS 1498 (2020) 012051, Fig. 5) shows measured ENC
of a few ×10³ e⁻ on a small MM prototype (~10 pF) as a function of gain,
and notes the DDF design targets input capacitance "much smaller than
200 pF". **No published ENC curve covers the P2 pad capacitance range
(75–170 pF) at our gain/peaking settings** — see §7.

## 5. Time resolution *(new 2026-08-06)*

Timing is now a first-class deliverable, not a Phase-4 afterthought: we want
a number to compare against the SPS run. Full write-up, contribution budget
and test-beam asks: **`../docs/research/TIME_RESOLUTION_NOTES.md`**. The
short version, from `time_resolution.py` (3 mm, gain 10⁴, 2 fC, ENC 3 k e⁻,
threshold-crossing timestamp, walk-corrected with the hit's own PDO):

| gas | σ_t @ t_p = 100 ns |
|---|---|
| Ar/iC₄H₁₀ 95/5 | 10.6 ns |
| Ar/CO₂/iC₄H₁₀ 93/5/2 | 12.1 ns |
| Ne/iC₄H₁₀ 90/10 | 15.0 ns |
| Ne/CH₄ 93/7 | 16.2 ns |
| Ne/CO₂/iC₄H₁₀ 95/3/2 | 25.0 ns |

**It is a gas measurement, not an electronics measurement**: σ_t ≈
(1.0–1.5)/(n_p·v_d), i.e. the arrival spread of the first primary clusters.
ENC jitter, TAC quantization, Polya and diffusion together contribute < 3 ns
in quadrature. Time-walk correction using the streamed (PDO, TDO) pair is
worth ~30 % and must be stated on both sides of any comparison. Note that
**µTPC is not available to us** — at ~1.1 pads per hit there is one timestamp
per hit, so gas, threshold and peaking time are the only knobs.

Emulator-side consequence: T1 already returns `time_ns`; what is still
missing is TDO quantization and the peak-vs-threshold timestamp mode (§8
item 3), both cheap and now worth doing.

## 6. How Stage B connects (interface contract)

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
use them; with the NL-off baseline they simply do not occur).

Recommended campaign scan knobs: `threshold_electrons` (the headline axis,
scan ~5–60 k e⁻ ≈ 1–10 fC), `peak_time_ns` ∈ {50, 100, 200} (100 ns is the
timing optimum at no efficiency cost — §5), `enc_electrons` ∈ {0, measured};
NL **off** by default, on only as an explicit variant with a *channel*
neighbor map — everything else frozen per run.

Downstream of the digitizer, clustering merges hit pads using
`nl_map.PadMap.geometric_neighbors` — independent of the NL setting, and the
place where the ~10 % of hits that do share charge get used.

## 7. Open items (NEEDS-DATA — ask BASKET/Saclay; none block coding)

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
4. **Which mapping revision does the DAQ use?** *(rewritten 2026-08-06 — the
   earlier "the two revisions disagree on 79/1280 pads" was a rounding
   artifact in the comparison script, now fixed. The pad plane is identical
   to 2.5 µm.)* The revisions disagree on the **readout order**: only
   11/1280 (connector, channel) pairs land on the same pad. That decides
   whether NL reads neighbors 12 mm away or 126 mm away (§3.1), and it is
   needed for any per-channel comparison with test-beam data, dead-channel
   masks included. Still open, still the same question to ask — but it is a
   *channel-assignment* question, not a geometry question.
   Sub-question, now closed: whether NL crosses the 64-channel chip boundary
   changes 0.8 % of pairs — **does not matter**.
5. **Will NL be enabled at all?** Working assumption (Dylan, 2026-08-06): no.
   Confirm with BASKET, and if the answer is "we could", the Phase-2 scan
   should quote what it buys (bounded at ~3–4 % of tracks gaining a partial
   second pad) against ~3× the channel occupancy.
6. **Timing configuration** *(new)*: TDO mode (threshold-crossing vs
   time-at-peak), TAC ramp (60/100/350/650 ns), BC clock, per-channel time
   calibration. Needed to compare against the SPS timing distributions —
   see `../docs/research/TIME_RESOLUTION_NOTES.md` §5.
7. **Streaming mode semantics**: which VMM readout mode will P2 stream
   (continuous mode presumably)? Affects dead-time modeling only at rate
   studies (Phase 4), not the per-event ADC cut.
8. **Pad-to-pad capacitive cross-talk**: ATLAS strip values (0.3/0.09)
   don't transfer to pads; measure or neglect with a stated bound. Note this
   is a *different* coupling from NL and does follow physical adjacency.
9. **Ion-tail suppression** (VMM3a mild/strong): changes the effective
   ballistic deficit; not modeled — relevant only if t_peak < 150 ns is
   considered.

## 8. Planned next steps (in order, all cheap)

1. ✅ *(done 2026-08-06)* **Pad adjacency helper** → `nl_map.py`, with the
   channel and geometric relations kept separate.
2. **Threshold-scan driver**: given a set of Stage-B event files, produce
   ε_e(threshold) and P_γ(threshold) curves + ROC per (gas, gap) — this is
   campaign plan §5 step 4; the emulator API is already shaped for it.
3. Add **TDO quantization + timestamp mode** (TAC ramp, BC clock,
   threshold-crossing vs time-at-peak) to T1 — promoted from "only if
   needed" now that timing is a deliverable (§5); `time_resolution.py`
   already models both, so this is porting, not research.
4. Validate T1 against T0 on Stage-B-like toy events (expected: near-equal
   at t_peak = 200 ns for our ≤150 ns drift spans — the T0 tier is then
   preferred for the big scans).
5. *(added 2026-08-05)* **Ballistic-deficit toy study** — δ-arrival
   (photon point deposit) vs uniform-arrival-over-T_drift (track) events
   through T1 across t_peak ∈ {25, 50, 100, 200} ns × an ENC grid:
   quantify the PDO separation factor and the ε_e cost. No Stage B or
   Geant4 needed; first concrete deliverable of
   `../docs/research/TIMING_PSD_NOTES.md` (§2, §8.1). *Partial input already
   available:* `time_resolution.py` measures the track-side efficiency cost
   of short peaking at fixed threshold (at 2 fC: 99 % → 78 % going from
   t_p = 100 → 25 ns in Ar/CO₂/Iso, and 97 % → 42 % in Ne/CH₄), which is the
   price side of that trade.
6. *(added 2026-08-06)* **σ_t re-run with Magboltz v_d/D_L** once P0.10
   lands, and σ_t added as a Phase-3 matrix cell value.
