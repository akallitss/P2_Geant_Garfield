# Expected time resolution — prediction for the SPS comparison

**Written:** 2026-08-06 (Dylan: "we are very interested in the time
resolution of VMM and would definitely want to compare what we get at SPS
with simulation … let's do some work to determine what we expect").
**Status: first quantitative pass done — numbers below are toy-MC
predictions, not yet anchored to any measurement.**
Tool: `vmm/time_resolution.py` (standalone, no Geant4, no Stage B, ~1 min per
configuration). Companion docs: `TIMING_PSD_NOTES.md` (timing as a
*discriminant*, a different question), `../TESTBEAM_PLAN.md` §2.6 (the data
side), `../../vmm/README.md` (the emulator).

---

## 0. Summary — what we expect to measure

Per pad hit, for a through-going MIP at normal incidence, 3 mm drift,
gain 10⁴, threshold 2 fC, ENC 3 k e⁻, **timestamp = threshold crossing,
time-walk-corrected offline using the PDO of the same hit**:

| gas | 1/(n_p·v_d) | σ_t, t_p = 50 ns | 100 ns | 200 ns |
|---|---|---|---|---|
| Ar/iC₄H₁₀ 95/5 | 7.1 ns | 11.5 ns | **10.6 ns** | 11.0 ns |
| Ar/CO₂/iC₄H₁₀ 93/5/2 | 8.1 ns | 13.5 | **12.1** | 12.5 |
| Ne/iC₄H₁₀ 90/10 | 12.1 ns | 16.7 | **15.0** | 15.0 |
| Ne/CH₄ 93/7 | 16.5 ns | 17.7 | **16.2** | 16.6 |
| Ne/CO₂/iC₄H₁₀ 95/3/2 | 24.1 ns | 28.3 | **25.0** | 24.4 |

(threshold 2 fC; robust Gaussian-core width = half the central 68 % interval;
`python3 vmm/time_resolution.py -n 3000 --thresholds 1 2 5 --enc 0 3000 6000`.)

**Headline: σ_t ≈ 10–13 ns per pad in argon mixtures and 15–25 ns in neon
mixtures**, and *these are not electronics numbers* — they are ionization
statistics. The scaling that explains the whole table is

> **σ_t ≈ (1.0 – 1.5) × 1/(n_p · v_d)**, the arrival-time spread of the first
> primary clusters, where n_p is the primary-cluster density and v_d the
> drift velocity.

Everything the front-end contributes — ENC jitter, TAC quantization, Polya
gain fluctuation, longitudinal diffusion — adds up to **< 3 ns in quadrature**
at a sane operating point. If the SPS data comes back at 10–15 ns in an
argon mix, the model works; if it comes back at 30+ ns, something outside
this model (reference-time resolution, uncorrected walk, per-channel offsets,
drift-field non-uniformity) dominates, and that is the interesting result.

**This is not the same number as the "sub-2.5 ns timing" quoted for the
VMM3a in the BASKET material** (`../P2_EXPERIMENT.md` §2). That is the
front-end's own contribution; ours is the detector hit, which is dominated by
where in the gas the first clusters happened to sit. Both are correct; only
the second one belongs in a coincidence-window or accidental-rate argument.

There is **no µTPC option here.** ATLAS NSW reaches sub-10 ns per plane by
fitting the arrival-time-vs-position slope across a multi-strip cluster; our
charge cloud fires ~1.1 pads, so each hit gives exactly one timestamp. Every
timing improvement has to come from the gas, the threshold, or the peaking
time.

## 1. What is actually being measured (definition matters)

The VMM streams one timestamp per channel — either the **threshold-crossing
time** or the **time at the peak**, a per-chip configuration. The two have
different physics and different systematics, and the test-beam configuration
must be recorded (`testbeam/TB_CONDITIONS.md` T0.4):

| | threshold crossing | time at peak |
|---|---|---|
| what it samples | arrival of the *first* charge above threshold | a charge-weighted centroid of the whole pulse |
| time walk | large, ~correctable with PDO | negligible by construction |
| threshold dependence | strong (see §3) | none |
| latency | early, ~30–70 ns after the crossing | ~t_p later |
| our predicted σ (Ar/CO₂/Iso, t_p = 100 ns) | 12.1 ns | 11.1 ns |
| model risk | low | **high** — depends on the induced-current shape (ion tail), which we approximate as delta charges (§4) |

Peak mode looks marginally better in the toy, but that margin is exactly
where the model is weakest, so **the prediction we stand behind is the
threshold-crossing one**.

What the data will report is
σ²_measured = σ²_pad + σ²_reference, so the reference-time resolution
(scintillator, telescope, or trigger BC) has to come with the data or the
comparison is meaningless — a 25 ns BC clock used as the reference alone
contributes 25/√12 = 7.2 ns.

## 2. Contribution budget

`python3 vmm/time_resolution.py --breakdown` (Ar/CO₂/Iso 93/5/2, 3 mm,
t_p = 200 ns, 2 fC, gain 10⁴), effects switched on cumulatively:

| stage | σ raw | σ walk-corrected |
|---|---|---|
| ionization statistics only | 16.3 ns | 11.5 ns |
| + longitudinal diffusion | 15.5 | 11.1 |
| + Polya gain fluctuations | 17.0 | 11.5 |
| + electronics noise (ENC 3 k e⁻) | 16.7 | 11.9 |
| + TDO quantization (100 ns TAC, 8 bit) | 16.8 | 11.9 |

Term by term:

1. **Primary-ionization statistics — everything.** The distance from the
   mesh to the first cluster is exponential with mean 1/n_p, so its arrival
   time has σ = 1/(n_p·v_d) before any electronics. In a 3 mm gap that is
   7 ns (Ar/Iso) to 24 ns (Ne/CO₂/Iso). A higher threshold waits for the
   *k*-th cluster, which adds a mean delay k/(n_p v_d) with spread
   √k/(n_p v_d): the mean slides, the width grows more slowly.
2. **Time walk is the second-largest term and it is free to remove.** The
   VMM streams PDO and TDO for the same hit, so an offline
   t → t − f(PDO) correction is always available. It is worth 17 → 12 ns
   here (30 %). **Data and simulation must state whether the quoted σ is
   walk-corrected** — this is the single easiest way to compare two
   different numbers by accident.
3. **Electronics jitter** σ = ENC/|dV/dt| at the crossing: ~1 ns at
   ENC = 3 k e⁻, ~3 ns at 6 k e⁻ on a ~30 fC pulse. Only matters if the real
   ENC at 75–170 pF turns out much worse than assumed, or if the threshold
   sits close to the noise floor (below ~4σ_ENC the crossing is noise-driven
   and σ_t blows up — visible in the toy at 0.25 fC with ENC 3 k e⁻).
4. **Polya gain fluctuations** matter only through the amplitude, i.e. they
   feed the walk; after correction they are ≲ 0.5 ns.
5. **Longitudinal diffusion**: 250 µm/√cm over the full 3 mm gap is 137 µm
   of spread, i.e. 3 ns at v_d ≈ 4.6 cm/µs — sub-dominant, and it *slightly
   helps* the raw number by smoothing the leading edge.
6. **TDO quantization**: TAC ramp / 2⁸ = 0.23 ns (60 ns ramp) to 2.5 ns
   (650 ns ramp); σ = LSB/√12 ≤ 0.7 ns. Irrelevant unless the longest ramp
   is used. The coarse BC counter is not a resolution term as long as the
   TAC interpolates within it.

## 3. Knobs and their trade-offs

**Peaking time.** σ is flat-to-slightly-better between 50 and 200 ns, and
clearly worse at 25 ns. But the efficiency cost of short peaking is severe
(ballistic deficit at fixed threshold): at 2 fC, Ar/CO₂/Iso keeps 99.4 % at
t_p = 100 ns, 96.7 % at 50 ns, **77.8 % at 25 ns**; the Ne mixes fall to
39–59 %. So the timing-optimal setting is **t_p = 100 ns**, which happens to
cost nothing in efficiency — worth proposing to BASKET if timing is a
priority, and worth checking against whatever the SPS run used.

**Threshold.** Raw σ degrades steadily with threshold (Ar/CO₂/Iso at
t_p = 200 ns, ENC 3 k e⁻: 15.6 ns at 1 fC → 17.4 at 2 fC → 22.7 at 5 fC)
because a higher threshold waits for more clusters; **walk-corrected σ is
remarkably flat (12.7 / 12.5 / 12.7 ns)** because the extra delay is mostly
amplitude-correlated. Two consequences: (i) timing does not push the threshold down,
efficiency does; (ii) any comparison with data that skips the walk
correction will show a spurious threshold dependence.

**Gas.** Timing follows 1/(n_p·v_d) — the same quantity that drives the
zero-cluster inefficiency floor (n_p) and the drift span (v_d). Argon
mixtures win on timing by ~2× over the slow neon ternary, which is *another*
entry on the argon side of the gas ledger (`SIM_CAMPAIGN_PLAN.md` §7.4).
Ne/CO₂/iC₄H₁₀ 95/3/2 is the worst case in this table on both counts.

**Gap.** Not scanned yet (§6). Expect σ ∝ roughly constant in the first
approximation — 1/(n_p v_d) does not contain the gap — but a thinner gap
raises the zero-cluster floor and truncates the tail, so the Phase-3 (gap ×
gas) matrix should carry σ_t as a cell value.

## 4. Model, and where it is weak

Chain per event (`vmm/time_resolution.py`): Poisson clusters along the track
with density n_p → electrons per cluster from p(n) ∝ n⁻² tuned to
n_tot/n_p → drift t = z/v_d with Gaussian longitudinal diffusion → Polya
gain (θ = 2) → all charge on one pad → Athena closed-form VMM shaper →
interpolated threshold crossing → analytic ENC jitter σ = ENC/slope →
TAC quantization → optional non-parametric walk correction on PDO.

Known limitations, in order of how much they could move the answer:

1. **Each electron enters the shaper as a delta charge** (Athena's model).
   The real Micromegas induced current is a fast electron spike plus a
   ~150 ns ion tail, so the true leading edge is slightly slower than
   modeled ⇒ **the electronics-jitter term is a lower bound** and peak-mode
   timing is optimistic.

   **This is much cheaper to fix than "Phase 4" implies** *(MX17 review,
   2026-08-06 — `../HANDOFF_MX17_RESPONSE.md` §2.1)*. MX17 needs a
   semi-spectral quasi-static solver because its stack is **resistive**;
   P2's is not, and for a non-resistive N-layer stack the weighting
   potential is **static** — no time dependence, no FEM. Riegler's
   closed-form layered solutions (JINST 11 (2016) P11002) are already
   implemented in Garfield++ `ComponentParallelPlate` (`AddPixel` /
   `AddStrip` on an N-layer stack). Concretely, replace each electron's
   delta charge with

   > `i(t) = Q_e·δ_fast(t) + Q_ion·i_ion(t; µ_ion, gap)`,
   > weighted by the pad's Ramo potential,

   with the ion term an analytic uniform-field current. That is an
   afternoon's work against the ~150 ns ion drift, not a phase, and it is
   the *only* thing standing between us and a trustworthy time-at-peak
   prediction — which the toy says is the better estimator (§1). It also
   retires `MM_ION_FLOW_TIME_NS = 150.0` in `vmm/vmm_shaper.py`: that is
   ATLAS NSW's gap and gas encoded as a ballistic-deficit constant, and
   P2's differ. **Drop it once real induction exists rather than retuning
   it.** Re-filed as **P0.17** (was `SIM_CAMPAIGN_PLAN.md` §8).
2. **v_d and σ_L are placeholders** (the drift times quoted in the campaign
   docs: 60–75 ns for Ar mixes, ~110 ns for Ne/CO₂ over 3 mm). Since
   σ_t ∝ 1/v_d, a 20 % error in drift velocity is a 20 % error in the
   prediction. **These must be replaced by the Magboltz tables (P0.10)
   before the numbers are quoted outside the collaboration.**

   And the Magboltz value will not be the whole story either: MX17's bench
   experience is that **water contamination dominates drift velocity** in
   real gas systems — they measure 36.6 µm/ns against a far higher dry
   Magboltz prediction in Ar/iso, and infer 1–2 % H₂O at the SPS
   (`../HANDOFF_MX17_RESPONSE.md` §2.4). Since our prediction scales as
   1/v_d, a dry table could be optimistic by tens of percent against the
   very SPS data we are trying to match. P0.10 now generates **wet variants
   (0.5/1/2 % H₂O) alongside dry**; quote the bracket, and expect the bench
   number to land inside it rather than on the dry line.
3. **n_p from PDG component values, volume-weighted** — good to ~10–20 % for
   mixtures, and the prediction is directly proportional to it.
4. No ion-tail suppression stages, no baseline wander, no pile-up, no
   per-channel threshold-trim spread (a fixed offset, calibrated out in data
   but a *systematic* if uncorrected: 5-bit trim on a ~10-bit global DAC).
5. Single pad, normal incidence. Inclined tracks (10–40°, radius-dependent
   per `NEEDED_INPUTS.md` §2) lengthen the path — more clusters, so slightly
   *better* timing — and spread charge over more pads, where the earliest
   pad time is the natural estimator. Not yet studied.
6. No δ-rays or Landau tail from Geant4 — the cluster-size tail is the
   analytic p(n) ∝ n⁻². Stage A/B will supply the real one.

## 5. What we need from the SPS data (asks)

Added to `TESTBEAM_PLAN.md` §2.6 and `testbeam/TB_CONDITIONS.md` T0.4/T0.5:

1. **TDO mode**: threshold-crossing or time-at-peak? (per-chip setting)
2. **TAC ramp** (60/100/350/650 ns) and the **BC clock frequency**, plus the
   TDO→ns calibration actually applied (per-channel TAC slope/pedestal).
3. **Reference time**: what defines t = 0 (scintillator? telescope? trigger
   BC?) and **its resolution** — without this, σ_pad cannot be extracted.
4. **Per-channel time offsets**: were they calibrated out? An uncalibrated
   channel-to-channel spread inflates the pooled σ and looks like poor
   detector resolution.
5. **Time-walk correction**: applied or not, and with what functional form.
   Please deliver **both** the raw and the walk-corrected distribution — we
   predict both.
6. The **(PDO, TDO) 2-D distribution** per run, not just the projections.
   That plot *is* the walk curve, and it is the sharpest single test of the
   whole Stage B+C chain — its slope tests drift velocity and shaper model
   together.
7. Peaking time per run (already asked as T0.4) — if more than one setting
   exists, σ_t vs t_p is a free model check against §3.

## 6. To-do

1. ✅ Toy MC + first predictions (this document).
2. Re-run with **Magboltz v_d, D_L per gas** once P0.10 lands — the single
   biggest input uncertainty. *(Cheap: minutes.)* Run the **wet variants
   too** (§4 item 2) and quote σ_t as a dry→2 %-H₂O bracket, not a line.
2b. **Static Ramo induction** (§4 item 1, P0.17): replace the delta-charge
   input with electron spike + analytic ion tail through a Riegler layered
   weighting potential. Retires the borrowed `MM_ION_FLOW_TIME_NS` constant
   and is what makes the time-at-peak prediction quotable. *(Afternoon.)*
3. Add **σ_t as a cell value in the Phase-3 (gap × gas) matrix**
   (`SIM_CAMPAIGN_PLAN.md` §6) — the machinery is a drop-in.
4. Angle scan (0–40°) and the **multi-pad case**: which estimator (earliest
   pad, charge-weighted) is best when 1.1 pads becomes 2+ at large angle.
5. Feed the σ_t prediction into the **3-layer coincidence window** arithmetic
   (`SIM_CAMPAIGN_PLAN.md` §5 step 5 and `TIMING_PSD_NOTES.md` §4): the
   window cannot be narrower than a few σ_t, which sets an accidental-rate
   floor. With σ_t ≈ 12 ns (Ar) vs 25 ns (Ne/CO₂ ternary), the achievable
   window differs by 2× and the accidental rate by 4× — this may matter more
   for the physics than the gas's photon-conversion difference does.
6. Once Stage B exists, cross-check the toy against the full chain on the
   same configuration (should agree to a few %; if not, the toy's cluster
   model is the suspect).
7. After the SPS comparison: if data and model agree, publish σ_t(gas, t_p,
   threshold) as a BASKET operating-point recommendation.
