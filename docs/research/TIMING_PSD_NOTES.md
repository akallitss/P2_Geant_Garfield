# Time-structure discrimination (pulse-shape analysis) — blue-sky notes

**Written:** 2026-08-05 (Dylan: "extremely interesting angle — explore the
possibilities"). Status: **ideas + study plan, nothing simulated yet.**
Companion docs: `../SIM_CAMPAIGN_PLAN.md` §5 step 6 (campaign hook),
`../TESTBEAM_PLAN.md` §2 (data-side asks), `../../vmm/README.md` (emulator
that runs these studies), **`TIME_RESOLUTION_NOTES.md` (2026-08-06 — how
*well* we can time a hit, as opposed to this document's question of what the
time structure can *discriminate*; the two share a tool and a physics
input, 1/(n_p·v_d))**.

> **Update 2026-08-06 — the price of handle 1 is now measured.** The
> short-peaking scheme costs MIP efficiency at fixed threshold far more than
> this note assumed: at 2 fC, going t_p = 100 → 50 → 25 ns costs
> 99.4 % → 96.7 % → **77.8 %** in Ar/CO₂/Iso and 92.6 % → 78.1 % → **42.1 %**
> in Ne/CH₄ (`TIME_RESOLUTION_NOTES.md` §3). That is before any ENC penalty,
> which was the cost this note flagged. Handle 1 therefore needs a
> simultaneous threshold reduction to survive, and the ROC must be evaluated
> at matched ε_e, not at matched threshold.

> **⚠ Caution added later on 2026-08-05, read before investing here.** The
> §8 kill criterion "if the conversion budget says wall through-goers
> dominate the fakes, all timing handles lose their target population" is
> now **provisionally met**: the analytic budget in
> `PHOTON_DISCRIMINATION_NOTES.md` §4 puts wall conversions at ~15–65×
> gas conversions in the 50–100 keV band, and wall photoelectrons crossing
> the gap are track-like in time — the row already marked "irreducible for
> this handle" in §1 below is the *dominant* row, not a corner case.
> That does not kill the ⁵⁵Fe-vs-⁹⁰Sr bench test (§7.1), which remains a
> cheap and decisive model check and is worth doing regardless. It does
> mean handles 1 and 5 should not be developed further until Phase 1
> measures the conversion budget with real mesh geometry (plan P0.13).

---

## 1. The physical difference we could exploit

All ionization drifts to the mesh at v_d, so the *time profile* of charge
arriving at the amp gap encodes the *z profile* of the ionization:

| event class | z profile | arrival-time profile |
|---|---|---|
| electron track (signal) | ~uniform across the full gap | ~uniform over T_drift = d/v_d ≈ **60–75 ns** (Ar/Iso, 3 mm) to **~110 ns** (Ne/CO₂) |
| gas photoabsorption / Compton (fake) | point blob at one z (≲1 mm for ≤10 keV deposits) | **δ-like**: one arrival time, smeared only by longitudinal diffusion (~few ns) |
| wall photoelectron crossing the gap | track-like | track-like — **irreducible for this handle** (bounded instead by the conversion budget) |
| amp-gap conversion | no drift at all | prompt spike |

So the signal/background separation in *time* is a factor ~10–30 in
duration — larger than the typical separation in charge. The question is
how much of it survives the VMM3a, which streams only **PDO** (peak of the
shaped signal) and **TDO** (one timestamp: threshold-crossing or
time-at-peak, a per-chip configuration) per channel — no waveforms.

## 2. Handle 1 — ballistic-deficit spectroscopy (peaking-time as a knob)

The shaper integrates charge over ~t_p. For t_p = 200 ns ≫ T_drift, PDO ∝
total charge for both classes (the campaign baseline assumption). For
**t_p = 25 ns ≪ T_drift**:

- point deposit: all charge inside the integration window → PDO ≈ Q;
- track: only a ~t_p slice of the arrival profile integrates → PDO
  suppressed by roughly t_p/T_drift ≈ **×3–4 down**.

A converted photon *already* deposits more charge than a MIP; short
peaking multiplies that separation by another ~3–4× — the upper cut of the
charge window gets correspondingly more headroom, or equivalently the
photon-acceptance at fixed ε_e drops. Costs and caveats:

1. **ENC rises at short t_p and large C_det** (series noise ~ C/√t_p; our
   pads are 75–170 pF, far above anything in published VMM noise curves).
   The electron signal is simultaneously *suppressed* ×3–4, so S/N may
   collapse — this is the make-or-break number, and it is measurable on
   the bench (pedestal σ vs t_p at real pad capacitance; already asked for
   in `TESTBEAM_PLAN.md` §2.1).
2. The VMM3a **ion-tail suppression** stages (mild/strong) and the ~150 ns
   ion-tail ballistic deficit (Athena s_bd factor) modify the arithmetic;
   the T1 emulator has the machinery to compute this exactly.
3. Gain headroom: suppressed electron PDO may want more gas gain →
   interacts with discharge/saturation budgets.

**Study**: T1-emulator toy — δ-arrival vs uniform-arrival events across
t_p ∈ {25, 50, 100, 200} ns × ENC grid; deliverable is separation factor
and ε_e vs ENC. Needs **no Stage B and no Geant4** — can run today.

## 3. Handle 2 — slow gas / long gap as a *feature*

Anything that stretches T_drift enlarges the time separation: slower gas
(Ne/CO₂ 90/10 at ~110 ns over 3 mm — the mix's "weakness" flips sign
here), lower drift field, or a larger gap. This couples the (gas, gap,
t_p) axes: the Phase-3 matrix should carry the time-separation figure per
cell, not evaluate timing at the baseline only. Constraint: T_drift sets
the coincidence window, so slower gas trades single-hit discrimination
against accidental-rate window size — the §5-step-5 arithmetic arbitrates.

## 4. Handle 3 — TDO clustering under the 3-layer coincidence

For a real track the *first-arriving* charge always comes from the
ionization nearest the mesh, so per-layer TDO sits early and tightly
clustered relative to the true crossing time (spread ~ nearest-cluster
statistics + slewing). A photon fake converts at uniform z → its TDO is
uniform over the full drift window. Consequences:

- the optimal coincidence window is **narrower than T_drift** (accidentals
  fall ∝ window², correlated fakes more slowly — quote both);
- an **inter-layer time-residual cut** on matched triples adds rejection
  even after the coincidence.

**Study**: track-TDO pdf from Stage B + T1 per gas; validated against the
SPS muon TDO distributions (`TESTBEAM_PLAN.md` §2.6). **Partly done
2026-08-06**: the track-side TDO width is now predicted —
σ_t ≈ 12 ns (Ar mixes) to 25 ns (Ne/CO₂/Iso), walk-corrected
(`TIME_RESOLUTION_NOTES.md`). This is the *floor* on the coincidence
window: it cannot be tightened below a few σ_t without losing real triples,
so the achievable window differs by 2× between argon and the slow neon
ternary — 4× in accidental doubles, 8× in triples. What remains is the
photon-side TDO pdf (uniform over the drift window), which needs Stage B.

## 5. Handle 4 (speculative) — heterogeneous layer configs

Peaking time and TDO mode are per-chip settings, and the three wheels are
independent detectors: nothing stops running **one wheel at t_p = 25 ns as
a "shape layer"** while the other two run 200 ns for charge. For a matched
3-layer track, PDO_short/PDO_long is a per-track *duration estimator* —
photon-conversion triples (accidental or correlated) show a ratio near 1,
electron tracks near t_p/T_drift. Zero hardware change; pure DAQ config.
Same trick with TDO modes (threshold-crossing on one layer, time-at-peak
on another → their difference samples the rise duration). Evaluate in sim
before proposing; the efficiency cost lives entirely on the shape layer.

## 6. Weak/rejected handles (for the record)

- **Threshold-to-peak interval per channel**: VMM streams one timestamp,
  not two — not available without config games (§5).
- **Neighbor-logic timestamps**: NL channels record their own (tiny-
  threshold) crossing times → crude early-edge samples on adjacent pads;
  at ~1.1 pads/track there is rarely a neighbor with charge. **Dead as of
  2026-08-06**, for a second and independent reason: NL fires on *chip
  channel* neighbors, which on this pad plane are the physical neighbor only
  74 % of the time in the best mapping revision and 7.7 % in the other, and
  P(the sharing partner is in the NL set) is 36 % / 2 % (`vmm/nl_map.py`).
  The "early-edge sample on the adjacent pad" would frequently be a sample
  from a pad ~126 mm away. Baseline is NL off; revisit only at large angles
  *and* with the mapping revision confirmed.
- **6-bit fast ADC**: trigger-path feature, not a waveform; ignore.
- **Full waveform readout**: not an option for P2 streaming; an APV/SRS
  parallel readout on a *prototype* remains the diagnostic fallback if we
  ever need real pulse shapes to settle a dispute with the emulator.

## 7. Bench + beam program (the reality checks)

1. **⁵⁵Fe vs ⁹⁰Sr peaking-time scan** on the prototype — *the* decisive
   cheap measurement, with sources in hand. ⁵⁵Fe = point deposit
   (5.9 keV, ~230 e⁻ blob), ⁹⁰Sr/⁹⁰Y β = through-going near-MIP track.
   Prediction: ⁵⁵Fe PDO peak position ~invariant as t_p drops 200→25 ns
   (up to the ion-tail s_bd), while the ⁹⁰Sr MPV slides down by the
   ballistic-deficit factor; the ratio-vs-t_p curve *is* the money plot,
   and the pedestal width vs t_p on the same pads answers the ENC
   question. Do first in Ar/Iso; repeat in Ne mixes when the mixer exists
   (doubles as the §7.2 gain-scan program).
2. **SPS muon data**: TDO distributions (drift-model check, handle 3) and
   PDO spectra at whatever t_p settings were run — if the campaign
   happened to include more than one peaking time, that is ballistic-
   deficit data for free (`TESTBEAM_PLAN.md` §2.9).
3. If both confirm the model, the emulator's t_p×window optimization can
   be trusted enough to propose operating settings (and possibly the
   heterogeneous-layer scheme) to BASKET with numbers attached.

## 8. Ordered to-do list

1. T1 toy study (§2) — no dependencies, run now.
2. Add the t_p axis + time-separation observable to the Phase-2/3 analysis
   spec (done in plan §5 step 6).
3. TDO-residual/coincidence-window study once Stage B exists (§4).
4. Heterogeneous-layer ROC (§5) if 1–3 look promising.
5. Bench scan (§7.1) when the prototype + sources are available; feed
   measured ENC(t_p, C_pad) back into everything.

**Kill criteria**: if bench ENC at 75–170 pF and t_p = 25–50 ns exceeds
the suppressed track signal (S/N ≲ 3), handles 1/5 die and only handle 3
(TDO windowing, which works at 200 ns) survives; if the conversion budget
says wall through-goers dominate the fakes, all timing handles lose their
target population and the money stays on gas choice + coincidence.
