# P2 wedge simulation campaign — plan

**Written:** 2026-08-05. Status: **planning — execution gated on Alexandra's
geometry review** (see `P2_MODEL.md` §"Still to confirm"). Companion docs:
`SIM_CAMPAIGN_BRIEF.md` (capabilities snapshot), `P2_EXPERIMENT.md` (rates,
kinematics), `P2_MODEL.md` (as-built geometry),
`research/TOOLCHAIN_NOTES.md` (Geant4/Garfield++/Magboltz/VMM citations),
`research/GAS_NOTES.md` (neon mixtures, flammability, procurement
citations), `TESTBEAM_PLAN.md` (SPS 200 GeV muon test-beam comparison —
the data-anchoring track that feeds tuned parameters back into this
campaign), `research/TIMING_PSD_NOTES.md` (time-structure discrimination
blue-sky, 2026-08-05), **`research/PHOTON_DISCRIMINATION_NOTES.md`
(2026-08-05, late — the 50–100 keV cross-section reality check, the
energy-deposition discrimination handle, and the wall-conversion floor;
§0, §1, §3, §4.2, §5, §7 and §10 below are amended from it)**,
**`research/TIME_RESOLUTION_NOTES.md` (2026-08-06 — predicted per-pad σ_t
and the SPS comparison protocol; §2, §5 step 5, new §5 step 7 and §6 are
amended from it)**.

> **Amendment 2026-08-06 (Dylan).** Two operating facts now shape the Stage C
> plan: (1) the charge cloud fires **~1.1 pads on average**, so inter-pad
> clustering is a ~10 % effect and **VMM neighbor logic — a strip-detector
> feature whose channel neighbors are not spatial neighbors on this pad
> plane — is off in the baseline** (§2, `vmm/README.md` §3.1); (2) **VMM
> time resolution is a headline deliverable**, to be compared against the SPS
> data (§2, §5 step 7, `research/TIME_RESOLUTION_NOTES.md`).

---

## 0. Purpose and key questions

Single-wedge response study: how well can one P2 wedge **detect ~30–160 MeV
electrons while rejecting the 50–150 keV photon background** that arrives at
~10³× the electron rate? The experiment plans trigger-less VMM3a streaming
with online 3-wheel track matching; single-photon hits are rejected by (a) the
per-hit **ADC/charge cut** and (b) the **3-layer coincidence**. This campaign
quantifies (a) directly and provides the per-layer inputs for estimating (b).

**Premise revised 2026-08-05 (late).** The campaign was built around
"switch from argon- to neon-based gas": MX17 studies suggested Ne mixtures
are far less sensitive to keV–100 keV photons while costing only ~30 % of
the electron primary-ionization signal. Colleagues have since fixed the
expected background band at **50–100 keV with a higher-energy tail**, which
is far above the Ar K-edge, and the analytic pass in
`research/PHOTON_DISCRIMINATION_NOTES.md` finds that premise does not
survive there:

- The Ar/Ne interaction-probability ratio falls from 5.4× at 50 keV to
  2.5× at 100 keV and asymptotes to 1.8× (the electron-density ratio).
- **Conversions in the walls — steel mesh and copper pad plane — outnumber
  gas conversions by ~15× (Ar) to ~65× (Ne)**, and that floor is
  gas-independent. Total fake sources differ between Ar/Iso and Ne/Iso by
  under 10 %, against a 30× worse zero-cluster floor for neon.

So **argon and neon are now co-baselines**: the gas scan still runs (neon's
verified advantages are max-gain headroom and discharge probability, not
photon blindness), but photon rejection is no longer the argument for it,
and the conversion-layer budget — not the sensitivity curve — is the
decision variable.

**Second thread, added the same day: energy-deposition discrimination.**
Electrons are MIPs and leave a Landau spectrum; both Compton and
photoelectric conversions should leave more. The note quantifies it at
3–30× the MIP deposit in the band (not the 10–100× that ⁵⁵Fe intuition
suggests, because at 50–100 keV the photoelectron's range is 1–8 cm of gas
and it *crosses* the 3 mm gap rather than stopping in it). The handle is
real and sets the charge window, but its target population is the minority
one — wall photoelectrons crossing the gap are both dominant in rate and
MIP-like in charge *and* in time. See §5 step 3a and P0.12.

**Measured 2026-08-05** (§1.1, after the P0.1 gas fix): the wall-dominance
prediction held — walls are **90 % of fakes in argon, 97 % in neon**, and
argon → neon removes only **(12 ± 3) %** of photon-induced hits. But the
wall photoelectrons sit at a **median 567 eV against a measured MIP median
of 738 eV** — MIP-sized, not above. So the *upper* charge cut, which §0.3
and §5.4 were built around, barely touches the dominant fake class; the
**lower threshold** is the cut that acts on it, and it pays in ε_e
directly. The upper cut still works on the ~10 % of fakes that are
point-like gas conversions and fluorescence blobs. Both cuts must be
scanned independently and reported per class. Note also that **argon
separates slightly better than neon** (wall fakes at 0.77× the Ar MIP
median but 1.05× the Ne one).

The campaign must answer, per (gas, drift gap):

1. **Electron efficiency** — P(detectable signal above ADC threshold), as a
   function of energy, angle, position.
2. **Photon response** — P(fake hit above the same threshold) per incident
   photon vs energy, *and which layer the conversion happened in* (gas choice
   only helps against gas conversions; wall conversions set a floor).
3. **Discrimination power** — the ROC of electron efficiency vs photon
   acceptance as the charge **window** (lower + upper cut) scans, per
   configuration. A converted photon typically deposits *more* charge than
   a crossing electron — **3–30× the MIP deposit in the 50–100 keV band**
   (revised down from "10–100×", which applies to a fully contained few-keV
   photoabsorption such as ⁵⁵Fe, not to this band) — so the upper cut, the
   time structure (§5 step 6) and the coincidence carry the rejection; the
   lower threshold only removes partial-clip deposits and noise. Report the
   ROC **split by conversion class** (gas-photoelectric, gas-Compton,
   fluorescence re-absorption, wall photoelectron), because the classes have
   very different rejectability and very different rates.
4. **Operational numbers** — pad multiplicity (streaming bandwidth), drift
   time span (coincidence window), gain requirements.

## 1. Physics framing — what to expect before simulating

Primary-cluster statistics dominate thin-gap efficiency. Component values
from **PDG 2024 Table 35.5** (n_p = primary clusters/cm at MIP; compilations
vary — Phase 1 re-derives these with our own tools):

| component | n_p [/cm] | n_tot [/cm] | W [eV] | dE/dx min [keV/cm] |
|---|---|---|---|---|
| Ne | 13 | 40 | 37 | 1.45 |
| Ar | 25 | 97 | 26 | 2.53 |
| CH₄ | 28 | 54 | 30 | 1.61 |
| C₂H₆ | 48 | 112 | 26 | 2.91 |
| iC₄H₁₀ | 90 | 220 | 26 | 5.67 |
| CO₂ | 35 | 100 | 34 | 3.35 |
| CF₄ | 52–63 | 120 | 35–52 | 6.38 |

Volume-weighted n_p and the **zero-cluster (hard inefficiency floor)**
P₀ = exp(−n_p·d) at normal incidence (CF₄ taken at n_p = 52):

| mixture | n_p [/cm] | P₀ @ 1 mm | @ 2 mm | @ 3 mm | @ 4 mm |
|---|---|---|---|---|---|
| Ar/Iso 95/5 | 28.3 | 5.9 % | 0.35 % | 0.021 % | 0.001 % |
| Ar/CO₂/Iso 95/3/2 | 26.6 | 7.0 % | 0.49 % | 0.034 % | 0.002 % |
| Ar/CF₄/Iso 88/10/2 | 29.0 | 5.5 % | 0.30 % | 0.017 % | 0.001 % |
| Ne/Iso 95/5 | 16.9 | 18 % | 3.4 % | 0.64 % | 0.12 % |
| Ne/Iso 90/10 | 20.7 | 13 % | 1.6 % | 0.20 % | 0.025 % |
| Ne/Iso 80/20 | 28.4 | 5.8 % | 0.34 % | 0.020 % | 0.001 % |
| Ne/CF₄ 90/10 | 16.9 | 18 % | 3.4 % | 0.63 % | 0.12 % |
| Ne/CO₂/Iso 95/3/2 | 15.2 | 22 % | 4.8 % | 1.05 % | 0.23 % |
| Ne/Iso/CO₂ 90/5/5 | 18.0 | 17 % | 2.8 % | 0.46 % | 0.076 % |
| Ne/CH₄ 93/7 | 14.1 | 24 % | 6.0 % | 1.5 % | 0.36 % |
| Ne/C₂H₆ 90/10 | 16.5 | 19 % | 3.7 % | 0.71 % | 0.14 % |

Take-aways that shape the campaign:

- **Gap and gas interact strongly.** Ne/Iso 95/5 at 3 mm has the same
  cluster-count floor as Ar/Iso at ~1.7 mm. Small gaps + neon is the risky
  corner; that is where cluster-statistics fidelity (PAI/Heed, not just
  edep/W) matters most.
- **Quencher fraction is a powerful knob**: Ne/Iso 80/20 recovers Ar-level
  primary statistics — but is deeply flammable. The flammability boundary
  (ISO 10156:2017: ~2.4 % isobutane in binary Ne for any-air-ratio
  non-flammability; details §7) is therefore a first-class constraint, not
  a detail.
- The real efficiency is worse than P₀ once the ADC threshold is applied
  (few-cluster events × gain fluctuations fall below cut) — exactly what
  Phase 2 computes.
- At 50–150 keV, photon interactions in the **solids** (Cu pads, steel mesh,
  Al on cathode) do not merely compete with gas photoabsorption — the
  analytic pass says they **dominate it by 15–65×**, now confirmed in
  simulation (details and the
  per-layer table in `research/PHOTON_DISCRIMINATION_NOTES.md` §4). So the
  per-layer conversion breakdown (Phase 1) decides how much a low-Z gas
  actually buys, and on present numbers the answer is "very little". The
  per-hit charge spectra (Phase 2) decide what the ADC cut can do; the
  3-layer coincidence handles what survives, and we feed it accidental-rate
  inputs.

### 1.1 Photon side — the analytic targets (added 2026-08-05)

Numbers the first photon runs must reproduce, from
`scripts/photon_budget.py` (NIST μ/ρ, Klein-Nishina, Katz-Penfold ranges).
Per incident photon, 3 mm gap, pure gas at 1 atm:

| E [keV] | P_int(Ar) | P_int(Ne) | ratio | Compton ⟨T⟩ | photoelectron range in Ar |
|---|---|---|---|---|---|
| 50 | 3.49×10⁻⁴ | 6.49×10⁻⁵ | 5.38 | 4.0 keV | 12 mm |
| 60 | 2.32×10⁻⁴ | 5.44×10⁻⁵ | 4.27 | 5.6 keV | 17 mm |
| 80 | 1.38×10⁻⁴ | 4.48×10⁻⁵ | 3.07 | 9.4 keV | 28 mm |
| 100 | 1.02×10⁻⁴ | 4.03×10⁻⁵ | 2.53 | 13.8 keV | 41 mm |
| 150 | 7.11×10⁻⁵ | 3.45×10⁻⁵ | 2.06 | 27.2 keV | 80 mm |

Wall injection per incident 60 keV photon: mesh 2.0×10⁻³, pad Cu 2.0×10⁻³,
mylar layers 5×10⁻⁵ — versus 2.3×10⁻⁴ (Ar gas) and 5.4×10⁻⁵ (Ne gas).
Plus a **fluorescence chain** (Cu-K 8.05 keV from the pad plane, Fe-K
6.4 keV and Cr-K 5.4 keV from the stainless mesh, re-absorbed in the gas as
a point-like ~8–10× MIP blob) at ~1.8×10⁻⁴ in Ar but ~1.7×10⁻⁵ in Ne — the
one genuinely gas-sensitive wall channel, and the one the charge window
kills best.

**✅ Confirmed by simulation, 2026-08-05** (4×10⁵ γ at 60 keV in each of
Ar/Iso and Ne/Iso, 3 mm, **after** the P0.1 gas fix;
`research/PHOTON_DISCRIMINATION_NOTES.md` §4.3a, `research/GAS_FIX_NOTES.md`).
Measured photon-induced hit rate **3.53×10⁻³ (Ar) vs 3.11×10⁻³ (Ne)** — a
**(12 ± 3) % reduction**, against 77 % if only the gas mattered.
Conversion-layer budget **Micromesh 2.10×10⁻³, PCB_Cu_F 1.03×10⁻³, DriftGas
3.10×10⁻⁴** (Ar): walls are **90 % of fakes in argon, 97 % in neon**, and
the mesh rate agrees between the two gases to 1.4 % — the internal check
that the floor is gas-independent. All fluorescence lines appear at the
right energies (Cu-K 7.98/8.01, Fe-K 6.36, Cr-K 5.38 keV).

**Two predictions were wrong and both matter.** (i) Wall photoelectrons
deposit a *median 567 eV* in argon against a measured MIP median of 738 eV
— they are **MIP-sized or smaller**, not the 3–5× MIP §0 assumed. So the
**upper** charge cut barely touches the dominant fake class; only the lower
threshold does, at direct cost to ε_e. (ii) **Argon separates marginally
better than neon**: wall fakes sit at 0.77× the Ar MIP median but 1.05× the
Ne MIP median, i.e. in neon the dominant background is indistinguishable
from signal by charge alone. The cleanly separable classes (gas
photoelectric ~8× MIP, fluorescence ~9×, gas Compton ~6×, all point-like)
are only ~10 % of argon fakes and ~3 % of neon's.

Consequences for reading Phase-1 results: Compton electrons below ~80 keV
photon energy are **contained** in 3 mm (practical range 0.14 mm at 60 keV
in Ar) and deposit their full energy; photoelectrons in this band are
**not** contained and deposit only 3–11 keV of their 50–100 keV. Rayleigh
scattering deposits nothing and must not be counted as an interaction.
⬜ **NEEDED: the actual incident photon spectrum** — soft (≈50 keV) vs hard
(≥100 keV) moves the Ar/Ne ratio between 5× and 2×.

## 2. Simulation chain — architecture

Three factorized stages, joined by files (each stage re-runnable alone):

```
Stage A  Geant4 (this repo)         geometry, transport, conversions
         → ClusterTree: ionization clusters (x,y,z, edep, nPrim, provenance)
Stage B  Response (new, Python + Garfield++/Magboltz gas tables)
         drift + diffusion + mesh transparency + Polya gain + pad summing
         → pad hits: (pad, charge [fC], time)
Stage C  Electronics + analysis (new, Python)
         threshold/ADC model (VMM-informed), efficiency/ROC/multiplicity
```

**Geant4 ↔ Garfield++ coupling decision** *(research pass done 2026-08-05;
citations below)*. Garfield++ does hook into Geant4: the official interface
(Pfeiffer, De Keukeleere, Schindler, Veenhof et al., *Interfacing Geant4,
Garfield++ and Degrad for the simulation of gaseous detectors*, NIM A 935
(2019) 121, arXiv:1806.05880; `Examples/Geant4GarfieldInterface`) attaches a
Garfield "fast-simulation model" to a G4Region over the gas, killing Geant4
tracks there and handing them to Heed/Magboltz classes. It is the right tool
for single-configuration validation runs, but it makes one monolithic binary
whose slowest stage (avalanching) throttles everything and re-does Geant4
transport at every parameter point. For a campaign, the **staged pipeline is
the established pattern** — it is literally how ATLAS NSW Micromegas
digitization and the CAST/IAXO Micromegas background chains work (Geant4
deposits → Garfield/Magboltz-*parameterized* drift, diffusion, Polya gain,
electronics):

- **Geant4 standalone (Stage A)** handles all materials, photon conversions,
  deltas — everything needing the full geometry. Use Livermore/Penelope-class
  EM with atomic deexcitation ON (`/process/em/fluo`, `auger`, `pixe`) so
  Cu-K 8.0 keV / Fe-K 6.4 keV fluorescence and Auger cascades are produced.
- **In-gas ionization statistics**: Geant4's default condensed-history EM is
  *not* adequate for cluster statistics in 1–4 mm of gas (right mean, wrong
  event-by-event distribution — the thin-absorber problem). Two remedies,
  used in tiers: (i) **PAI model** in a G4Region over the gas volumes —
  benchmarked against Heed in the Pfeiffer paper, cheap (~2× EM cost);
  (ii) for the threshold/efficiency studies, **re-ionize gas crossings with
  Heed** (`TrackHeed::NewTrack()`) in Stage B from exported track segments
  (entry point, direction, momentum, PID) — Heed is the MPGD standard for
  cluster counts and cluster-size distributions, and it is what the Pfeiffer
  paper recommends for relativistic particles in gas. So Stage A exports
  *both* clusters (edep/W, PAI-informed) and gas-crossing track segments;
  Stage B chooses fidelity per study.
- **Magboltz gas tables, once per gas** (`MediumMagboltz::GenerateGasTable`
  → `WriteGasFile`): drift velocity, D_L/D_T, Townsend, attachment on a
  log E-grid spanning drift and amplification fields. Accuracy: ncoll≈10 is
  high precision; cost ~minutes per E-point for quenched mixtures at 1 atm
  (order an hour per 20-point table per core — benchmark, no official
  numbers published). Cache `.gas` files; never regenerate per job.
  **Penning transfer must be enabled** — without it simulated MM/GEM gains
  are low by large factors. Published r values (Şahin/Kowalski/Veenhof
  series, JINST 5 (2010) P05002; JINST 16 (2021) P03026; Garfield++
  built-ins): Ar/iC₄H₁₀ (10 %) r = 0.40, Ar/CO₂ 93:7 r = 0.405,
  **Ne/CO₂ 90:10 r = 0.577**, Ar/CF₄ r ≈ 0 (CF₄ IP above Ar* levels).
  **Gaps: no published r for Ne/iC₄H₁₀, Ne/CF₄, or any ternary** — the
  argument-free `EnablePenningTransfer()` works for binaries only. For our
  Ne mixes: apply the binary noble-gas r manually, bracket ±, and treat
  absolute gain as calibration-anchored (below).
- **Stage B is a lightweight sampler**, not per-electron microscopic
  avalanches: drift each electron with Gaussian diffusion σ_{L,T}·√z from
  the tables, mesh transparency ε_mesh(E_drift/E_amp) (flat ≈100 % for field
  ratios ≲0.01, steep drop above), Polya-fluctuated gain (θ ≈ 1.5–2.5 for
  Ar-mix Micromegas, Zerguerras et al., NIM A 772 (2015) 76), pad summing.
  **Electrons born inside the amp gap** (pad-Cu/mesh photoelectrons, direct
  amp-gas conversions) must get z-dependent partial gain from their birth
  depth — birth at the pad surface ⇒ ~no gain, mesh side ⇒ full gain; this
  population is a plausibly dominant wall-fake channel and a flat full-gain
  treatment misestimates it.
  Cost: µs per electron → full campaign trivial.
  **AvalancheMicroscopic is used once per (gas, HV) point** on ~10³–10⁴
  single electrons to *extract* ḡ, θ, transparency (cost scales linearly
  with avalanche size — seconds per avalanche at gain 10⁴, not viable as
  workhorse; cf. Garfield++ CUDA paper arXiv:2509.15377), then frozen into
  the sampler.
- **Absolute gain calibration**: integrated Penning-corrected Townsend,
  anchored to measured curves (BASKET prototype, published bulk-MM data).
  Report results vs *gain* and vs threshold in primary-electron equivalents
  so the gain-scale uncertainty stays factorized.
- **Degrad** (full shell-resolved photoabsorption + Auger + shake-off in
  gas; Fano-exact to <3 % below 20 keV) is needed only if the exact
  electron-count distribution of few-keV gas photoabsorption becomes
  decision-relevant (escape-peak topology, Fano-limited resolution). Held as
  a Phase-4 one-time validation of the 3–10 keV response; Geant4 deexcitation
  + Heed `TransportPhoton()`/W+Fano covers the campaign need.

**VMM3a modeling (Stage C)** *(chip: 64 ch, 3rd-order DDF shaper, peaking
25/50/100/200 ns, gains 0.5–16 mV/fC, 10-bit peak ADC (~8-bit effective),
neighbor logic; Iakovidis et al., NIM A 955 (2020) 163306)*. Two published
response models exist: the **ATLAS Athena `VMM_Shaper` closed-form h(t)**
(pole/residue form, ~20 lines of code, includes peak detection, ballistic-
deficit scale factor vs the ~150 ns ion tail, neighbor-logic emulation), and
the ESS/CERN CR-RC³ semi-Gaussian proxy (Backis, Piscitelli, Pfeiffer et
al., arXiv:2510.01981). **Charge-integration per pad + threshold is a
defensible baseline** (at t_p = 200 ns essentially the whole MM signal is
integrated, so peak ≈ ∝ collected charge) with three required corrections:
(i) threshold set in ENC units on the shaped amplitude — noise at the P2 pad
capacitance **75–170 pF** sets the floor (several-thousand-e⁻ ENC → few-fC
practical thresholds; BASKET expects ~21 fC signals, range 62 fC–2 pC);
(ii) ballistic-deficit factor if t_p < 150 ns is ever modeled;
(iii) ~~neighbor logic~~ — **re-scoped 2026-08-06, see below**.
Also simulate explicitly: **ADC saturation** — at high gas gain a 6–8 keV
photoabsorption (~250 e⁻ point deposit) can saturate the ~8-bit range while
MIP pads stay in range; saturation interacts with the discrimination cut.
Full waveform simulation with the Athena h(t) is Phase 4, needed only for
time-walk modeling beyond the (PDO, TDO) pair or streaming pile-up studies.
**Implementation started (2026-08-05): `vmm/`** — Athena shaper ported to
Python and unit-tested, channel emulator with T0 (charge+threshold) and T1
(shaper + neighbor logic + ENC + ADC) tiers, Stage-B interface contract and
NEEDS-DATA list; see `vmm/README.md`.

**Neighbor logic is a strip feature and mostly does not apply here**
*(2026-08-06, from Dylan's operating input + `vmm/nl_map.py`)*. Two facts:

- **The charge cloud fires ~1.1 pads on average.** Inter-pad clustering is a
  ~10 % minority effect — real, and to be used optimally, but not a design
  driver. It is *not* the multi-strip cluster that NL, µTPC and cluster-size
  cuts were designed for; anything inherited from ATLAS NSW that assumes a
  several-channel cluster has to be re-derived, not transplanted.
- **NL fires on chip-channel neighbors, and on a pad plane those are not
  spatial neighbors.** Measured over all 1280 pads: in the 9-column mapping
  revision chan ±1 is a touching pad 73.8 % of the time (median separation
  11.9 mm) and covers 25 % of the 5.78-pad physical neighborhood; in the
  7-column revision it is **7.7 %** and **126 mm median** — NL there would
  triple the data volume reading pads across the wedge. P(the actual
  charge-sharing partner is in the NL set) = **36 % / 2 %** respectively.

Consequences adopted for the campaign: **baseline is NL off** (the experiment
does not plan to use it, and the upside is bounded at ~10 % of tracks × 36 %
coverage ≈ 3–4 % of tracks gaining a partial second pad, against ~3× channel
occupancy); **clustering is an offline geometric step** on pads that each
passed threshold, which works with NL off and is where the real inter-pad
sharing gets used; and if NL is ever enabled it must be simulated in
*channel* space with the revision stated. Modeling NL as "geometric
neighbors" — the placeholder in the original plan — simulates a detector we
do not have and would bias recorded cluster size in both directions.

**Timing is now a Stage C deliverable, not a Phase-4 option** *(2026-08-06 —
"we are very interested in the time resolution of VMM and would definitely
want to compare what we get at SPS with simulation")*. First quantitative
pass done: `vmm/time_resolution.py` +
`research/TIME_RESOLUTION_NOTES.md`. Predicted per-pad σ_t (3 mm, gain 10⁴,
2 fC, ENC 3 k e⁻, threshold-crossing timestamp, walk-corrected with the same
hit's PDO): **10.6 ns Ar/Iso, 12.1 ns Ar/CO₂/Iso, 15.0 ns Ne/Iso 90/10,
16.2 ns Ne/CH₄, 25.0 ns Ne/CO₂/Iso 95/3/2** at t_p = 100 ns. The result is
**σ_t ≈ (1.0–1.5)/(n_p·v_d)** — ionization statistics, not electronics:
ENC jitter, TAC quantization, Polya and diffusion together contribute < 3 ns
in quadrature. The cost is one toy run per configuration (~1 min), so σ_t
can be carried through the Phase-3 matrix as a cell value rather than
computed once. Note **µTPC is unavailable at 1.1 pads/hit** — one timestamp
per hit, so gas/threshold/peaking time are the only knobs.

## 3. Phase 0 — preparation (code + infrastructure)

Do now (none of it depends on Alexandra's answers):

| # | task | notes |
|---|---|---|
| P0.1 | **Fix gas mixture composition bug** ✅ **done 2026-08-05** (see `research/GAS_FIX_NOTES.md`) | `makeMix2/3` feeds volume fractions to `AddMaterial()` which takes **mass** fractions. Negligible-ish for Ar/Iso (5 % vol = 7.3 % mass) but fatal for Ne mixes (Ne/Iso 95/5 vol = 86.5/13.5 by mass). Rewrite builder to take vol% + component molar masses → mass fractions; recompute densities from ideal-gas mixing. Also verify the 2.67e-3 g/cm³ isobutane density (ideal-gas at 20 °C gives 2.42e-3). |
| P0.2 | **Add campaign gases + W values** ✅ **done 2026-08-07** | All seven §7.4 mixtures are built and run: `ArCO2Iso9352`, `NeIso9010`, `NeIso8020`, `NeCO2Iso9532`, `NeIsoCO29055`, `NeCH4937`, `NeC2H69010` (+ `ArIso`). Composition, density, mass fractions and W now come from **one table** (`include/GasMixtures.hh`) instead of a composition list in `DetectorConstruction` and an unrelated hand-maintained W map in `SteppingAction` — the two used to be keyed by the same string with nothing checking they agreed, i.e. the same class of error as the 2026-08-05 mass/volume-fraction bug. `mm_sim --list-gases` prints the table; an unknown gas now fails at argument parsing with the list rather than core-dumping after geometry construction. **Three deviations from this task as written:** (i) mixture W is **stopping-power weighted**, not the harmonic *volume* weighting specified here — the correct weight is the energy fraction each species absorbs, which scales with molecular electron count, and the two differ by ~7 % for Ne/iC₄H₁₀ 90/10; (ii) CF₄'s W is a CLI parameter (`--w-cf4`) so the 35–52 eV spread is bracketed as run points rather than hardcoded, as asked; (iii) **Ne/CO₂/N₂ 90/10/5 was not added — those fractions sum to 105 %.** Someone needs to say what the intended mixture is before it can be built. |
| P0.2c | ⬜ **Re-source the pure-gas W values against PDG** *(new 2026-08-07)* | The component W table in `src/GasMixtures.cc` is inherited from the old `SteppingAction` map and has **not** been checked against a primary source. Isobutane is carried at 26.0 eV where the commonly quoted value is nearer 23.4 eV — a 10 % shift that would propagate into every `nPrimary` in every argon run. Deliberately left unchanged rather than silently moved, so previously produced numbers stay explicable. Cheap; do it before Phase 2 and record a source per component. Note the ordering: this matters *less* than Penning transfer, which is not modelled at all and makes every Ne-mixture W an upper bound. |
| P0.3 | **Gun angle control** ✅ **done 2026-08-07** | `--gun-theta` (deg from the wedge normal), `--gun-phi` (tilt azimuth), `--gun-standoff`. The beam **pivots about the aim point at the drift mid-plane**, so an angle scan re-illuminates the same pads instead of walking across the wedge with θ; the standoff auto-raises to clear the front window bulge at large θ. Verified against a 30° run: primary-track slope dx/dz = 0.5779 vs tan 30° = 0.5774, aim point stable to 16 µm. Unblocks the §4.1/§4.2 angle scans. `--beam-spread` also landed (2026-08-07). **The aim point has to be chosen against BOTH polar coordinates.** The historical default (r = 355 mm) sat 35 µm from a radial pad boundary; moving it to a ring centre then put it at φ = 30.0000°, which is 61 µm from the centre of a 126 µm *azimuthal* inter-pad gap. Measured with the real pad artwork: only **25.7 %** of events deposited in pad copper at all, against **100 %** at a true pad centre (⟨edepPadCu⟩ 4740 eV vs 20383 eV). The default is now a pad centre in both coordinates (305.385, 169.399), and `mm_sim` checks it against the gerber-derived table in `include/P2PadMap.hh` rather than a radial heuristic — it now catches the radial gap, the azimuthal gap, and off-pad-field aim points separately. |
| P0.4 | **Cluster provenance in output** ✅ **done 2026-08-05** | per cluster: creator process + logical volume where the depositing track (or its ionizing ancestor) was born. This is what turns photon runs into a per-layer conversion budget. Add per-event "first interaction volume" for γ runs. **Needs a `TrackingAction` + `G4VUserTrackInformation`, neither of which exists in the repo yet** (`ActionInitialization` registers only primary/run/event/stepping) — one-level `parentID` is not enough, since a δ-ray of a photoelectron reads as `eIoni`. Carry {origin process, origin volume, origin energy, ancestor trackID} and copy it to every secondary at creation. |
| P0.4b | **Gas-crossing track-segment export** | per charged track crossing DriftGas: entry point, direction, momentum, PID (a small `SegmentTree`). This is the Heed handoff for Stage B's high-fidelity mode (re-ionize with `TrackHeed::NewTrack()`), per the Pfeiffer-paper recommended split. Cheap to add now, enables the threshold studies later. |
| P0.4c | **Exit-particle export** | per event: every particle leaving the wedge envelope (PID, energy, direction, exit surface, creator process + birth volume, ancestry link to the primary). In photon runs this bounds **correlated multi-wheel fakes** (a Compton-scattered 100–150 keV photon keeps most of its energy and can convert again in wheels 2/3 — accidental-rate arithmetic misses this entirely) by pure wheel-to-wheel convolution, no 3-wheel geometry needed yet; it is also the input the eventual full 3-wheel sim (§8) will want. Decision 2026-08-05: save full interaction histories — track everything that produces signal, where it came from, and how it interacted. |
| P0.5 | **`--skip-empty` output flag** ✅ **done 2026-08-05** | photon runs are ≥99.9 % empty events; write only events with ≥1 cluster (keep the total-thrown count in metadata for normalization). |
| P0.6 | **PAI model in gas regions + deexcitation** 🟡 **deexcitation + fine-cut region done 2026-08-05; PAI still open** | G4Region over DriftGas+AmpGas with PAI/PAIphot (~2× EM cost, thin-layer straggling fixed); switch physics list to Livermore-class EM with fluorescence/Auger/PIXE on (needed for Cu-K/Fe-K fluorescence in photon runs). CLI switch to compare with default EM (one validation study early in Phase 1). **⚠ Includes a bug fix, do not skip:** `PhysicsList::SetCuts()` sets a **0.1 mm γ production cut**, which in copper is ~5–10 keV — i.e. the **8.05 keV Cu-K fluorescence line is at or below threshold and is probably being suppressed today.** That is the dominant *gas-sensitive* wall channel (§1.1). Set `G4EmParameters::SetDeexcitationIgnoreCut(true)` and add a fine-cut (~1 µm) `G4Region` over mesh + pad Cu + gas; the current 10 µm e⁻ cut is ~1 keV in gas but ~80 keV in Cu. |
| P0.7 | **Condor scripts refresh** | `submit_condor*.py` still speak MX17 modes; teach them the p2 flags + campaign manifest (one JSON/CSV row per run point). Also `collect_results.py`, which still parses MX17-style `{gas}_{particle}_{E}MeV` filenames and will not recognise the §9 naming scheme. **Half of this got easier on 2026-08-07** — see P0.21: the manifest fields the submit script would have had to record are now written into the ROOT file itself, so the script's job is to schedule and name runs, not to be the sole custodian of what a file is. |
| P0.21 | **Run provenance inside the output file** ✅ **done 2026-08-07** | New `RunMeta` tree, one row per worker (`hadd` concatenates; `thrown` sums). Carries the **thrown count** — the normalization denominator that previously existed *only* in stdout, making any file separated from its job log unnormalizable and, strictly, unidentifiable — plus the code git hash and a **dirty flag**, a hash + readable digest over every geometry-affecting parameter, the resolved gas composition/density/W, and the full beam configuration. `collect_results.py: read_run_meta()` sums thrown across workers and **warns** when rows disagree on git or geometry hash (an accidental merge across two run points) or when a run came from a dirty tree. This is what makes §9's "no un-manifested runs" enforceable rather than aspirational. Schema documented in `OUTPUT_FORMAT.md` §0. |
| P0.8 | **First lxplus build + geometry validation** ✅ **done 2026-08-05** | compile, overlap check, 100-event smoke run, one event display. Can be done with provisional geometry — worth doing *before* Alexandra's answers so her changes land in verified code. |
| P0.9 | **Stage B skeleton** ✅ **skeleton done 2026-08-07** | `stage_b/` runs end to end on real Stage A output: drift depth from the recorded frame → attachment → drift time → longitudinal and transverse diffusion → mesh transparency → Polya gain → polar pad lookup → per-pad charge and earliest arrival, written as a `PadHitTree` plus a `StageBMeta` tree carrying every assumption the run was made under. 16 closed-form self-tests (`python3 -m stage_b.selftest`) run without Geant4, so a Magboltz table is a genuine drop-in: `TransportTable` is the only interface to implement. **First result: 1.070 pads/event** (373 single / 26 double / 1 triple, 400 muons, Ar/iso, 3 mm) — which *independently reproduces* the ~1.1 pads/hit figure previously derived from the geometric analysis in `vmm/nl_map.py`. Two independent routes to the number the NL-off baseline rests on. Median pad charge 37.6 fC, drift span 0–60 ns as the placeholder predicts. **Four placeholders remain and are recorded in every output file**: gas transport (P0.10), induced current (P0.17), pillars (P0.16), constant mesh transparency (P0.13) — see `stage_b/README.md` §3 before quoting anything. *(On lifting MX17's digitizer: their package is real and good, but it is built around a 30 mm gap, 512×512 square pads, a resistive stack and DREAM. Our pad plane is polar, our gap is 3 mm and our charge cloud is ~1 % of a pad, so the parts that would have been lifted verbatim — the induction kernel and the electronics — are exactly the two plug-ins that differ. What was taken instead is their **hard-won knowledge**: the drift-frame trap, nPrimary-not-edep/W, the packet-vs-electron tradeoff, and recording assumptions in the product. Their §2.2 offer stands for the induction kernel when P0.17 happens.)* |
| P0.10 | **Magboltz gas tables** 🟡 **dry tables done 2026-10-08** for the six October-2026 gases (ArIso, ArIso9010, NeIso, NeIso8515, ArCO2Iso9352, ArCF4Iso; 40 E points 100 V/cm–90 kV/cm, ncoll 10), loaded in Stage B as `MagboltzTable` (default, `--drift-field`). **Still open: wet variants and Penning rates (needed for the gain).** | generate on lxplus per gas (log E-grid spanning drift ~100 V/cm–1 kV/cm and amplification ~30–60 kV/cm; ncoll = 10), `WriteGasFile` + ion mobility files, cache in repo or EOS. Order an hour per gas per core; embarrassingly parallel. Garfield++ comes from cvmfs LCG views on lxplus. Enable Penning transfer per §2 (manual r for Ne/iC₄H₁₀ and ternaries — document the chosen value and bracket). **Generate wet variants (0.5/1/2 % H₂O) alongside dry** *(MX17, 2026-08-06 — §2.4)*: MX17's bench finds water contamination dominates drift velocity — 36.6 µm/ns measured against a far higher dry prediction in Ar/iso, 1–2 % H₂O inferred at the SPS. Cheap (same job, one extra component) and it matters directly: σ_t ∝ 1/v_d, so a dry-only table can be optimistic by tens of percent against the SPS data we are trying to match. Quote a dry→2 % bracket rather than a line. |
| P0.11 | **Source mode** | isotropic point/disc source gun with line energies — ⁵⁵Fe 5.9/6.5 keV, ²⁴¹Am 59.5 keV, ¹⁰⁹Cd 88 keV, ⁵⁷Co 122/136 keV — plus the ⁹⁰Sr/⁹⁰Y β spectrum (CSV already in `sr90_calibration/`); configurable standoff. Powers the §5a source-validation predictions and the ⁵⁵Fe-vs-⁹⁰Sr timing bench test (`research/TIMING_PSD_NOTES.md`). |
| P0.12 | **Per-event interaction classification + edep by region** ✅ **done 2026-08-05** | the branch that makes the requested money plot a one-liner. Per event: `interactionClass` ∈ {none, gas-photoelectric, gas-Compton, gas-Rayleigh, wall-photoelectric, wall-Compton, fluorescence-reabsorption, MIP-crossing, other}, `conversionVolume`, `conversionZ`, plus **edep in every region, not just DriftGas/AmpGas** — `FrontGas`, `DriftCathode_Gas` and `BackGas` are currently unscored, so "converted somewhere that produces no signal" is indistinguishable from "did not convert". Also record amp-gap cluster z *relative to the mesh* explicitly (Stage B needs it for partial gain, §2). Then "event-by-event edep separated by interaction type" is one `TTree::Draw` per gas. Builds on P0.4's TrackingAction; do them together. |
| P0.13 | **Cylindrical-wire mesh cross-check** *(downgraded 2026-08-05 — see note)* | the mesh **is** the largest single fake source (measured 1.98×10⁻³ per 60 keV photon, 9× the argon gas rate), but the existing effective-density slab is **not** wrong about the rate: its areal density π·d²/(2·pitch) = 8.46 µm solid-equivalent is *exactly* the plain-weave value (verified). Only the escape geometry is approximate, and the diluted-slab and solid-wire estimates differ by ~1.4×, not 3×. So this is a **cross-check worth doing, not a blocker**. |
| P0.14 | **Analytic cross-check gate** ✅ **done 2026-08-05** | `scripts/photon_budget.py` reproduces the §1.1 targets; the first Geant4 photon run matched it (3.4×10⁻³ measured vs 4.3×10⁻³ predicted total, mesh dead on). Re-run this comparison after the P0.1 gas-mixture fix, since a mass/volume-fraction error shows up here first. |
| P0.15 | **Phase-space gun** *(new 2026-08-05, from `NEEDED_INPUTS.md` §1)* | read primaries (position, direction, energy, PID, weight) from a file at the wheel plane instead of the fixed pencil beam, so results can be weighted by the real spatial/angular distribution of signal electrons and background photons. Cheap to write; **blocked on the collaboration supplying the file**, which is the single highest-value external ask we have. |
| P0.17 | ⬜ **Realistic induced current — static Ramo, not a Phase-4 project** *(promoted out of §8 on 2026-08-07; MX17 §2.1)* | Stage B currently feeds the VMM shaper **delta charges** (Athena's model). The real Micromegas signal is a fast electron spike plus a ~150 ns ion tail, and that shape sets the leading-edge slope — the main model risk in the σ_t prediction and the only reason the time-at-peak estimator is untrustworthy. **P2's stack is non-resistive, so the weighting potential is static**: no time dependence, no FEM, no solver. Riegler's closed-form layered solutions (JINST 11 (2016) P11002) are already in Garfield++ `ComponentParallelPlate` (`AddPixel`/`AddStrip`). Replace each electron's delta charge with `i(t) = Q_e·δ_fast(t) + Q_ion·i_ion(t; µ_ion, gap)` weighted by the pad's Ramo potential. **Effort: an afternoon.** Also retires `MM_ION_FLOW_TIME_NS = 150.0` in `vmm/vmm_shaper.py` — an ATLAS NSW gap/gas constant that should be *dropped* once real induction exists, not retuned. Details: `research/TIME_RESOLUTION_NOTES.md` §4.1. **⚠ Read `HANDOFF_MX17_RESPONSE.md` §2.6 before implementing.** MX17 measured that grounding the inter-pad gaps in the weighting solve — which is exactly what `ComponentParallelPlate::AddPixel` assumes, since it treats the pixel as a patch of a continuous conductor — costs **+27.2 % on prompt captured charge** and invents a 4.5× sub-pad amplitude swing. It is a DC error that does *not* decay with distance from the plane. **Checked for P2, 2026-08-08: our gap fraction is 2.13 %** (metal fraction 0.9787, computed over all 42 rings from the gerber pad table in `include/P2PadMap.hh`; independently 0.9781 from that file's own header). Against MX17's 24 %, both of §2.6's scalings predict a **1.8–2.4 %** error for us — §2.6's own "≲2 %, fine at the percent level" branch. The polar geometry is what saves us: ~11 mm pads with ~126 µm gaps. So this route stands, with the ~2 % under-estimate of absolute captured charge quoted as a known-sign systematic. Note the affected quantities are *not* the ones P0.17 exists for: §2.6 reports ratios between pads and arrival times as barely affected (peak kernel +0.2 %), and σ_t and the time-at-peak estimator are what we are fixing. Revisit with MX17's `response/solver/v6_pad_gaps.py` only if absolute normalization or a sub-pad efficiency map becomes a deliverable. |
| P0.18 | ⬜ **Zone the readout copper — P2 has this error today** *(new 2026-08-07; MX17 §2.4 / `NEEDED_INPUTS.md`)* | `SimConfig.hh` documents `p2_fcu_coverage`/`p2_bcu_coverage` as measured **over the active area**, but `EffCu` applies them to a copper slab spanning the **whole wedge board** — so the model puts near-solid copper out in the fan-out and periphery, where real coverage is far lower. MX17 had the mirror-image bug (one board-wide average everywhere: 26 % too little copper in the active area, ~10× too much outside) and fixed it with two volumes and no measurable CPU — see `BuildReadoutZone` in `MX17_Geant/shared/MX17ModuleGeometry.hh`. **Why this is not cosmetic:** copper is ~90 % of the photon fake budget in argon (§1.1, and the first Geant4 photon run agreed), so mis-placing it distorts the conversion-layer budget that the entire gas decision rests on. **MEASURED 2026-08-07** (gerbers, exact vector geometry): inner margin r 95–115 mm F.Cu **0.123**; pad field r 115–594.9 mm F.Cu **0.913**, B.Cu 0.171; fan-out r 594.9–650 mm F.Cu **0.180**, B.Cu 0.235; whole board F.Cu **0.787** (was 0.937) and B.Cu **0.180** (was 0.181). So the old board-wide F.Cu sheet carried **25 % more copper than the artwork has**, most of it dumped in the inner margin and fan-out where there is almost none — while B.Cu's 0.174 was never materially wrong, only `P2_GEOMETRY.md`'s prose. Zoning is now implemented as 10 radial bands. **The coverage script had to be fixed first**: `analyze_cu_coverage.py` mis-measures F_Cu because `ApertureDef.size` returns `max(params)`, which for KiCad `RotRect` macro apertures picks up the *rotation angle* (e.g. `RotRect,0.3X1.8X302.763` → "size" 302.763 mm), painting ~151 mm-radius discs at the ten connector footprints; it also counts `NonConductor`/`Profile` apertures as copper. So `p2_fcu_coverage = 0.983` is itself inflated and must not be used as-is (`NEEDED_INPUTS.md`, verified 2026-08-07). |
| P0.19 | ⚠ **Real pad/strip artwork instead of a homogenized slab — CONTESTED, decide before building** *(new 2026-08-07; MX17 §2.4 item 2)* | The original suggestion: if the artwork is periodic it can be built as real geometry with nested `G4PVReplica` — MX17's 786 432-feature pattern costs **15 extra volumes and ~3 % CPU**, not 786 432 placements. It only affects backscatter off the board (the signal copper is downstream of both gas gaps), so it already ranked below P0.18. **But MX17's own plan has since concluded the opposite for their detector**: `RESPONSE_SIM_PLAN.md` §6 item 4 reads *"Do NOT model strips/pads/coverlay as Geant4 volumes — material budget is unchanged at the level that matters and the response chain owns that geometry."* That is a direct reversal of the advice we were given, from the people who tried it. **Before more effort goes in, someone has to decide which applies to P2.** The arguments do not transfer automatically: P2's drift gap is 3 mm against MX17's 30 mm, so a backscattered electron from the board is a much larger relative perturbation here, and P2's whole campaign question is a *fake-hit budget* dominated by wall conversions — the case MX17's "material budget is unchanged" claim is least likely to cover. Counter-argument for P2 is therefore real, but it is an argument, not a measurement: **the cheap resolution is to build it once, compare the conversion-layer budget against the homogenized slab, and keep the pattern only if the difference is visible.** *(Implemented 2026-08-07: 1280 annular-sector pads on the gerber-exact polar grid, 42 `G4PVReplica` rings, `--homogenized-readout` to switch back. 0 overlaps in both modes.)* **FIRST COMPARISON RUN, 2026-08-07 — split the question in two, as the implementing session proposed:** (a) **zoning** (one board-wide average → 10 radial bands) is *not* contested — the old F.Cu number was simply wrong, 0.937 against a measured 0.787, i.e. 25 % too much copper, most of it dumped where the artwork has almost none. (b) **patterning** (band-average sheet → real pads) is the part MX17 §6 item 4 argues about, and the measurements so far are: **muon ionization — no difference.** 2500 muons each: ⟨edepDrift⟩ 925 ± 28 eV patterned vs 937 ± 30 homogenized, medians matching to 0.4 eV; pad multiplicity 1.120 vs 1.117. (An apparent 8.6 % gap at 400 events was pure Landau statistics — a warning about small-sample comparisons here.) **Photon conversion budget — suggestive, not resolved.** 400 k γ at 60 keV each: total gas-ionization rate 3.500 ± 0.094 ×10⁻³ vs 3.378 ± 0.092 ×10⁻³ (0.9σ). But the difference is **confined to the one volume that changed**: `PCB_Cu_F` gives 1.025 ± 0.051 ×10⁻³ patterned against 0.918 ± 0.048 ×10⁻³ homogenized (**+12 %, 1.5σ**), while `Micromesh` (2.07×10⁻³, +0.1σ), `DriftGas` (+0.1σ) and every other layer are identical. A coherent story in the right place, at 1.5σ — which is exactly the level at which one should not conclude. **EXPLAINED, and the ~4 M-photon run is cancelled as uninformative.** The +12 % is not escape physics — it is **local areal density**. The bands are radial-only, so band 4 (r 294.95–354.94) carries fCu = 0.906 measured over *board ∩ annulus*, which includes the azimuthal margins of the wedge where there are no pads at all. The pad field itself is 0.978 of its annulus. The beam sits in the pad field, so the homogenized sheet simply has less copper under it. Quantitatively, for the `--beam-spread 11.4286` run above the beam spans bands 4 and 5 with weights 0.804/0.196 → homogenized fCu under the beam = 0.9081, against a patterned fill of 0.9781. **Predicted ratio 1.077 (+7.7 %), measured 1.117 ± 0.081 — 0.49σ apart.** Independently confirmed with muons: ⟨edepPadCu⟩ MPV ratio at a pad centre = 1.0959 over 6000 events, against 1/0.906 = 1.104 predicted. Thin-layer conversion probability is linear in areal density, so this fully accounts for the observation with no escape physics invoked. **And escape is bounded below the level that matters.** CSDA range in g/cm² is density-independent, and so is layer thickness in g/cm² (pad 0.0161, sheet 0.0146), so forward/backward escape fractions are nearly identical. Sideways escape into the etched grooves affects only perimeter × range / area = **0.35 %** of a pad (a 60 keV photoelectron has ~10 µm CSDA range in Cu on an 11.3 × 11.8 mm pad). Resolving a sub-percent effect at 3σ would need ~10⁹ photons — so escape is settled analytically and should not be simulated. **The corrected conclusion, which sharpens P0.19 rather than closing it:** radial-only banding is *too coarse where the beam is*. It dilutes pad copper into azimuthal margins the beam never illuminates, under-stating the local conversion rate by ~8–10 % at any pad centre and over-stating it in the margins. So: **board-integrated quantities agree** (total gas-ionization 0.9σ, pad multiplicity 1.120 vs 1.117, muon drift ionization identical), while **anything local to the beam spot differs by ~10 % with a known sign** — a systematic, not a fluctuation, and one the homogenized model cannot fix by re-measuring. Fixing it would require azimuthal zoning, at which point the pads have been re-implemented badly. That is a stronger argument for keeping the pattern than the original 1.5σ was. **Confirmed by a uniform-illumination run (2026-08-07):** with the beam spread over the whole wedge (500 k γ each, 250 mm spread), patterned and homogenized agree to **0.0σ** on `PCB_Cu_F` (7.700 vs 7.680 ×10⁻⁴) and −0.3σ on the mesh — exactly as the areal-density picture predicts, since the two models carry the same total copper mass by construction. So the effect is *purely* local: ~10 % where the beam sits on the pad field, zero when averaged over the board. **And the back copper is irrelevant either way:** `PCB_Cu_B` produced **0 conversions in 900 k photons** across both illumination patterns. It sits behind 18 µm Cu + 200 µm FR4, which a 60 keV photoelectron cannot traverse to reach the gas. So B.Cu staying a homogenized sheet (its 30 789 routed traces are not periodic and cannot be replicated) costs nothing for the fake budget — it matters only as material budget, where the zoned sheet has the mass right. |
| P0.20 | ⬜ **Resistive layer: confirm it is actually uniform** *(new 2026-08-07; MX17 §2.4 item 3)* | `DetectorConstruction.cc` places `ResistivePaste` as a **100 µm slab of full-density (1.4 g/cm³) paste** over the whole wedge — no coverage scaling, no structure. MX17's is now 515 discrete 550 µm ESL strips on a 0.8 mm pitch with real chamber gas in the grooves. This matters more than the pad pattern: the resist is the first solid the avalanche region sees and at 100 µm it is ~5× a P2 copper layer. **Needed input:** is the P2 resistive layer strips, a uniform DLC/paste sheet, or pads? If uniform, the current slab is right and only thickness/density need confirming — but that should be **stated, not assumed**. |
| P0.16 | 🟡 **Get the Dynamask pillar map — before any production run** *(new 2026-08-06, Dylan)* **Geant4 part done 2026-10-06:** `--pillars saclay` (default; V1 insulation mask, det1–det4, 3.80 %) or `cern` (det5, 4.90 %), placed in the amp gap from `include/P2Pillars.hh`, which `scripts/gerber/extract_pillars.py` generates with P2_basket_analysis' own `load_pillars`. **Still open: pillar dead spots in Stage B** (electrons landing on a pillar get no gain). | Fetch it from **`P2_Basket_Analysis`** (or wherever it lives) and put pillars in the geometry. We already have the *design* pattern from the CERN bulk mask (`P2_Mask2.gbr`: Ø 0.5 mm, pitch 2.000 mm exact, 41 366 pillars, **4.8 % of the amp gap**) and the material (**Dynamask** dry film, ρ ≈ 1.2–1.4 g/cm³ — not kapton/FR4); what we want is the map the *analysis* actually uses, so sim and data share one definition (same argument as the mapping revision, §10.6). Feeds: pillar dead spots in Stage B (electrons landing on a pillar are lost, no amplification there), the 4.8 % effective-gas correction, material budget. **Gating**: retrofitting pillars after a production run invalidates per-pad efficiency and the amp-gap scoring. Details: `NEEDED_INPUTS.md` addendum. |

Blocked on Alexandra (re-run affected points if answers move defaults): mylar
thicknesses/aluminisation, back-frame depth, frame opening/material,
overpressure/sag, mesh spec (wire diameter/pitch → transparency model),
pad-mapping revision.

## 4. Phase 1 — Geant4 sensitivity scans (MX17-style)

Goal: primaries produced in the drift volume (and elsewhere) per incident
particle, vs energy — the direct analog of the MX17 sensitivity plots, now
per gas and with provenance. Runs at baseline geometry (3 mm gap) plus the
gas list; **deliberately wider energy ranges than strictly needed**.

### 4.1 Electron scan

- **Energies** [MeV]: 1, 2, 5, 10, 20, 30, 50, 75, 100, **119**, 130, **144**,
  155, 175, 200 (log-ish below 30; dense around the elastic lines).
- **Angles**: 0°, 10°, 20°, 30°, 40° (full set at 30/119/155 MeV; {0°, 10°}
  elsewhere).
- **Position**: default mid-active (r = 355 mm); position axis deferred to a
  small dedicated scan (inner r = 150, outer r = 550, near-edge φ) at
  baseline energy only.
- **Statistics**: 10⁴/point general, 10⁵ at the baseline points (Landau
  tails, zero-cluster tail).
- **Observables**: nPrimDrift distribution (mean, width, P(0), P(<n));
  edep spectrum; cluster z-profile; δ-ray fraction reaching pads/amp gap;
  same for AmpGas (direct amp-gap ionization = anomalously large signals).

### 4.2 Photon scan

- **Energies** [keV]: 3, 5, 8, 10, 15, 20, 30, 40, **50**, 60, **80**, **100**,
  120, **150**, 200, 300, 500. (3–10 keV covers fluorescence/K-edge
  structure: Ar K 3.2 keV, Cu K 8.05 keV, Fe K 6.4 keV from the mesh; bold =
  the dominant background band.)
- **Angles**: 0°, 10°, 30°. **Statistics**: 10⁶–10⁷/point (with P0.5 the
  output stays small; most photons pass through untouched).
- **Observables**:
  - P(≥1 primary in drift gas) vs E — *the* photon sensitivity curve;
  - **conversion-layer budget — now the decision variable of the whole gas
    question** (§0, §1.1): fraction of gas-ionizing events whose electron
    originated in {drift gas, cathode mylar/Al, mesh, pad Cu, FR4, windows,
    front/back gas} — per gas, this separates "gas choice helps" from
    "wall floor". Run it **both** with the effective-density mesh slab and
    with the P0.13 wire geometry, and quote both; the analytic estimate
    (mesh 2.0×10⁻³ vs argon gas 2.3×10⁻⁴ per 60 keV photon) says this
    single number decides Ar vs Ne;
  - **per-event edep spectrum split by interaction class** (P0.12) — the
    plot this scan exists to produce. Overlay, per gas and per photon
    energy: gas-photoelectric / gas-Compton / fluorescence-reabsorption /
    wall-photoelectric / wall-Compton, against the electron-run MIP Landau
    from §4.1 on the same axis. **Rayleigh must be counted separately and
    excluded from "interactions"** — it deposits nothing;
  - nPrim spectrum *given* a conversion, split the same way (this is the
    fake-hit charge spectrum seed for Phase 2). Expect gas conversions at
    3–30× the MIP deposit and wall photoelectrons at ~3–5× — i.e. the wall
    population sits *inside* the MIP Landau tail, which is the pessimistic
    case for a charge cut and must be measured, not assumed;
  - **containment check**: fraction of the conversion electron's energy
    actually deposited in the gap vs its initial energy. §1.1 predicts
    Compton electrons are contained below ~80 keV while photoelectrons are
    not; if the sim disagrees, the range/geometry model is wrong;
  - escape-peak / fluorescence structure in Ar vs Ne mixes — including the
    **Cu-K 8.05 keV re-absorption chain** (§1.1), which needs P0.6's
    `deexcitationIgnoreCut` fix to exist at all;
  - **exit-particle inventory** (P0.4c): energy/direction spectra of
    photons and electrons leaving the wedge per incident photon — the
    convolution input for the correlated multi-wheel fake estimate
    (§5 step 5).

### 4.3 Muon validation runs — SPS test-beam anchor *(added 2026-08-05)*

An SPS test with **200 GeV muons** has just been taken (data analysis in
progress); `TESTBEAM_PLAN.md` is the full comparison plan. This subsection
adds the sim-side run points to Phase 1 so they ride the same
infrastructure and bookkeeping:

- **Particle/energy**: μ (sign per as-run beam), 200 GeV — CLI is already
  capable: `-p muon -e 200000` (energy flag takes MeV; `muon+` exists for
  μ⁺). **Normal incidence** = the existing default gun direction (+z
  pencil beam), aimed via `--gun-x/y` at the pad region illuminated in the
  test; no P0.3 angle work needed for this point. A small beam-spot/
  divergence option gets added once measured beam parameters arrive
  (`TESTBEAM_PLAN.md` §3).
- **Configurations**: (i) the **as-run test-beam configuration** (gas, HV,
  geometry deltas — from the `TESTBEAM_PLAN.md` §1 conditions checklist)
  for the direct data comparison; (ii) the **campaign reference gases**
  (§7.4 list) at baseline 3 mm, so the test-beam-tuned parameters can be
  transferred per gas and the tuned-vs-untuned shift quoted.
- **Statistics**: 10⁵ general; 10⁶ where efficiency-vs-threshold is
  compared (plateau inefficiency is a ≲1 % effect).
- **Observables**: same set as §4.1 (nPrim, edep, z-profile, δ-ray
  fraction, AmpGas direct ionization), plus through Stage B/C the
  data-comparable set: Landau MPV/width, pad multiplicity, efficiency vs
  threshold/HV (`TESTBEAM_PLAN.md` §4 has the observable↔parameter
  matrix).
- **Physics note**: 200 GeV μ sits on the relativistic-rise plateau —
  same MIP regime as our 30–160 MeV electrons in mm of gas, so the anchor
  transfers to the electron-efficiency side of the ROC nearly directly;
  it does *not* constrain the photon side. Keep muon radiative processes
  on (rare hard events populate the charge-spectrum tail seen in data).

### 4.4 Phase-1 deliverables

Sensitivity figure pair per gas (electron + photon curves overlaid across
gases), the conversion-layer stacked-bar figure, and a short note fixing the
nPrim→charge conventions for Phase 2. Also the PAI-vs-default-EM comparison
(P0.6) at 119 MeV / 3 mm / ArIso + NeIso.

## 5. Phase 2 — response chain and ADC discrimination

Apply Stage B+C to Phase-1 baseline samples (3 mm, all gases; 119 MeV e⁻ at
10° vs 50/80/100/150 keV γ):

1. Gas tables per gas (P0.10); pick the drift-field operating point per gas
   (near velocity plateau where one exists; document choice — this also
   yields the **drift-time span** = coincidence-window input).
2. Gain: set per-gas mesh voltage for a common reference gain (e.g. 10⁴,
   anchored per §2); record the voltage — operational headroom matters for
   stability/discharge arguments later.
3. Produce per-event **pad charge maps**; per event: max-pad charge, summed
   charge, pad multiplicity, hit-pad pattern length (a crude "trackness"
   proxy for oblique electrons). Note the discrimination is **topological as
   much as total-charge**: a 6–8 keV fluorescence photoabsorption is a
   point-like ~230–300 e⁻ deposit on 1–2 pads, while a MIP crossing 3 mm
   leaves ~30–50 e⁻ *spread along the track* — so diffusion and pad sharing
   must be modeled honestly, and both per-pad and per-event
   charge cuts should be scanned. High-energy wall photoelectrons ranging
   through the gas are the hard case (electron-like charge); flag their
   fraction separately. Expectation check (Dylan, 2026-08-05): at 10° a
   track averages only **~1.1 pad hits**, so topology is a *secondary*
   discriminant at the nominal angle — quantify what it adds, but the
   working assumption is charge deposition + 3-layer coincidence carry the
   weight. *(2026-08-06)* Clustering of the ~10 % of hits that do share
   charge is an **offline geometric step** (`vmm/nl_map.py`
   `geometric_neighbors`), run on pads that each passed threshold; the
   **VMM neighbor logic is off in the baseline** and is a separate,
   channel-space model (§2). Report pad multiplicity with the NL state
   stated — the two are not comparable numbers.
   **3a. Edep → charge proportionality audit** *(new 2026-08-05; the "does
   the deposited energy actually become ADC" question, full list in
   `research/PHOTON_DISCRIMINATION_NOTES.md` §5)*. Before any ROC is
   believed, quantify each proportionality breaker at the working point and
   state which ones are deposit-dependent (i.e. distort the discrimination)
   rather than a constant scale factor:
   - **z-dependent partial gain in the amp gap** —
     G(z) ≈ G_full^((z_mesh−z)/d_amp), four decades within one volume, and
     it governs the pad-Cu photoelectron population specifically. Modelling
     these at flat full gain is qualitatively wrong, not just imprecise.
   - **space-charge gain suppression on dense deposits** — a contained
     8 keV blob is ~300 primaries in ~100 µm; at gain 10⁴ that approaches
     the Raether limit and the gain sags. This **compresses exactly the
     separation the discrimination relies on**, and it is also the spark
     mechanism. No measured curve exists for our geometry ⇒ bracket it and
     flag any conclusion that depends on the compressed region.
   - **ballistic deficit** (§5 step 6) — breaks proportionality in our
     favour, but not against wall photoelectrons, which are track-like in
     time too.
   - **charge sharing + per-pad threshold**, **ADC saturation**, **pillars
     (4.8 % dead area)**, **mesh transparency**, **attachment** — all
     smaller; quantify once, then carry as a fixed systematic.
   - **Penning transfer / effective W** — gas-dependent, and it silently
     mis-scales Ne against Ar if a single tabulated W per gas is used.
     Also fix the `SteppingAction` W treatment: `floor(edep/W)` + uniform
     remainder gives Bernoulli rather than Fano fluctuations, which
     distorts the *width* the ROC integrates over.
   Deliverable: a one-page waterfall from ⟨edep⟩ to ⟨PDO⟩ per gas, with the
   deposit-dependent terms separated from the constant ones.
4. **ROC curves**: electron efficiency vs photon acceptance as a per-pad
   (and per-event max-pad) **charge window** scans, per gas — lower cut
   against noise and partial-clip deposits, **upper cut** against
   conversions (a gas conversion in the band deposits 3–30× the MIP charge
   — see §1.1; the "10–100×" figure applies only to fully contained few-keV
   photoabsorption such as ⁵⁵Fe — so a lower threshold alone cannot reject
   it; the ε_e price of the upper cut is the Landau/δ-ray tail, which the
   sim gives directly). **Break the ROC out by conversion class**: the
   dominant wall-photoelectron population lands at ~3–5× MIP, inside the
   Landau tail, and may be largely irreducible by charge alone — if so,
   say so plainly and let the coincidence carry it.
   Headline FOM: **ε_e at 90/99/99.9 % single-hit photon rejection**, and
   the window in fC (and in primary-electron equivalents) where it sits.
5. **Coincidence extrapolation** (analysis-only): with per-layer photon hit
   probability p_γ(E, threshold) and the experiment's 15–20 kHz/pad γ rate
   and drift-time window, estimate accidental 3-layer fake-track rates vs
   threshold — the bridge from single-wedge sim to the online-tracking
   question. Pure arithmetic once Phase 2 numbers exist. Add the
   **correlated-fake bound** from the P0.4c exit-particle record
   (Compton punch-through convolved wheel-to-wheel) — accidentals and
   correlated fakes scale differently with the coincidence window, so
   report both. **The window has a floor set by the per-hit time
   resolution** *(2026-08-06)*: it cannot be tightened below a few σ_t
   without losing real triples, and σ_t is 12 ns (Ar mixes) vs 25 ns
   (Ne/CO₂/Iso) — a 2× difference in window, i.e. **4× in accidental
   doubles and 8× in triples**, which is larger than the Ar-vs-Ne
   difference in photon conversion probability. Take σ_t per gas from
   `research/TIME_RESOLUTION_NOTES.md` §0 and quote the window as
   max(k·σ_t, drift-span term).
6. **Time-structure discriminants** *(blue-sky pass 2026-08-05 — details
   and the ⁵⁵Fe-vs-⁹⁰Sr bench proposal in `research/TIMING_PSD_NOTES.md`)*:
   a track's charge arrives spread over the full drift-time span
   (~60–110 ns depending on gas) while a gas conversion is a point deposit
   arriving at one time. Scan peaking time 25–200 ns in the T1 emulator
   (short t_p suppresses track PDO but not point-deposit PDO — it *widens*
   the charge-window separation, at an ENC cost that the noise data must
   adjudicate), and use the clustering of real-track TDOs to tighten the
   3-layer coincidence window below the full drift span (accidentals fall
   quadratically with the window).
   **Price side now measured** *(2026-08-06,
   `research/TIME_RESOLUTION_NOTES.md` §3)*: at a fixed 2 fC threshold the
   ballistic deficit costs MIP efficiency 99.4 % → 96.7 % → **77.8 %** going
   t_p = 100 → 50 → 25 ns in Ar/CO₂/Iso, and 92.6 % → 78.1 % → **42.1 %** in
   Ne/CH₄. Any short-peaking scheme must either lower the threshold in step
   or pay that ε_e directly — fold this into the handle-1 ROC before
   proposing it.

7. **Time resolution per pad hit** *(new step, 2026-08-06)*. Deliverable:
   σ_t(gas, gap, t_p, threshold, ENC) with the contribution budget, both raw
   and time-walk-corrected, plus the predicted (PDO, TDO) correlation — the
   quantity to compare against the SPS muon data. First pass is already done
   with a standalone toy (`vmm/time_resolution.py`, no Geant4/Stage B
   needed): **σ_t ≈ 10–13 ns in argon mixtures, 15–25 ns in neon mixtures**,
   dominated by primary-ionization statistics as σ_t ≈ (1.0–1.5)/(n_p·v_d).
   Re-run through the full chain once Stage B exists and once Magboltz v_d
   replaces the placeholder drift velocities (P0.10) — σ_t ∝ 1/v_d, so that
   is the leading input uncertainty. Feeds step 5 (coincidence window) and
   `TESTBEAM_PLAN.md` §2.6.

Validation inside Phase 2: a handful of full Garfield++ AvalancheMC/
microscopic events vs the sampler; sensitivity of ROC to Polya θ, ε_mesh,
and gain ±50 %. **Data anchoring**: once the SPS muon comparison
(`TESTBEAM_PLAN.md`) produces `tb_tuned_params.json` (gain anchor, θ,
ε_mesh, measured ENC), re-derive the Phase-2 headline numbers with tuned
parameters and quote tuned-vs-untuned as the modeling systematic.

### 5a. Source-validation predictions — bench anchor for the photon side

The muon test beam anchors only the electron/MIP side; the photon side
gets its bench anchor from lab sources, which map directly onto the
background band: **⁵⁵Fe 5.9 keV** (gain anchor *and* the canonical
point-deposit reference — our preferred source), **²⁴¹Am 59.5 keV**,
**¹⁰⁹Cd 88 keV**, **⁵⁷Co 122 keV**, plus **⁹⁰Sr/⁹⁰Y β** as the track-like
reference (spectrum CSV already in `sr90_calibration/`). Simulate each
source on the prototype geometry (P0.11 source mode), per gas: predicted
pad-charge/PDO spectra, fake-hit probability above the working charge
window, escape-peak and fluorescence structure. Deliverable: one
predicted-vs-measured figure per source; generate the predictions *before*
the bench data lands so the comparison stays blind-ish. The same source
set doubles as the Ne-mixture gain scan of §7.2 once the mixer exists.

## 6. Phase 3 — drift gap × gas matrix

The decision matrix. Axes:

- **Drift gap**: 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 mm (user-specified grid).
- **Gas**: the §7.4 list — A1–A3, N1–N8 + N10 (all ~13 through Phase 1 +
  Magboltz; down-select to ~8 for the full matrix; N9 ALICE-mix contingent
  on its drift time).
- **Kinematics subset**: e⁻ {30, 119, 155} MeV × {10°} (+ {0°, 30°} at
  119 MeV); γ {50, 80, 100, 150} keV × {10°}.
- Geant4 runs: 7 gaps × ~8 gases × ~9 kinematic points ≈ **500 runs**
  (10⁴–10⁵ e⁻ / 10⁶–10⁷ γ each) — routine HTCondor scale given the thin
  detector. Stage B/C re-runs are cheap (Python over cluster files).

Deliverables: the (gap, gas) heatmaps of — ε_e at fixed γ rejection; P₀
floor; pad multiplicity; drift-time span; **per-hit time resolution σ_t**
(2026-08-06 — cheap to add, and it feeds the coincidence window of §5 step 5
where it competes with the gas's photon-rejection difference); required mesh
voltage for reference gain. Plus a recommendation memo: operating point(s) to carry into
prototype tests, with the flammability/premix constraints from §7 applied.

## 7. Gas candidates — physics, flammability, procurement

*(Research pass done 2026-08-05; full citations in
`research/GAS_NOTES.md`. Physics numbers already merged into §1.)*

### 7.1 Flammability — the binding constraint (ISO 10156:2017)

Non-flammable in air *in any proportion* iff Σ(A′ᵢ/Tcᵢ) ≤ 1 after converting
inerts to N₂-equivalents (Kk: **Ar 0.55, Ne 0.7**, CO₂ 1.5, N₂ 1; Tci:
**iC₄H₁₀ 3.4 %, C₂H₆ 4.5 %, CH₄ 8.7 %**). Computed limits for us:

- binary **Ar/iC₄H₁₀: ≤ ~1.9 %** iso; **Ne/iC₄H₁₀: ≤ ~2.4 %** (Ne is a
  slightly better inertant than Ar); Ne/C₂H₆ ≤ ~3.2 %; Ne/CH₄ ≤ ~6.2 %.
- **All our ≥5 % isobutane candidates are flammable**, full stop — usable
  (COMPASS ran flammable Ne/C₂H₆/CF₄ for 20 years) but they buy sniffers,
  zoning and interlocks at MESA.
- CO₂ rescues ternaries: NSW Ar/CO₂/iC₄H₁₀ 93/5/2 passes at ratio 0.97
  (barely — the folkloric "<3 % iso is non-flammable" traces to the
  superseded 1996 edition); **Ne/CO₂/Iso 95/3/2 passes at 0.83, more margin
  than NSW itself.** ISO also allows certifying a marginal mix by spark-tube
  test. Governing CERN rule: Flammable Gas Safety Code G.

### 7.2 Operating experience with Ne in Micromegas — premise check

- Best direct comparison (microbulk, Iguaz et al. JINST 7 (2012) P04007):
  **max gain before sparking ~5× higher in Ne** mixtures, better energy
  resolution — but Ne needed slightly *higher* amplification field at equal
  gain. **The "Ne = same gain at lower voltage" lore is unverified for our
  150 µm bulk geometry**; the solid, verified Ne advantages are max-gain
  headroom and orders-of-magnitude lower discharge probability at fixed
  gain (Gasik et al.; ALICE GEM experience). Plan accordingly: don't sell
  the voltage argument until measured.
- **No published gain curve exists for Ne/iC₄H₁₀ or Ne/CO₂ in a 128–150 µm
  bulk MM** — a ⁵⁵Fe/⁹⁰Sr gain scan on the P2 prototype would be both
  necessary calibration input (§2 gain anchoring) and publishable.
- **Closest existing dataset: the PICOSEC eco-gas campaign (DRD1, Oct
  2025)** — 128 µm MM measured with Ne/Iso 95/5→75/25 *and* Ne/Iso/CO₂
  90/5/5, 85/5/10, 90/3/7: gains >10⁵, the 95/5, 90/10 and 90/5/5 curves
  nearly identical, CO₂ variants declared non-flammable. Strong validation
  that our candidate space is sane.
- Penning: Ne* (16.6–16.7 eV) Penning-ionizes *every* candidate quencher
  (Ar* cannot ionize CO₂/N₂) — r published only for Ne/CO₂ (0.577 @90:10);
  **Ne/iC₄H₁₀, Ne/CH₄, and all ternaries need assumed r, bracketed** (§2).

### 7.3 Procurement realities

- **No supplier stocks Ne detector premixes** (unlike Ar P10/ArCO₂); all of
  Air Liquide/Linde/Messer/Nippon do them as certified *custom* mixes,
  2–4 week lead, surcharge. Pure Ne 4.5/5.0 is stock everywhere.
- **Isobutane vapor pressure (3.1 bar abs @21 °C) kills full-pressure
  premixes above ~1.5 % iso** — a Ne/Iso 95/5 bottle can only be filled at
  reduced pressure (few liters of gas per bottle). Any serious isobutane
  option ⇒ **on-site mixing**. CO₂/CH₄/C₂H₆/CF₄ premixes are fine at
  200 bar.
- **Mixer**: 2–3 channel MFC racks are completely standard (CERN EP-DT
  mixes all LHC gases on-site; RD51 lab likewise); Bronkhorst/Vögtlin/MKS
  class accuracy holds the quencher fraction to ±0.05–0.1 points —
  a few-% gain effect, same order as P/T drift. **Recommendation: get a
  3-channel rack (~2–4 k€/channel)** — it costs about one MFC more than a
  2-channel and unlocks every ternary candidate plus in-situ quencher
  scans, which we want anyway given the missing published gain curves.
- Neon cost: **~50–100× argon per m³** in normal times; 2022 Ukraine-war
  shortage (~4× spot spike) has normalized, but budget via quotes. At
  single-wedge flush rates (~l/h) this is noticeable, not prohibitive;
  favors modest flow or recirculation at experiment scale.

### 7.4 Campaign gas list (amended)

| # | mixture | flammable? | premix @200 bar? | why |
|---|---|---|---|---|
| A1 | Ar/Iso 95/5 | **yes** | no (iso) | current default, reference |
| A2 | Ar/CO₂/Iso 93/5/2 † | no (0.97, barely) | no (iso; mix on-site) | NSW heritage reference |
| A3 | Ar/CF₄/Iso 88/10/2 | no | no (iso) | timing reference |
| N1 | Ne/Iso 95/5 | **yes** | no | the original proposal |
| N2 | Ne/Iso 90/10 | **yes** | no | quencher-recovery point |
| N3 | Ne/Iso 80/20 | **yes** | no | full Ar-equivalent n_p |
| N4 | Ne/CF₄ 90/10 | no | yes | fast, if timing decisive |
| N5 | Ne/CO₂/Iso 95/3/2 | **no (0.83)** | no (iso; on-site) | NSW-analog, best-margin ternary |
| N6 | Ne/Iso/CO₂ 90/5/5 | no (per PICOSEC) | no (iso; on-site) | PICOSEC eco-mix, measured gain >10⁵ |
| N7 | Ne/CH₄ 93/7 | no (calc.) | **yes** | simplest non-flammable binary, fast drift, NEWS-G heritage |
| N8 | Ne/C₂H₆ 90/10 | **yes** (limit ~3.2 %) | yes (ethane doesn't condense) | COMPASS chemistry sans CF₄; PICOSEC has it queued, no data yet |
| (N9) | Ne/CO₂/N₂ 90/10/5 | no | on-site | ALICE fallback; slow unsaturated drift (~2.8 cm/µs) — include in Magboltz pass, sim only if drift time acceptable |
| N10 | Ne/CO₂ 90/10 | no | **yes** | *added 2026-08-05*: the only Ne mix with a **published** Penning r (0.577) — anchors the r-bracketing for all other Ne mixes; NA49/CERES heritage; non-flammable and 200-bar premixable (simplest procurement of any Ne candidate). Slow drift (~110 ns over 3 mm) is acceptable at these gaps and *helps* the time-structure discriminant (§5 step 6) |

† amended from 95/3/2 to the NSW-standard 93/5/2 — the extra CO₂ is what
makes the ISO arithmetic pass; keep 95/3/2 as a Magboltz-only variant if
desired.

**Argon row promoted 2026-08-05 (late).** The A1–A3 entries were written as
"reference" mixtures — the incumbents the neon candidates were expected to
beat. Per §0/§1.1 they are now **co-baselines on equal footing**: at
50–100 keV the argon mixtures lose almost nothing on photon rejection
(the wall floor dominates and is gas-independent) while keeping a 30×
better zero-cluster floor and a decade of operational heritage. Two
consequences for the run matrix: **A2 (Ar/CO₂/Iso 93/5/2) is a serious
operating candidate, not just an NSW-heritage reference** — it is
non-flammable and has published Penning data — and a **fourth argon point,
Ar/CO₂ 93/7 (A4), should be added** to mirror N10: published r = 0.405,
200-bar premixable, no isobutane, and the closest thing to a
zero-unknowns baseline the campaign has. Argon's one real penalty in this
band is the Cu-K 8.05 keV re-absorption chain (§1.1), which is 10× larger
than in neon but is also the most easily rejected fake class.

**Down-select logic**: all ~13 (14 with A4) go through Magboltz tables +
Phase 1 Geant4 (cheap); the Phase 3 full matrix keeps ~8, **now with at
least three argon points, not one**; the final recommendation should name
(i) a non-flammable operating mixture, (ii) a flammable
maximum-performance variant if the gain justifies the infrastructure, and
(iii) the premix-vs-mixer procurement route for each — and must state
explicitly what the neon option buys once the wall floor is included,
since on present analytic numbers that is "max-gain headroom and discharge
probability, and essentially no photon rejection".

**Simulation caveats per gas** (carry through all result plots): assumed
Penning r for Ne/Iso, Ne/CH₄, ternaries; no measured 150 µm bulk gain
curves for any Ne mix (until prototype scan); CF₄ W-value spread (35–52 eV)
is a genuine systematic — quote both ends.

## 8. Phase 4 — optional refinements (only if earlier phases demand)

- ~~VMM3a waveform-level simulation with a **realistic induced current**~~
  **Moved out of Phase 4 on 2026-08-07 → P0.17.** MX17's review
  (`HANDOFF_MX17_RESPONSE.md` §2.1) pointed out the cost estimate was wrong
  for *our* stack: the expensive machinery is only needed for a **resistive**
  detector, and P2's is not, so the weighting potential is static and
  closed-form (Riegler / Garfield++ `ComponentParallelPlate`). An afternoon,
  not a phase — and it is the main model risk in the σ_t prediction we plan
  to compare against SPS, so it should not be sitting in "optional".
  *(µTPC is **not** on this list either way: it needs a multi-channel cluster
  and we have ~1.1 pads.)* Still Phase 4 from the original bullet: streaming
  pile-up and baseline wander.
- Mesh transparency from field maps / measured curves instead of a constant.
  *(2026-08-07: likely free — MX17 is solving a woven-mesh unit cell in
  Garfield++ neBEM, a desktop-scale job, for transparency ε(E_d/E_a),
  funneling and ion endpoints. P2's mesh differs only in weave parameters,
  so the same script with different constants covers this item and the
  measured-transparency escalation path in `TESTBEAM_PLAN.md` §6. MX17 will
  link the script path in `HANDOFF_MX17_RESPONSE.md` §2.3 when it exists —
  **ask before writing our own.** Related: P0.13's woven-mesh cross-check.)*
- Pillar geometry refinement in Geant4. *(Woven-mesh geometry moved out of
  Phase 4 on 2026-08-05 — it is now P0.13 and blocking, since the mesh is
  the single largest fake source in the analytic budget.)*
- Degrad spot-checks for 3–10 keV gas photoabsorption cascades.
- Beam-spot/angular-spread and radiative-tail-weighted source spectra.
- Overpressure/sag and mapping-revision sensitivity checks.
- **Full 3-wheel geometry with realistic coincidence counting** (wanted
  eventually — Dylan 2026-08-05; until then the P0.4c exit-particle
  convolution is the interim estimate of correlated fakes).
- **Photon-shield trade study**: thin high-Z front shield vs the
  multiple-scattering cost against the ±1° tracking requirement (~2° for
  0.5 mm Pb at 119 MeV; CEA 2024 internship precedent). *Sidelined
  2026-08-05 (baseline experiment plan mostly fixed) — revisit if the
  conversion budget or the correlated-fake rate demands it.*

## 9. Bookkeeping

- **Run manifest**: one CSV/JSON row per run point (runID, mode, gas, gap,
  particle, E, θ, position, N, seed, code git-hash, geometry-parameter hash)
  — written by the submit script, joined by analysis. No un-manifested runs.
- **Naming**: `p2_<gas>_<gap>mm_<part><E><unit>_th<deg>_r<mm>` (+ seed).
- **Run manifest MUST carry the thrown-event count.** With `--skip-empty`
  the trees hold only events with a gas cluster, so every per-incident-photon
  probability is `entries / thrown`, and `thrown` exists only in the run
  summary printout. Any analysis that divides by `GetEntries()` is silently
  wrong by ~300×. (`scripts/collect_results.py` now carries a warning.)
- **Output schema**: `docs/OUTPUT_FORMAT.md` — EventTree (incl.
  `interactionClass` and the per-layer sums), ClusterTree (provenance),
  VolumeTree (every volume), plus the caveats on `nPrimary` and `edepAmp`.
- **Output volume**: measured, not estimated — see `docs/BENCHMARKS.md`.
  γ runs with `--skip-empty` are ~6 bytes/thrown photon at 60 keV; e⁻ runs
  ~1 kB/event. Whole campaign comfortably < 50 GB; keep on EOS, analysis
  pulls summaries.
- **Statistics floor**: size γ samples so the *rarest* reported quantity
  (P(fake hit) after threshold, could be 10⁻⁵) has ≥100 events → up to 10⁷
  thrown γ at the critical energies; cheap because empty events are skipped.
- Keep `P2Wedge.hh`/`SimConfig.hh` ↔ `p2_model.py` ↔ docs in sync on any
  geometry change (standing rule).

## 10. Risks / provisional items

1. **Geometry provisional until Alexandra's review** — do not launch
   Phases 1–3 production; Phase 0 is safe to do now.
2. **Gain model absolute scale** — no published 150 µm bulk gain curves for
   any Ne mixture, and no published Penning r for Ne/Iso, Ne/CH₄ or any
   ternary. Mitigated by reporting vs gain/threshold in primary-equivalents,
   bracketing r, and anchoring to prototype measurements (a ⁵⁵Fe/⁹⁰Sr gain
   scan in Ne mixes is needed calibration *and* publishable). The "Ne runs
   at lower voltage" premise is unverified — the defensible claims are
   ~5× max-gain headroom and far lower discharge probability.
   **New mitigation (2026-08): the SPS 200 GeV muon test beam** — its
   as-run gas gets a data-anchored gain scale, Polya θ, ε_mesh and ENC via
   `TESTBEAM_PLAN.md`; transfer to Ne mixes still rides the Penning-r
   assumption until the Ne source scan exists.
3. **Solid-slab mesh — checked 2026-08-05, and it is sound.** The mesh *is*
   the largest single fake source (§1.1), but the effective-density slab
   reproduces the plain-weave areal density **exactly**
   (π·d²/(2·pitch) = 8.46 µm solid-equivalent), so the conversion rate is
   right — and the simulated value landed on the analytic one. Only the
   escape geometry is approximate (~1.4× between diluted-slab and
   solid-wire estimates). P0.13 is a cross-check, not a blocker. Pillars
   remain a genuine second-order caveat. *(An earlier version of this entry
   and of the photon note claimed the slab over-counted interactions; that
   was an arithmetic error in the note, not a bug in the code.)*
4. **The wall-floor conclusion itself is analytic, not simulated.** §0/§1.1
   rest on a uniform-depth, straight-line-range escape model good to a
   factor 2–3. It would take a factor ~20 error to restore the neon
   argument, which is why the conclusion is stated at all — but no gas
   decision should be *finalized* before Phase 1 confirms it with real
   geometry. Falsification criteria are listed in
   `research/PHOTON_DISCRIMINATION_NOTES.md` §8.
5. **We do not know the spectrum *or the spatial/angular distribution* of
   either species** — signal electrons or background photons. Full
   treatment in **`NEEDED_INPUTS.md`** (raised by Dylan 2026-08-05). Short
   version: the "15–20 kHz/pad" and "10³ γ/e⁻" figures are area averages,
   and accidental coincidence rates go as the *square or cube* of the local
   hit probability, so a factor ~5 radial non-uniformity is a factor 25–125
   in accidental triples — larger than the entire gas effect. Incidence
   angle is also correlated with radius rather than being the free global
   parameter the plan treats it as, which makes the "~1.1 pads/track"
   assumption behind demoting topology radius-dependent. ⬜ The single
   highest-value ask is a **phase-space file at the wheel plane** for both
   species (P0.15 consumes it). Until then run the §4.2 energy grid and
   present results per energy, not folded.
6. **The two pad-mapping revisions disagree on the readout order, not the
   geometry** *(restated 2026-08-06 — the earlier "79/1280 pad positions
   disagree" was a rounding artifact in `analyze_p2_readout.py`, now fixed;
   matched by nearest neighbor the two pad planes agree to 2.5 µm, i.e.
   print precision)*. The real discrepancy is larger and different: only
   **11/1280 (connector, channel) pairs land on the same pad**. Geometric
   pad-level observables (multiplicity, cluster size, position resolution)
   are therefore *safe*; everything channel-level is not — VMM neighbor
   logic (chan ±1 is a 12 mm neighbor in one revision and a 126 mm one in
   the other, `vmm/nl_map.py`), dead/noisy-channel masks, and any
   per-channel comparison with test-beam data. Ask which revision the DAQ
   uses (`testbeam/TB_CONDITIONS.md` T0.4).
7. **W-value/PAI approximations** — validated internally (P0.6) and by the
   Ar-vs-Ne *ratio* being robust even if absolutes shift. Note the Penning
   caveat in §5 step 3a: a single tabulated W per gas mis-scales Ne against
   Ar, and the current `floor(edep/W)` treatment gets the *width* wrong.
8. **Space-charge gain suppression is unmodelled and unmeasured** for our
   geometry, and it acts precisely on the large deposits the discrimination
   depends on (§5 step 3a). Any ROC conclusion that leans on the far upper
   tail of the charge window should be flagged as gain-model-limited.
9. **Air world** — fine for this study; revisit only for full-wheel work.
