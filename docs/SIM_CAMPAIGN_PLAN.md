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
§0, §1, §3, §4.2, §5, §7 and §10 below are amended from it)**.

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
(iii) **neighbor logic** — without it recorded cluster size/charge is
underestimated, which directly biases the electron-vs-photon topology cut.
Also simulate explicitly: **ADC saturation** — at high gas gain a 6–8 keV
photoabsorption (~250 e⁻ point deposit) can saturate the ~8-bit range while
MIP pads stay in range; saturation interacts with the discrimination cut.
Full waveform simulation with the Athena h(t) is Phase 4, needed only for
µTPC timing, time-walk, or streaming pile-up studies.
**Implementation started (2026-08-05): `vmm/`** — Athena shaper ported to
Python and unit-tested, channel emulator with T0 (charge+threshold) and T1
(shaper + neighbor logic + ENC + ADC) tiers, Stage-B interface contract and
NEEDS-DATA list; see `vmm/README.md`.

## 3. Phase 0 — preparation (code + infrastructure)

Do now (none of it depends on Alexandra's answers):

| # | task | notes |
|---|---|---|
| P0.1 | **Fix gas mixture composition bug** ✅ **done 2026-08-05** (see `research/GAS_FIX_NOTES.md`) | `makeMix2/3` feeds volume fractions to `AddMaterial()` which takes **mass** fractions. Negligible-ish for Ar/Iso (5 % vol = 7.3 % mass) but fatal for Ne mixes (Ne/Iso 95/5 vol = 86.5/13.5 by mass). Rewrite builder to take vol% + component molar masses → mass fractions; recompute densities from ideal-gas mixing. Also verify the 2.67e-3 g/cm³ isobutane density (ideal-gas at 20 °C gives 2.42e-3). |
| P0.2 | **Add campaign gases + W values** | the §7.4 list: Ar/CO₂/Iso 93/5/2; Ne/Iso 90/10, 80/20; Ne/CO₂/Iso 95/3/2; Ne/Iso/CO₂ 90/5/5; Ne/CH₄ 93/7; Ne/C₂H₆ 90/10 (+ Ne/CO₂/N₂ 90/10/5 optional); recompute all W-values as harmonic vol-weighted mixes; document each. CF₄ W spread (35–52 eV) → parameterize, don't hardcode. |
| P0.3 | **Gun angle control** | `--gun-theta` (deg from wedge normal), `--gun-phi` (tilt azimuth), keeping the aim point on the wedge fixed; optional `--beam-spread` later. |
| P0.4 | **Cluster provenance in output** ✅ **done 2026-08-05** | per cluster: creator process + logical volume where the depositing track (or its ionizing ancestor) was born. This is what turns photon runs into a per-layer conversion budget. Add per-event "first interaction volume" for γ runs. **Needs a `TrackingAction` + `G4VUserTrackInformation`, neither of which exists in the repo yet** (`ActionInitialization` registers only primary/run/event/stepping) — one-level `parentID` is not enough, since a δ-ray of a photoelectron reads as `eIoni`. Carry {origin process, origin volume, origin energy, ancestor trackID} and copy it to every secondary at creation. |
| P0.4b | **Gas-crossing track-segment export** | per charged track crossing DriftGas: entry point, direction, momentum, PID (a small `SegmentTree`). This is the Heed handoff for Stage B's high-fidelity mode (re-ionize with `TrackHeed::NewTrack()`), per the Pfeiffer-paper recommended split. Cheap to add now, enables the threshold studies later. |
| P0.4c | **Exit-particle export** | per event: every particle leaving the wedge envelope (PID, energy, direction, exit surface, creator process + birth volume, ancestry link to the primary). In photon runs this bounds **correlated multi-wheel fakes** (a Compton-scattered 100–150 keV photon keeps most of its energy and can convert again in wheels 2/3 — accidental-rate arithmetic misses this entirely) by pure wheel-to-wheel convolution, no 3-wheel geometry needed yet; it is also the input the eventual full 3-wheel sim (§8) will want. Decision 2026-08-05: save full interaction histories — track everything that produces signal, where it came from, and how it interacted. |
| P0.5 | **`--skip-empty` output flag** ✅ **done 2026-08-05** | photon runs are ≥99.9 % empty events; write only events with ≥1 cluster (keep the total-thrown count in metadata for normalization). |
| P0.6 | **PAI model in gas regions + deexcitation** 🟡 **deexcitation + fine-cut region done 2026-08-05; PAI still open** | G4Region over DriftGas+AmpGas with PAI/PAIphot (~2× EM cost, thin-layer straggling fixed); switch physics list to Livermore-class EM with fluorescence/Auger/PIXE on (needed for Cu-K/Fe-K fluorescence in photon runs). CLI switch to compare with default EM (one validation study early in Phase 1). **⚠ Includes a bug fix, do not skip:** `PhysicsList::SetCuts()` sets a **0.1 mm γ production cut**, which in copper is ~5–10 keV — i.e. the **8.05 keV Cu-K fluorescence line is at or below threshold and is probably being suppressed today.** That is the dominant *gas-sensitive* wall channel (§1.1). Set `G4EmParameters::SetDeexcitationIgnoreCut(true)` and add a fine-cut (~1 µm) `G4Region` over mesh + pad Cu + gas; the current 10 µm e⁻ cut is ~1 keV in gas but ~80 keV in Cu. |
| P0.7 | **Condor scripts refresh** | `submit_condor*.py` still speak MX17 modes; teach them the p2 flags + campaign manifest (one JSON/CSV row per run point). |
| P0.8 | **First lxplus build + geometry validation** ✅ **done 2026-08-05** | compile, overlap check, 100-event smoke run, one event display. Can be done with provisional geometry — worth doing *before* Alexandra's answers so her changes land in verified code. |
| P0.9 | **Stage B skeleton** | pad-map loader exists (`analyze_p2_readout.py`); write the drift/gain/pad-summing sampler with pluggable gas tables; unit-test with fake tables so Magboltz becomes a drop-in. |
| P0.10 | **Magboltz gas tables** | generate on lxplus per gas (log E-grid spanning drift ~100 V/cm–1 kV/cm and amplification ~30–60 kV/cm; ncoll = 10), `WriteGasFile` + ion mobility files, cache in repo or EOS. Order an hour per gas per core; embarrassingly parallel. Garfield++ comes from cvmfs LCG views on lxplus. Enable Penning transfer per §2 (manual r for Ne/iC₄H₁₀ and ternaries — document the chosen value and bracket). |
| P0.11 | **Source mode** | isotropic point/disc source gun with line energies — ⁵⁵Fe 5.9/6.5 keV, ²⁴¹Am 59.5 keV, ¹⁰⁹Cd 88 keV, ⁵⁷Co 122/136 keV — plus the ⁹⁰Sr/⁹⁰Y β spectrum (CSV already in `sr90_calibration/`); configurable standoff. Powers the §5a source-validation predictions and the ⁵⁵Fe-vs-⁹⁰Sr timing bench test (`research/TIMING_PSD_NOTES.md`). |
| P0.12 | **Per-event interaction classification + edep by region** ✅ **done 2026-08-05** | the branch that makes the requested money plot a one-liner. Per event: `interactionClass` ∈ {none, gas-photoelectric, gas-Compton, gas-Rayleigh, wall-photoelectric, wall-Compton, fluorescence-reabsorption, MIP-crossing, other}, `conversionVolume`, `conversionZ`, plus **edep in every region, not just DriftGas/AmpGas** — `FrontGas`, `DriftCathode_Gas` and `BackGas` are currently unscored, so "converted somewhere that produces no signal" is indistinguishable from "did not convert". Also record amp-gap cluster z *relative to the mesh* explicitly (Stage B needs it for partial gain, §2). Then "event-by-event edep separated by interaction type" is one `TTree::Draw` per gas. Builds on P0.4's TrackingAction; do them together. |
| P0.13 | **Cylindrical-wire mesh cross-check** *(downgraded 2026-08-05 — see note)* | the mesh **is** the largest single fake source (measured 1.98×10⁻³ per 60 keV photon, 9× the argon gas rate), but the existing effective-density slab is **not** wrong about the rate: its areal density π·d²/(2·pitch) = 8.46 µm solid-equivalent is *exactly* the plain-weave value (verified). Only the escape geometry is approximate, and the diluted-slab and solid-wire estimates differ by ~1.4×, not 3×. So this is a **cross-check worth doing, not a blocker**. |
| P0.14 | **Analytic cross-check gate** ✅ **done 2026-08-05** | `scripts/photon_budget.py` reproduces the §1.1 targets; the first Geant4 photon run matched it (3.4×10⁻³ measured vs 4.3×10⁻³ predicted total, mesh dead on). Re-run this comparison after the P0.1 gas-mixture fix, since a mass/volume-fraction error shows up here first. |
| P0.15 | **Phase-space gun** *(new 2026-08-05, from `NEEDED_INPUTS.md` §1)* | read primaries (position, direction, energy, PID, weight) from a file at the wheel plane instead of the fixed pencil beam, so results can be weighted by the real spatial/angular distribution of signal electrons and background photons. Cheap to write; **blocked on the collaboration supplying the file**, which is the single highest-value external ask we have. |

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
   leaves ~30–50 e⁻ *spread along the track* — so diffusion, pad sharing and
   neighbor logic must be modeled honestly, and both per-pad and per-event
   charge cuts should be scanned. High-energy wall photoelectrons ranging
   through the gas are the hard case (electron-like charge); flag their
   fraction separately. Expectation check (Dylan, 2026-08-05): at 10° a
   track averages only **~1.1 pad hits**, so topology is a *secondary*
   discriminant at the nominal angle — quantify what it adds, but the
   working assumption is charge deposition + 3-layer coincidence carry the
   weight.
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
   report both.
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
floor; pad multiplicity; drift-time span; required mesh voltage for
reference gain. Plus a recommendation memo: operating point(s) to carry into
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

- VMM3a waveform-level simulation: convolve induced currents with the ATLAS
  Athena `VMM_Shaper` closed-form h(t) (pole/residue form, ~20 lines) —
  needed only for µTPC timing, time-walk, streaming pile-up.
- Mesh transparency from field maps / measured curves instead of a constant.
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
6. **Two pad-mapping revisions disagree (79/1280)** — pad-level observables
   carry that caveat until the DAQ mapping is confirmed.
7. **W-value/PAI approximations** — validated internally (P0.6) and by the
   Ar-vs-Ne *ratio* being robust even if absolutes shift. Note the Penning
   caveat in §5 step 3a: a single tabulated W per gas mis-scales Ne against
   Ar, and the current `floor(edep/W)` treatment gets the *width* wrong.
8. **Space-charge gain suppression is unmodelled and unmeasured** for our
   geometry, and it acts precisely on the large deposits the discrimination
   depends on (§5 step 3a). Any ROC conclusion that leans on the far upper
   tail of the charge window should be flagged as gain-model-limited.
9. **Air world** — fine for this study; revisit only for full-wheel work.
