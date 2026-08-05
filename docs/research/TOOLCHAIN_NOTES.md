# Research notes — Geant4 ↔ Garfield++ ↔ VMM simulation chain

**Compiled 2026-08-05** (web research pass for the campaign plan; summarized
into `../SIM_CAMPAIGN_PLAN.md` §2 — this file keeps the full detail and
citations). Unverified items flagged inline.

---

## 1. Geant4 ↔ Garfield++ coupling

**Reference paper:** D. Pfeiffer, L. De Keukeleere, C. Azevedo, F. Belloni,
S. Biagi, V. Grichine, L. Hayen, A. R. Hanu, I. Hřivnáčová, V. Ivanchenko,
V. Krylov, H. Schindler, R. Veenhof, *"Interfacing Geant4, Garfield++ and
Degrad for the simulation of gaseous detectors"*, NIM A 935 (2019) 121-134,
[arXiv:1806.05880](https://arxiv.org/abs/1806.05880). Division of labor:
Geant4 does primary generation and all interactions in non-gas material;
Garfield++ (Heed/SRIM/Degrad for primary ionization, Magboltz-based
transport) does electron/ion drift, avalanche, induced signals.

**Two coupling patterns:**

**(a) Embedded "Garfield model"** (fast-simulation parameterization): a
`G4Region` over the gas volume + a `G4VFastSimulationModel` subclass
(`IsApplicable()`, `ModelTrigger()`, `DoIt()`) kills the Geant4 track in the
region and hands position/direction/energy to Garfield++
(`TrackHeed::NewTrack()`, then drift/avalanche). A `GarfieldPhysics`
singleton maps particle + KE window → Heed or Geant4. Shipped as
`Examples/Geant4GarfieldInterface`; tutorial:
<https://garfieldpp.web.cern.ch/tutorials/g4/>; extended Degrad version:
<https://github.com/lennertdekeukeleere/Geant4GarfieldDegradInterface>.
Pros: one executable, consistent geometry, seamless handoff of secondaries
born in gas. Cons: monolithic job throttled by the avalanche stage;
ROOT+Geant4+Garfield co-linked build; can't re-run only detector response
when scanning HV/gas; unit conversions (mm/MeV vs cm/eV) fiddly. The
official tutorial itself warns it will "hardly be fast".

**(b) Staged pipeline** (export Geant4 clusters/segments; standalone
digitizer downstream). This is what large experiments do:

- ATLAS NSW Micromegas digitization: Geant4 deposits → *parameterizations of
  Garfield/Magboltz* (drift velocity, D_L/D_T, Lorentz angle) → Polya gain →
  VMM transfer function.
  [Springer chapter "Simulation of the ATLAS New Small Wheel System"](https://link.springer.com/chapter/10.1007/978-981-13-1316-5_25).
- CAST/IAXO Micromegas background chains:
  [arXiv:1310.3391](https://arxiv.org/pdf/1310.3391),
  [BabyIAXO arXiv:2509.08138](https://arxiv.org/html/2509.08138);
  institutionalized in REST-for-Physics.

**Recommendation for the P2 campaign: (b).** Geant4 pass runs once per
beam/background config; cheap digitization re-runs per (field, gain,
threshold, gas) point from the same files. Garfield++ available on
lxplus/Condor via cvmfs LCG views (`/cvmfs/sft.cern.ch/lcg/views/...` then
`share/Garfield/setupGarfield.sh`; [install docs](https://garfieldpp.web.cern.ch/install/)).
Use (a) selectively for one-off validation (e.g. x-ray response with in-gas
Geant4→Heed handoff).

## 2. Primary ionization fidelity in 1-4 mm gas

Geant4 default condensed-history EM + edep/W is **not** adequate for cluster
statistics in a thin gap: right *mean*, wrong event-by-event distribution
(Urban fluctuation model; classic thin-absorber problem — J. Apostolakis et
al., NIM A 453 (2000) 597). The zero/low-cluster tail drives threshold
efficiency.

Fidelity ladder:
1. **G4PAIModel / G4PAIPhotModel** per-region: correct thin-layer straggling
   and delta spectrum; benchmarked vs Heed in the Pfeiffer paper. Timing
   there (10⁴ × 1 GeV e⁻ in 1 cm Ar/CO₂): Heed 5.8 s, combined Geant4/Heed
   14.7 s, Geant4-PAI 17.9-34.7 s — all trivial next to avalanching.
2. **Heed** (I. Smirnov, NIM A 554 (2005) 474) via `TrackHeed`: the MPGD
   standard — correct cluster count, cluster-size distribution (incl. rare
   large clusters), per-electron positions.

**Recommended practice (Pfeiffer et al., for relativistic particles): Heed
does charged-particle ionization inside the gas; Geant4 everything else.**
In the staged pipeline: export the gas-crossing track segment (entry,
direction, momentum, PID) and re-run `TrackHeed::NewTrack()` in the
digitizer; or minimally validate edep/W clusters against Heed and keep a
Heed mode for threshold studies. If mixing Geant4 + Heed ionization in one
volume: Geant4 production cut in the gas must sit between the lowest
ionization potential and W (empirically ~21-24 eV for Ar/CO₂ 70/30, ~19.5 eV
for He/iC₄H₁₀ 70/30 — Pfeiffer paper), coordinated with the low-energy
electron limit, else W and Fano come out wrong. keV photoelectrons are below
Heed's PAI validity → handled via its delta-electron transport; the paper's
combined Geant4-PAI + Heed scheme uses a 0.1-2 keV transfer threshold.

Shortcut bookkeeping: nPrimary = edep/W fine for the mean; add Fano smearing
(F ≈ 0.2-0.3 for Ar mixtures) — still misses cluster clumping, which matters
much more at 1 mm than 4 mm.

## 3. X-ray response (3-10 keV fluorescence; 50-150 keV incident)

**Degrad** (S. Biagi, with [Magboltz](https://cern.ch/magboltz)): photons to
2 MeV, electrons to 4 MeV, shell-resolved photoabsorption, Auger +
Coster-Kronig, **shake-off**, fluorescence; noble-gas Fano to <3 % below
20 keV. Needed only when the exact electron-count distribution of a gas
photoabsorption must be right (⁵⁵Fe/Cu-K resolution, escape-peak topology).
Coupled to Geant4 by file I/O in the Pfeiffer paper.

Lighter standard: Geant4 photoelectric with deexcitation on
(`/process/em/fluo`, `auger`, `pixe`; since 10.2), or Heed
`TrackHeed::TransportPhoton()` (returns conversion-electron list — the usual
Garfield++ x-ray route). Gets mean yield and escape/fluorescence branchings
right; lacks shake-off, less exact Fano.

**Hierarchy for 50-150 keV incident photons: surroundings dominate the gas.**
At 60 keV, μ/ρ(photo, Ar) ~10⁻² cm²/g → direct gas conversion in 3 mm NTP
~10⁻⁵-10⁻⁴, falling to 150 keV where Compton dominates (numbers from memory
— **recompute with NIST XCOM per gas**). Detected background dominated by:
- photo-/Compton electrons ejected from Cu pads, steel mesh, cathode into
  the gas (short, high-dE/dx stubs, energy up to E_γ);
- secondary fluorescence: **Cu Kα 8.05 keV (anode), Fe Kα 6.4 keV / Cr Kα
  5.4 keV (mesh), Ar escape at E−3.2 keV** — exactly what CAST Micromegas
  background sims observed ([arXiv:1310.3391](https://arxiv.org/pdf/1310.3391)).
  These few-keV photons photoabsorb efficiently in Ar and deposit fully —
  the population the ADC cut must handle.

Discrimination observable: 6-8 keV photoabsorption = point-like ~230-300 e⁻
on 1-2 pads; MIP over 3 mm = ~30 e⁻ spread along the track → separation is
topological as much as total charge; diffusion + pad sharing must be modeled.

## 4. Magboltz gas tables & Penning transfer

Standard flow:

```cpp
MediumMagboltz gas("Ar", 93., "CO2", 7.);
gas.SetPressure(760.); gas.SetTemperature(293.15);
gas.SetFieldGrid(100., 100.e3, 20, true);   // log-spaced E; B, angle optional
gas.GenerateGasTable(10);                    // ncoll × 1e7 collisions/point
gas.WriteGasFile("ar93co2_7.gas");
gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
```

Tabulated per (E,B,θ): v_drift, D_L/D_T, Townsend α, attachment η, and
per-level excitation/ionization rates → Penning rescaling **without
regenerating the table**. `MergeGasFile()` combines tables.
[Gas files tutorial](https://garfieldpp.docs.cern.ch/tutorials/gasfiles/).

CPU: ncoll = 2-5 → ~1 % for quenched mixtures; ncoll ≥ 10 for <0.5 % and
pure nobles; error ∝ 1/√ncoll, CPU ∝ ncoll·nE·nB·nθ. No official timing
published (flag). Datapoint: pure Xe @10 bar, 19 E-points, ncoll=20 → 1-2
days ([ROOT forum](https://root-forum.cern.ch/t/time-to-generate-gas-file-using-magboltz/63729));
quenched Ar/Ne at 1 atm ≈ minutes/E-point, ~hour per 20-point grid/core.

**Penning:** `EnablePenningTransfer(r, lambda[, species])`; argument-free
variant (since late 2021) applies literature parameterizations — **binary
mixtures only**. Each excitation above ε_ion ionizes with probability r;
Townsend rescaled; AvalancheMicroscopic converts stochastically; λ = 0
normally. Without it, simulated gains low by large factors (e.g.
[arXiv:1606.04852](https://arxiv.org/abs/1606.04852) needed r = 0.57 for
Ar/30%CO₂).

Published r (Şahin/Kowalski/Veenhof):
- [JINST 5 (2010) P05002](https://iopscience.iop.org/article/10.1088/1748-0221/5/05/P05002),
  1 atm, 10 % quencher: **Ar/CH₄ 0.212, Ar/C₂H₆ 0.31, Ar/C₃H₈ 0.43,
  Ar/iC₄H₁₀ 0.40**, Ar/C₂H₂ ~0.72; Ar/CO₂ r(c) = (a₁c+a₃)/(c+a₂), a₁≈0.62.
- NIM A 768 (2014) 104 (Ar-CO₂ precision); JINST 12 (2017) C01035
  (pressure-dependent Ar/CO₂); JINST 11 (2016) P01003 (Ne-CO₂ gains);
  JINST 16 (2021) P03026 (Ne-CO₂, Ne-N₂ model → Garfield++ built-ins).
- Built-in values @1 atm: **Ar/CO₂ 93:7 → 0.405; 90:10 → 0.456; 80:20 →
  0.539; 70:30 → 0.574; Ne/CO₂ 90:10 → 0.577; Ne/N₂ 95:5 → 0.537**;
  Ar/CO₂ 93:7 → ~0.34 @3 atm
  ([RD51 Nov 2021, Alsamak/Schindler/Veenhof/Şahin](https://indico.cern.ch/event/1071632/contributions/4602415)).
- **Gaps (flagged):** no published r for Ar/iC₄H₁₀ @5 % (Garfield++ reuses
  0.40 with a warning); **Ar/CF₄ r ≈ 0** (CF₄ IP 15.9 eV above Ar* levels);
  **Ne/CF₄ allowed energetically but no published r; ternaries none** (RD51
  2021: "still being modeled"; common practice: apply binary value to
  noble-gas excitations manually). Example anchors: `Examples/Gem` r = 0.51
  (Ar/CO₂ 80:20); G4 interface example 0.57 (70:30);
  `Examples/ResistiveMicromegas` (Ar/CO₂ 93:7) uses the automatic
  parameterization.

## 5. Avalanche modeling, 150 µm gap

- **AvalancheMicroscopic**: gold standard (gain, Polya shape, mesh
  transparency, timing); cost ∝ avalanche size
  ([CUDA paper arXiv:2509.15377](https://arxiv.org/abs/2509.15377): gain 10⁴
  already impractical for CPU campaigns; GPU 60-100× at 10⁶-10⁷). O(s-10s)
  per electron at gain 10⁴ (flag: no authoritative published number).
- **AvalancheMC**: 10-100× faster; gain number OK; no mesh-funneling
  microscopics.
- **Parameterized (campaign workhorse)**: Polya P(g) ∝ (g/ḡ)^θ e^(−(1+θ)g/ḡ),
  θ ≈ 1.5-2.5 for Ar-CO₂ MM (Zerguerras et al., NIM A 772 (2015) 76,
  [CDS PDF](https://cds.cern.ch/record/2153708/files/nima772-76.pdf));
  ḡ from measured curves or Penning-corrected Townsend integral; × mesh
  transparency T(E_drift/E_amp), flat ≈100 % for ratio ≲0.01, steep drop
  above (Kuger/NSW studies; exchangeable-mesh paper
  [NIM A 2015](https://www.sciencedirect.com/science/article/abs/pii/S016890021501373X)).
  Standard pattern: NSW digitization, PICOSEC phenomenological model
  ([arXiv:2010.13535](https://arxiv.org/pdf/2010.13535)), CAST/IAXO, fast GEM
  digitization ([arXiv:1904.06142](https://arxiv.org/pdf/1904.06142)).

Tiering: (i) AvalancheMicroscopic once per (gas, HV) on 10³-10⁴ single
electrons in a unit-cell field map (ANSYS/Elmer/neBEM) → extract ḡ, θ,
transparency, current profile; (ii) freeze into fast digitizer; (iii) anchor
ḡ to calibration data. (P2 is non-resistive bulk MM, so no resistive-layer
RC spreading needed; Garfield++ `Examples/ResistiveMicromegas` + Riegler's
extended Ramo-Shockley exist if that changes.)

## 6. VMM3a electronics

Chip (verified): 64 ch; CA + **3rd-order shaper, one real pole + complex
pair, "Delayed Dissipative Feedback"** (not plain CR-RC³); t_p
25/50/100/200 ns; gains 0.5/1/3/4.5/6/9/12/16 mV/fC; global 10-bit threshold
DAC + 5-bit/ch trim; peak-sensing 10-bit ADC (~8-bit effective, two-stage
domino) + 8-bit TAC + 6-bit fast ADC; ~250 ns conversion, ~4 MHz/ch;
**neighbor logic** (reads below-threshold neighbors, crosses chip
boundaries); VMM3a adds ion-tail suppression. Sources: G. Iakovidis et al.,
NIM A 955 (2020) 163306 ([OSTI full text](https://www.osti.gov/servlets/purl/1603751));
De Geronimo et al., J.Phys.Conf.Ser. 1498 (2020) 012051
([ATL-MUON-PROC-2019-009](https://cds.cern.ch/record/2693463)); shaper
concept: De Geronimo & Li, IEEE TNS 58 (2011) 2382; SRS: Scharenberg et al.,
NIM A 1031 (2022) 166548, [arXiv:2109.10287](https://arxiv.org/abs/2109.10287).

Published response models:
1. **ATLAS Athena `MM_Digitization/VMM_Shaper.cxx`**
   ([GitLab](https://gitlab.cern.ch/atlas/athena/-/blob/main/MuonSpectrometer/MuonDigitization/MM_Digitization/src/VMM_Shaper.cxx),
   constants credited to G. Iakovidis, thesis
   [CERN-THESIS-2014-148](https://cds.cern.ch/record/1955475) §7.1.3):
   h(t) = Q·A·[K₀ e^(−p₀t) + 2|K₁| e^(−Re p₁ t) cos(−Im p₁ t + arg K₁)],
   K₀ = 1.584, K₁ = −0.792 − 0.115i, p₀ = 1.263/a, p₁ = (1.149 ± 0.786i)/a,
   a = t_peak/1.5 [ns]. Includes peak detector, threshold-to-peak timing,
   ballistic-deficit scale factor vs the ~150 ns ion tail, neighbor-logic
   emulation (re-run neighbors with lowered threshold).
2. **ESS/CERN CR-RC³ proxy**: F(t) = g·e³·(t/t_p)³·e^(−3t/t_p), t_p ≈ 200 ns
   (Backis, Piscitelli, Pfeiffer et al.,
   [arXiv:2510.01981](https://arxiv.org/abs/2510.01981)) — adequate for
   amplitude-level observables.

Charge-integration-per-pad + threshold: defensible baseline (no published
A/B comparison — flag) with corrections: (i) ballistic-deficit factor if
t_p < ~150 ns; (ii) threshold in ENC units on shaped amplitude (noise vs t_p
and pad capacitance 75-170 pF); (iii) neighbor logic (else cluster
size/charge underestimated). Breaks for: per-hit timing, time-walk,
streaming pile-up → use Athena h(t). Watch **ADC saturation**: high gas gain
× 6-8 keV photoabsorption (~250 e⁻) can saturate ~8-bit range while MIP pads
stay linear — interacts with the discrimination cut.

## Unverified-items recap

Absolute CPU for Magboltz per point and AvalancheMicroscopic per avalanche
(only scaling laws + one forum datapoint); XCOM numbers from memory; Penning
r for Ne/CF₄, Ne/iC₄H₁₀, Ar/iC₄H₁₀ 95/5, all ternaries; no published
charge-integration vs full-shaper comparison; some ATLAS proceedings
verified via mirrors/Athena source rather than CDS PDFs.
