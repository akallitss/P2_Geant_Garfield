# P2 @ MESA and the BASKET backward-angle program — research notes

**Compiled 2026-08-05** from public sources (deep web dive; every claim linked).
Status legend: ✅ well-sourced · ⚠️ inferred/computed (method shown) · ❌ not
public. Purpose: physics context and beam-condition inputs for this simulation.

---

## 1. P2 experiment basics

| Quantity | Value | Status | Source |
|---|---|---|---|
| Physics goal | sin²θW from PV asymmetry in elastic e-p; ⟨A⟩ = −39.94 ppb, ΔA = 0.56 ppb → Δsin²θW/sin²θW ≈ 0.13–0.15% | ✅ | [arXiv:1802.04759](https://arxiv.org/abs/1802.04759) (EPJA 54 (2018) 208) |
| Beam energy | **155 MeV** | ✅ | arXiv:1802.04759 |
| Beam current | 150 µA | ✅ | arXiv:1802.04759 |
| Polarization | 85% longitudinal, ΔP/P ≤ 0.5% | ✅ | arXiv:1802.04759; [arXiv:2402.01027](https://arxiv.org/html/2402.01027v1) |
| Target | Liquid H₂, **600 mm long**, r = 25 mm; L ≈ 2.4×10³⁹ cm⁻²s⁻¹ | ✅ | arXiv:1802.04759 §5.1.2 |
| Solenoid | Superconducting, **B_z ≈ 0.6 T**; FOPI field map used in design studies; real magnet delivered 11/2024 | ✅ | arXiv:1802.04759 §5; [P2 news](https://www.blogs.uni-mainz.de/fb08p2/) |
| Forward acceptance | θ = 25°–45°, 2π azimuth, ⟨Q²⟩ ≈ 6×10⁻³ (GeV/c)² | ✅ | arXiv:1802.04759 |
| Forward detector | 82 wedged fused-silica bars, r = 450–900 mm, integrating readout | ✅ | arXiv:1802.04759 §5.2 |
| Forward tracker | HV-MAPS, 4 layers / 2 double planes | ✅ | arXiv:1802.04759 §5.5 |
| Run time | 10,000 h; production from ~2026 | ✅ | arXiv:1802.04759; [Baunack HP2030 talk](https://indico.ijclab.in2p3.fr/event/10641/contributions/35364/attachments/24117/35046/Baunack_Hadron_Physics_2030.pdf) |

## 2. BASKET — the backward-angle program

**BASKET = "Backward Scattering Electron Tracker"** (IRFU/CEA-Saclay,
M. Boonekamp et al.). ✅

- **Physics:** PV asymmetry at backward angles → strange magnetic FF **G_M^s**
  and effective axial FF **G_A^{p,Z}**; reduces the hadronic-structure input
  uncertainty of the forward sin²θW extraction by ~4×. Options in the 2018
  paper: parallel running at 200 MeV (Q² = 0.1 GeV², A ≈ 7.5 ppm) or dedicated
  2×1000 h H/D runs at 150 MeV (Q² = 0.06 GeV²). ✅ arXiv:1802.04759 §7.3
- **Acceptance: θ = 140°–150°, full azimuth.** ✅ (same)
- **Configuration:** "three detection layers, each housing **7,680 pads**
  (~23,000 channels total)", trigger-less streaming, **15–20 kHz/pad**,
  ~**100 MHz** total backward-electron rate, design capability 100 GHz. ✅
  [MPGD 2026 abstract, Kallitsopoulou et al.](https://indico.cern.ch/event/1473492/contributions/7105990/);
  [P2 site](https://www.blogs.uni-mainz.de/fb08p2/the-p2-experiment/)
  - ⚠️ 6 wedges × 1,280 pads = 7,680 pads/layer, ×3 = 23,040 — **the
    3-wheels-of-6-wedges picture is arithmetically confirmed**; wheel
    z-positions and the r = 95–650 mm span are not public (internal design).
- **Technology:** low-material-budget **metallic (non-resistive) bulk
  Micromegas**; background rejection by 3-layer hit coincidence; **VMM3a**
  front-end (62 fC–2 pC range, sub-2.5 ns timing, expected signal ~21 fC,
  pad capacitance 75–170 pF); size-1 prototype tested at CERN SPS H4; first
  beam test at MAMI 06/2023. ✅ (MPGD 2026; P2 site; Baunack talk)
  - ⚠️ **"sub-2.5 ns timing" is the *front-end* figure, not the detector
    time resolution.** Our prediction for a pad *hit* is **σ_t ≈ 10–13 ns
    (argon mixtures) / 15–25 ns (neon mixtures)**, set by primary-ionization
    statistics as σ_t ≈ (1.0–1.5)/(n_p·v_d); the electronics terms (ENC
    jitter, TAC quantization) add < 3 ns in quadrature. The two numbers are
    not in conflict, but they must not be quoted interchangeably — and any
    coincidence-window or rate argument has to use the detector number.
    See `research/TIME_RESOLUTION_NOTES.md`.
- **Requirements** (ANR project BASKET-P2, ANR-23-CE31-0025): ±1° scattering
  angle, **2% momentum resolution**, 100 MHz readout, ~10¹⁴ events in ~1000 h.
  ✅ [ANR page](https://anr.fr/Projet-ANR-23-CE31-0025)
- Thesis to watch: M. Gailliard (U. Paris-Saclay, 2024–2027, Boonekamp/
  Vandenbroucke) [theses.fr/s405038](https://theses.fr/s405038); first BASKET
  publication expected as MPGD 2026 proceedings. ❌ no NIM/JINST paper yet.

## 3. Kinematics at the wedges (cross-check of "~100 MeV at ~10°")

Elastic recoil E′ = E/(1 + (E/M)(1−cosθ)) (⚠️ computed, endpoints reproduce
the paper's Q² values exactly):

| Beam | Target | θ=140° | 145° | 150° |
|---|---|---|---|---|
| 155 MeV | p | 120.0 | 119.2 | 118.5 MeV |
| 200 MeV | p | 145.3 | 144.1 | 143.1 MeV |
| 155 MeV | ¹²C | 151.3 | 151.2 | 151.1 MeV |

→ **"~100 MeV electrons" confirmed (119–145 MeV).**

**Incidence angle** (⚠️ computed): in a uniform 0.6 T solenoid the angle to
the beam axis is *conserved* at 180°−θ = **30°–40°**, with Larmor radius
R_L ≈ 380 mm and radial excursions up to 2R_L ≈ 630–920 mm (consistent with
the wedge span). A **~10° incidence requires fringe-field collimation**:
sin α = sin α₀·√(B/B₀) (adiabatic invariant), so α = 10° ⇒ **B at the wheels
≈ 45–70 mT** (B/B₀ ≈ 0.07–0.12), i.e. the wheels sit outside/behind the
solenoid bore — consistent with 2021 meeting notes preferring mounting
outside the vacuum chamber
([indico.in2p3.fr/event/23632](https://indico.in2p3.fr/event/23632/?note=2669)).
Caveats: adiabaticity is marginal (R_L ~ gradient scale), the angle has a
large spread with θ and vertex-z over the 60 cm target — **use the real field
map; confirm wheel z and local B with the collaboration.** ❌ not public.

## 4. Backgrounds & rates (Geant4 campaign inputs)

- Dominant backward background: **photons, 50–150 keV** signals; mitigated by
  layer coincidence + a photon shield (2024 CEA internship studied ~10 keV
  X-rays + shield). ✅ MPGD 2026; [internship PDF](https://irfu.cea.fr/Phocea/stages/pdf.php?id_stage=405)
- Rate budget: ~100 MHz real electrons vs 100 GHz capability → designed for
  ~10³ photon/electron ratio; 15–20 kHz/pad ⇒ occupancy ~10⁻³ per 100 ns (⚠️).
- Forward-region benchmarks for normalizing photon generators
  (✅ arXiv:1802.04759 Table 12): 1.54×10¹² hits/s on the quartz ring
  (1.14×10¹² bg photons vs 7.1×10¹⁰ signal e⁻); first tracker plane sees up to
  10⁶× more bremsstrahlung photons than signal; quartz dose ~80 Mrad/10kh.
  Backward-specific dose: ❌ not published.
- **Møller electrons** are kinematically forward (θ_lab < 90°) and cannot reach
  140°–150° directly (⚠️); low-energy Møllers spiraling along field lines +
  target bremsstrahlung are what to generate for backgrounds.

## Bottom line for this simulation

> **Working assumption adopted 2026-08-05 (Dylan):** the magnetic field
> exists around the target only; the backward wheels are **outside** it.
> So no B-field in this simulation — electrons arrive as straight tracks,
> centrally at ~10° from the wedge normal, and the campaign scans the
> incidence angle 0°–40° (and energy) to bracket the uncertainty.
> See `SIM_CAMPAIGN_BRIEF.md` §2.

- Primary particle: **electrons, ~119 MeV** (155 MeV beam baseline; 144 MeV
  for the 200 MeV option).
- **Incidence angle is the key open kinematic parameter**: ~10° if the wheels
  sit in the fringe field (B ≈ 50 mT) as believed, but 30–40° inside the bore.
  Plan the campaign to scan incidence 0°–40°; get wheel z-positions + field
  map to pin it down.
- Rates: 15–20 kHz/pad signal, photon background (50–150 keV) at ~10³× —
  a photon-response simulation of the wedge is a natural campaign component.
- Not public (ask Saclay): wheel z, gas mixture, drift/amp details (we have
  those confirmed internally: 3 mm / 150 µm), photon-shield design.
