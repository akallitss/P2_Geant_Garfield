# Research notes — neon-based gas mixtures for the P2 Micromegas

**Compiled 2026-08-05** (web research pass for the campaign plan; summarized
into `../SIM_CAMPAIGN_PLAN.md` §7 — this file keeps the full detail and
citations). Unverified items flagged inline.

---

## 1. Neon premix availability

**Bottom line: pure Ne 4.5/5.0 is a stock catalog item everywhere; no major
supplier lists a neon detector premix (Ne/CO2, Ne/CH4, Ne/iC4H10) as an
off-the-shelf product — all handle them as certified custom mixtures
(2–4 week lead time, mixing surcharge).** Argon premixes (P10, Ar/CO2) are
stock by contrast.

- **Air Liquide**: pure Ne multiple grades
  ([product page](https://uk.airliquide.com/gases-and-products/neon)); custom
  via ALPHAGAZ/CRYSTAL + [mixture configurator](https://mixtureguide.airliquide.de/emixture/gas-mixture);
  Airgas (US) demonstrably blends Ne-balance mixtures.
- **Linde**: HiQ Neon 4.5/5.0
  ([datasheet](https://static.prd.echannel.linde.com/wcsstore/DK_REN_Industrial_Gas_Store/pdf/Linde_Neon_product_datasheet_EN_tcm575-579344.pdf));
  custom via [HiQ specialty mixtures](https://www.linde-gas.com/what-we-offer/gases/specialty-gases/specialty-gas-mixtures).
- **Messer**: Neon 5.0 ([page](https://specialtygases.messergroup.com/neon));
  [custom mixture configurator](https://www.messer.de/individuelle-spezialgas-gemische).
- **Westfalen**: custom mixing to 6.0 purity; neon stock unverified.
- **Nippon Gases**: explicitly "prepares neon mixtures"
  ([page](https://eu.nipponsanso.com/it/en/gases/neon/)). Air Products: pure
  Ne + excimer premixes; detector mixes custom.

**Neon mixtures used at scale:**

| Experiment | Mixture | Reference |
|---|---|---|
| ALICE TPC (~90 m³) | Ne/CO₂/N₂ 90-10-5 (Run 1, Run 3+) | [arXiv:2012.09518](https://arxiv.org/abs/2012.09518), [arXiv:1001.1950](https://arxiv.org/abs/1001.1950) |
| COMPASS Micromegas (2002-22) | Ne/C₂H₆/CF₄ 80/10/10 | [NIM A 469 (2001) 133](https://www.sciencedirect.com/science/article/abs/pii/S0168900201007690) |
| PICOSEC | Ne/C₂H₆/CF₄ 80/10/10 | [arXiv:2303.18141](https://arxiv.org/abs/2303.18141), [arXiv:1712.05256](https://arxiv.org/abs/1712.05256) |
| NA48/KABES | Ne/C₂H₆/CF₄ 79/11/10 @ 3×10⁵ Hz/mm² | Titov review [physics/0403055](https://arxiv.org/abs/physics/0403055) |
| NEWS-G | Ne + 0.7 % CH₄ @3.1 bar; Ne + 10 % CH₄ (S140) | [arXiv:1706.04934](https://arxiv.org/abs/1706.04934), [arXiv:2205.15433](https://arxiv.org/abs/2205.15433) |
| NA49 / CERES | Ne/CO₂ 90/10, 80/20 | [NIM A 430 (1999)](https://www.sciencedirect.com/science/article/abs/pii/S0168900298001107), [arXiv:0802.1443](https://arxiv.org/abs/0802.1443) |
| sPHENIX TPC | Ne/CF₄ 50/50 | [arXiv:2110.02082](https://arxiv.org/abs/2110.02082) |

Corrections: LHCb RICH uses C₄F₁₀/CF₄ (no Ne); RPC eco-mixes are He/CO₂/HFO;
T2K TPC is Ar-based. CERN EP-DT mixes all LHC gases on-site from pure
components with MFC modules ([gas systems](https://ep-dep-dt.web.cern.ch/gas-systems));
ALICE conserves Ne by closed-loop recirculation; "CERN neon recovery plant"
unconfirmed/likely false (recuperation targets CF₄/C₄F₁₀/SF₆/R134a,
[NIM A 2024](https://www.sciencedirect.com/science/article/pii/S0168900224007150)).

**Supply/cost**: 2022 shortage real — Ingas (Mariupol) + Cryoin (Odesa)
~45-54 % of semiconductor-grade neon halted 03/2022
([Reuters](https://www.business-standard.com/article/international/ukraine-halts-half-of-world-s-neon-output-for-chips-clouding-outlook-122031200039_1.html);
[CSIS](https://www.csis.org/blogs/perspectives-innovation/russias-invasion-ukraine-impacts-gas-markets-critical-chip-production));
China spot ~$63 → ~$252/m³; normalized by 2024. Bulk Ar ~$1/m³ vs pre-shortage
Ne ~$60/m³ → **Ne ~50-100× Ar per m³**; cylinder Ne 5.0 more. For a wedge
flushing ~1-5 l/h: real but not prohibitive; argues for modest flow or
recirculation. No reliable 2024-26 price series — get quotes.

## 2. Flammability (ISO 10156:2017)

LFL in air: iC₄H₁₀ 1.5 vol% (EU/EN 1839) / 1.8 % (US NFPA); C₂H₆ 2.4/3.0 %;
CH₄ 4.4/5.0 % (method difference, not conflict).

**Governing standard: ISO 10156:2017** ([ISO](https://www.iso.org/standard/66752.html);
coefficients in [EIGA Doc 169](https://www.eiga.eu/uploads/documents/DOC169.pdf)).
Inert → N₂-equivalent via Kk; non-flammable in air in any proportion iff
Σ(A′ᵢ/Tcᵢ) ≤ 1 (Tcᵢ = max fraction of flammable i in N₂ that stays
non-flammable). 2017 coefficients: **Kk: N₂ = 1, CO₂ = 1.5, Ar = 0.55,
Ne = 0.7, He = 0.9** (Ne slightly better inertant than Ar); **Tci:
iC₄H₁₀ = 3.4 %, C₂H₆ = 4.5 %, CH₄ = 8.7 %** (1996 edition was looser:
5.7/7.6/14.3 with noble Kk = 0.5).

Max quencher for non-flammability, **calculated** from the verified
coefficients (the standard prints only coefficients):

| Mixture | max quencher (ISO 2017) | (1996 ed.) |
|---|---|---|
| Ar/iC₄H₁₀ | **~1.9 %** | ~2.9 % |
| Ne/iC₄H₁₀ | **~2.4 %** | ~2.9 % |
| Ar/C₂H₆ | ~2.5 % | ~3.9 % |
| Ne/C₂H₆ | **~3.2 %** | ~3.9 % |
| iC₄H₁₀ in CO₂ | ~5.0 % | ~8.2 % |
| CH₄ in Ne | ~6.2 % | — |

Consequences:
- **Ar/iC₄H₁₀ 95/5 and all Ne/iC₄H₁₀ ≥5 % are unambiguously flammable**
  (operable — COMPASS ran flammable gas 20 years — but needs flammable-gas
  infrastructure: sniffers, zoning, interlocks).
- **ATLAS NSW Ar/CO₂/iC₄H₁₀ 93/5/2 is non-flammable — barely**: N₂-equiv
  inert = 0.55·93 + 1.5·5 = 58.65; A′ = 2/60.65 = 3.30 %; 3.30/3.4 = 0.97.
  The 5 % CO₂ does real work. ATLAS's "<3 % isobutane is non-flammable"
  statement ([arXiv:2411.17202](https://arxiv.org/abs/2411.17202);
  [ATL-MUON-PROC-2021-009](https://cds.cern.ch/record/2783915)) is uncited
  and matches the superseded 1996 edition.
- **Ne/CO₂/iC₄H₁₀ 95/3/2 passes with more margin**: A′ = 2/(2+0.7·95+1.5·3)
  ≈ 2.8 %, ratio 0.83.
- CERN rule: **Flammable Gas Safety Code G**
  ([index](https://safetyguide.web.cern.ch/SafetyGuide/Part1/03.6.html));
  internals SSO-blocked (unverified). ISO 10156 allows a spark-tube **test
  to override the calculation** for marginal mixtures.
- RPC precedent for ratio arguments: Daya Bay/BaBar ternary diagrams
  ([note](http://kirkmcd.princeton.edu/dayabay/rpc_gas_system2.pdf)).

## 3. Gain / operation of Ne mixtures in Micromegas

- **Iguaz et al., JINST 7 (2012) P04007, [arXiv:1201.3012](https://arxiv.org/abs/1201.3012)**
  (50/25 µm microbulk, Ar vs Ne): **max gain before sparking ~5× higher in
  Ne** (10⁵ vs 2×10⁴ at 50 µm); gain >5×10⁴ for iC₄H₁₀ >7 %; better energy
  resolution in Ne (10.5 vs 11.7 % FWHM @5.9 keV). **Nuance: Ne needed
  slightly *higher* amplification field for equal gain** (75 vs 72 kV/cm at
  gain 10⁴). → The "same gain at 40-80 V lower voltage" lore is **not
  verified** for 128-150 µm bulk; the solid Ne advantages are max-gain
  headroom and spark margin. No published Ne/iC₄H₁₀ or Ne/CO₂ gain curve for
  a 128-150 µm bulk exists — **a ⁵⁵Fe/⁹⁰Sr gain scan would be publishable**.
- **PICOSEC eco-gas campaign (Brunoldi, 6th DRD1 meeting, Oct 2025,
  [slides](https://events.camk.edu.pl/event/124/contributions/1187/attachments/813/2190/drd1_presentation.pdf))**:
  128 µm two-stage MM; tested **Ne/iC₄H₁₀ 95/5, 94/6, 90/10, 85/15, 75/25
  and Ne/iC₄H₁₀/CO₂ 90/5/5, 85/5/10, 90/3/7**: gains >10⁵ at 320-460 V,
  Ne/iso 95/5 ≈ 90/10 ≈ Ne/iso/CO₂ 90/5/5 gain curves nearly identical,
  CO₂-containing ones declared non-flammable. Closest existing dataset to
  the P2 candidate list. COMPASS/PICOSEC heritage: Ne/C₂H₆/CF₄ 80/10/10,
  gain 10⁵-10⁶ at 450-525 V ([arXiv:1712.05256](https://arxiv.org/abs/1712.05256)).
- **Discharges**: at fixed gain, discharge probability in Ne/CO₂ orders of
  magnitude below Ar/CO₂ (GEMs); critical charge 7.3×10⁶ e (Ne-CO₂ 90-10) vs
  4.7×10⁶ e (Ar-CO₂ 90-10) — Gasik et al.
  [arXiv:1704.01329](https://arxiv.org/abs/1704.01329); ALICE Ne-CO₂-N₂
  discharge prob ~10⁻¹⁰ at gain 2000
  ([arXiv:1807.02979](https://arxiv.org/abs/1807.02979)). Driver — lower
  charge density in Ne — applies to the x-ray-heavy P2 environment.
- **Penning**: Ne* 16.6-16.7 eV ionizes every candidate quencher (iC₄H₁₀
  10.6, C₂H₆ 11.5, CH₄ 12.6, CO₂ 13.8, N₂ 15.6 eV); Ar* (11.6 eV) cannot
  Penning-ionize CO₂/N₂. r_P(Ne/CO₂ 90/10) = 0.577 vs r_P(Ar/CO₂ 90/10) =
  0.456 @1 atm ([Garfield++ Penning tutorial](https://garfieldpp.docs.cern.ch/tutorials/penning/);
  Şahin/Veenhof JINST 5 (2010) P05002, JINST 16 (2021) P03026). **No
  published r_P for Ne/iC₄H₁₀ or Ne/CH₄** — simulations need assumed/fitted
  values.

## 4. Primary ionization numbers (PDG 2024, Table 35.5)

Source: [PDG 2024 §35.6 (Sauli & Titov)](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-particle-detectors-accel.pdf);
caption warns values vary between compilations.

| Gas | W_I [eV] | dE/dx min [keV/cm] | n_p [/cm] | n_t [/cm] |
|---|---|---|---|---|
| Ne | 37 | 1.45 | 13 | 40 |
| Ar | 26 | 2.53 | 25 | 97 |
| CH₄ | 30 | 1.61 | 28 | 54 |
| C₂H₆ | 26 | 2.91 | 48 | 112 |
| iC₄H₁₀ | 26 | 5.67 | 90 | 220 |
| CO₂ | 34 | 3.35 | 35 | 100 |
| CF₄ | 35-52 | 6.38 | 52-63 | 120 |

Per mm volume-weighted n_p: Ar 95/5 ≈ 2.8; Ne 95/5 ≈ 1.7; Ne 90/10 ≈ 2.1;
Ne 80/20 ≈ 2.8 → **Ne/iso 80/20 fully recovers Ar 95/5 primary statistics;
90/10 recovers most**. PICOSEC cross-check: ~21 clusters/cm for
Ne/C₂H₆/CF₄ 80/10/10. The "~30 % fewer primaries" penalty is
quencher-recoverable — the trade is flammability (§2).

## 5. CF₄

- **Aging/etching**: F radicals/HF (catalyzed by H₂O) etch glass, FR4,
  Al/Cu; keep H₂O < ~1000 ppm, no glass/Si in gas path (Titov
  [physics/0403055](https://arxiv.org/abs/physics/0403055); DESY 2001 aging
  workshop [CDS 546220](https://cds.cern.ch/record/546220)). LHCb saw FR4
  etching in Ar/CF₄/CO₂ — **a bulk MM is an FR4 board**. LHCb GEMs survived
  ~200 mC/mm² in Ar/CO₂/CF₄ ([Alfonsi et al.](https://www.sciencedirect.com/science/article/abs/pii/S0168900204003407)).
  Micromegas-specific CF₄ corrosion under-characterized (Titov).
  Counterpoint: CMS CSC sees carbon deposits *appear* when CF₄ < 2 % —
  anti-polymerization role real ([arXiv:2402.04181](https://arxiv.org/abs/2402.04181)).
- **Regulatory**: GWP 7,380, EU F-gas Reg. (EU) 2024/573
  ([text](https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:32024R0573));
  containment/reporting, no HFC-style quota; lab pathway exists (article
  unverified). CERN: −28 % Scope 1 target, CF₄ recuperation (CMS ~70 %)
  ([EP News](https://ep-news.web.cern.ch/content/reducing-greenhouse-gas-emissions-particle-detection-systems-cern-status-and-perspectives)).
- **When justified**: fast/saturated low-field drift (sPHENIX Ne/CF₄ 50/50:
  ~8 cm/µs @400 V/cm), ps timing, anti-aging at high charge. For P2: only
  the streaming drift-speed argument; a small T2K-style ~3 % CF₄ fraction
  buys most of the velocity with less baggage.

## 6. Candidate ranking additions

1. **Ne/CO₂/iC₄H₁₀ 95/3/2** — NSW-analog; passes ISO with margin (0.83);
   Penning-active on both quenchers; 3-channel mixing.
2. **Ne/iC₄H₁₀/CO₂ 90/5/5 or 90/3/7** — PICOSEC eco-mixtures; only
   published-adjacent non-flammable Ne mixes with measured MM gain curves;
   arguably strongest new candidate.
3. **Ne/CO₂/N₂ 90/10/5** — ALICE; boring/reliable/non-flammable, 25 y
   heritage; N₂ quenches Ne ~17 eV states, immunizes against air-leak N₂
   ([NIM A 622 (2010)](https://www.sciencedirect.com/science/article/pii/S0168900210008910)).
   Drawback: slow unsaturated drift (2.83 cm/µs @400 V/cm,
   [NIM A 548 (2005)](http://www-alice.gsi.de/tpc/papers/ref/tpc_nim548.pdf))
   → ~110 ns over 3 mm, v_d sensitive to E/T/composition. Binary Ne/CO₂
   90/10 shares this + marginal quenching alone.
4. **Ne/CH₄ ~93/7** — proportional-counter classic, NEWS-G heritage, cheap,
   fast low-field drift; 93/7 passes ISO calculation cleanly (10 % is near
   the ~6.2 % limit → would need test-route or supplier cert — flagged).
5. **Ne/C₂H₆ 90/10** — COMPASS chemistry minus CF₄; on PICOSEC's announced
   list, no published binary results yet; flammable above ~3.2 %; ethane
   premixes fine at 200 bar (no condensation).
6. **Ne/CF₄ 90/10** — non-flammable, fast; §5 baggage; keep only if drift
   speed proves decisive.

Context surveys: eco-MPGD [arXiv:1512.08542](https://arxiv.org/abs/1512.08542);
LHC mixture transport compendium [arXiv:1110.6761](https://arxiv.org/abs/1110.6761).

## 7. Gas mixers & premix constraints

- 2-3 channel MFC racks standard (CERN GDD/RD51 lab
  [arXiv:1806.09955](https://arxiv.org/abs/1806.09955); EP-DT mixes all LHC
  gases on-site). Bronkhorst EL-FLOW ±(0.5 % Rd + 0.1 % FS)
  ([datasheet](https://products.bronkhorst.com/media/cp3lfvo1/el-flow-select_en.pdf));
  Vögtlin red-y ±0.3 % FS; MKS/Brooks ~1 % class. Sized near full scale →
  quencher stable to ±0.05-0.1 points at 10 % (few-% gain effect, like P/T
  drifts); calibrate as ratio. ~2-4 k€/channel.
- **Isobutane vapor pressure 3.1 bar abs @21 °C
  ([Air Liquide encyclopedia](https://encyclopedia.airliquide.com/isobutane)):
  in a 200 bar cylinder iC₄H₁₀ condenses above ~1.5 % → no full-pressure
  premix for ≥5 % iso mixtures** (reduced-pressure fills only). CO₂, CH₄,
  C₂H₆, CF₄ premixes fine at full pressure.
- Cylinder stratification is a myth (LLNL-TR-812248,
  [OSTI](https://www.osti.gov/servlets/purl/1637588)); real caveats: initial
  homogenization, condensation. Certified mixes per ISO 6142-1; "±2 %
  relative" tolerance is contractual, not standard (flagged).
- **Decision**: 3-channel mixer ≈ cost of 2-channel + one MFC and unlocks
  the whole ternary space + in-situ quencher scans (wanted anyway given no
  published 150 µm bulk Ne gain curves). Premix only for a frozen
  non-condensable binary (Ne/CO₂, Ne/CH₄ ≤ few %).

## Key flags

1. "Ne = lower operating voltage": unverified; microbulk data show opposite
   in field terms. Verified: ~5× max-gain headroom, much lower discharge
   probability.
2. No published 128-150 µm bulk gain curves for Ne/iC₄H₁₀ or Ne/CO₂.
3. No published Penning r_P for Ne/iC₄H₁₀, Ne/CH₄; no citable low-field v_d
   tables for Ne/iC₄H₁₀, Ne/C₂H₆ → Magboltz.
4. CERN Code G internals SSO-blocked; §2 percentages are calculations from
   verified coefficients; ATLAS "<3 %" uncited/1996-edition-like.
5. Unverified: Westfalen Ne stock, EU F-gas lab-exemption article, EU CF₄
   pricing, Ne-balance CH₄ CGA classification, ±2 % premix tolerance.
