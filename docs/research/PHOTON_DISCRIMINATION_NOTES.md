# Energy-deposition discrimination against the 50–100 keV photon background

**Written:** 2026-08-05 (Dylan, after discussion with colleagues: the
expected background band is **50–100 keV with a higher-energy tail**, which
is far above the Ar K-edge, so the "neon rejects photons" premise needs
re-examining; the proposed replacement handle is **per-event energy-deposition
discrimination** — electrons are MIPs with a Landau spectrum, converted
photons should deposit more).

**Status:** analytic pass done **and confirmed in Geant4** (2026-08-05,
4×10⁵ photons at 60 keV in each of Ar/Iso and Ne/Iso on lxplus — see §4.3a).
The arithmetic here was written as the target the campaign had to reproduce
or contradict; it reproduced. Numbers in §4.3a were re-measured after the
P0.1 gas-composition fix (`GAS_FIX_NOTES.md`); the pre-fix pass is
superseded. Reproduce the analytic side with
`python3 scripts/photon_budget.py` (NIST attenuation data transcribed in
the script header, Klein-Nishina integrated numerically, Katz-Penfold
electron ranges).

Companion docs: `../SIM_CAMPAIGN_PLAN.md` (§0, §1, §4.2, §5 amended from
this note), `TIMING_PSD_NOTES.md` (the time-structure handle, which has the
same blind spot identified in §4 below), `GAS_NOTES.md` (gas physics and
procurement).

---

## 1. Summary — the three findings

1. **The premise is directionally right but the crossover is at ~80 keV.**
   Argon is *not* "the same as neon" at 50 keV — it still has 5.4× the
   interaction probability. But the ratio falls fast: 4.3× at 60 keV, 3.1×
   at 80 keV, 2.5× at 100 keV, 2.1× at 150 keV. Above ~100 keV the two
   gases differ only by their electron densities (Compton), which is a
   factor 1.8 and cannot be improved on.

2. **Energy-deposition discrimination is real, but weaker than the plan
   currently assumes, and its target population is the *minority* one.**
   A gas conversion in the band deposits ~3–30× the MIP charge — good
   separation, but not the "10–100×" the plan inherited from ⁵⁵Fe
   intuition, because at 50–100 keV the photoelectron's range is
   **1–8 cm of gas**, so it crosses the 3 mm gap instead of stopping in it.

3. **The finding that reframes the campaign: conversions in the *walls*
   outnumber conversions in the gas by ~15× (Ar) to ~65× (Ne),** and those
   are exactly the events that look most MIP-like. The copper pad plane and
   the steel mesh each convert ~1–3 % of incident 60 keV photons, and a few
   percent of those photoelectrons escape into the gas. Since that floor is
   nearly gas-independent while the *signal* is ~40 % weaker in neon, the
   honest conclusion is:

   > At 50–100 keV, switching argon → neon buys almost no photon rejection
   > (**measured: (12 ± 3) % fewer photon-induced hits**, against 77 % if
   > only the gas mattered) while costing real electron efficiency. **Argon
   > should be treated as the baseline, not the fallback.** The rejection
   > has to come from the 3-layer coincidence, with the charge window as a
   > secondary polish.

   **This was written as an analytic conclusion resting on a crude escape
   model (§4.2); Geant4 confirmed it on 2026-08-05** — walls are 90 % of
   fakes in argon and 97 % in neon, the mesh conversion rate agrees between
   the two gases to 1.4 % (the internal check that the floor really is
   gas-independent), and the measured conversion budget matches the
   analytic estimate within a factor ~2 (§4.3a). §6 lists the recording
   changes that made the measurement possible; they are all implemented.

---

## 2. Interaction probability — argon vs neon in the band

Per incident photon, 3 mm gap, pure gas at 1 atm / 20 °C. μ/ρ from NIST
(Hubbell & Seltzer); the photoelectric fraction is
`(μ/ρ − μ_incoh/ρ)/(μ/ρ)` with the incoherent part from Klein-Nishina, so
it lumps coherent (Rayleigh) scattering in with "photoelectric" and is an
**upper bound** — Rayleigh deposits nothing at all.

| E [keV] | P_int (Ar) | P_int (Ne) | ratio | f_PE (Ar) | f_PE (Ne) |
|---|---|---|---|---|---|
| 50 | 3.49×10⁻⁴ | 6.49×10⁻⁵ | 5.38 | ≤78 % | ≤35 % |
| 60 | 2.32×10⁻⁴ | 5.44×10⁻⁵ | 4.27 | ≤68 % | ≤25 % |
| 80 | 1.38×10⁻⁴ | 4.48×10⁻⁵ | 3.07 | ≤49 % | ≤13 % |
| 100 | 1.02×10⁻⁴ | 4.03×10⁻⁵ | 2.53 | ≤35 % | ≤8 % |
| 150 | 7.11×10⁻⁵ | 3.45×10⁻⁵ | 2.06 | ≤16 % | ≤3 % |

Notes:

- The Ar photoelectric coefficient scales as E⁻³ across the whole band
  (checked: 0.30 cm²/g at 60 keV → 2.53 predicted vs 2.53 measured at
  30 keV, 0.070 predicted vs 0.070 at 100 keV). So there is no structure
  to exploit — no edge, no resonance, just a steep monotone fall.
- **Argon is still majority-photoelectric at 50–60 keV.** The colleagues'
  "past the K-edge so they're the same" is an over-correction at the bottom
  of the band; it becomes accurate around 100 keV.
- Ne's advantage is capped at 1.83× (the ratio of electron densities at
  equal pressure) once Compton dominates. Nothing in the gas menu changes
  that — it is set by Z/A and density.
- **Where the tail lands matters a lot.** If the spectrum is soft (peaked
  near 50 keV) the Ar/Ne difference is a real factor ~5; if it is hard
  (100 keV+) it is a factor ~2. Getting the actual background spectrum from
  the BASKET/MESA side is now a higher-value input than another gas.
  ⬜ **NEEDED: the incident photon spectrum, not just the 50–150 keV band
  label.**

## 3. What a conversion actually deposits in 3 mm

MIP reference, **measured** (16 000 e⁻ at 119 MeV through 3 mm, §4.3a):
**Ar/Iso 1067 eV mean / 738 eV median; Ne/Iso 603 eV mean / 354 eV median.**
Above the PDG dE/dx|min figures (774 / 475 eV at these densities) because
119 MeV electrons sit on the relativistic rise and the Landau tail pulls the
mean up; the median is the better MPV proxy, and the *spread* is what the
ROC actually integrates over.

Compton electron energies at these photon energies are small: mean transfer
is 9.4 % of the photon energy at 60 keV, 13.8 % at 100 keV (Klein-Nishina,
verified against σ = 0.2865 b/e⁻ at 511 keV and 0.2112 b at 1 MeV). So:

| channel | at 60 keV | at 100 keV | behaviour in a 3 mm gap |
|---|---|---|---|
| Compton | ⟨T⟩ = 5.6 keV, T_max = 11.4 keV | ⟨T⟩ = 13.8 keV, T_max = 28.1 keV | **contained** below ~80 keV (practical range 0.14 mm in Ar at 60 keV) → deposits the full electron energy, ~7× MIP at 60 keV, ~14–18× at 100 keV |
| photoelectric | 60 keV e⁻, practical range **16.6 mm** | 100 keV e⁻, range **40 mm** | **crosses the gap** → deposits only 3–11 keV (Ar), i.e. 4–14× MIP, not 80× |
| Rayleigh | — | — | deposits **nothing**; must not be counted as an interaction |

Two consequences the current plan gets wrong:

- The plan's §0/§5 language — "a contained photoabsorption is 10–100× the
  MIP charge, so a lower threshold alone cannot reject it" — is right for
  ⁵⁵Fe (5.9 keV, genuinely contained, ~230 e⁻ point blob) and **wrong for
  the actual background band**, where the photoelectron leaves the gap.
  The real separation is more like **3–30×**, with the low end overlapping
  the MIP Landau tail. The upper charge cut still works; it just has less
  headroom than assumed, and ADC saturation is much less of a concern than
  the plan implies.
- **The Compton spectrum reaches zero.** dσ/dT is roughly flat from 0 to
  T_max, so at 60 keV about 1/11 of Compton events transfer under 1 keV and
  are indistinguishable from (or dimmer than) a MIP. These are irreducible
  by any charge cut. They are, however, mostly *below threshold* as well,
  so they cost acceptance rather than creating fakes — the sim has to
  resolve which.

## 4. The wall floor — why the gas may not matter

### 4.1 The numbers

Per incident 60 keV photon, using the stack as built in
`DetectorConstruction::ConstructP2()`:

| layer | P(interact) | f_PE | practical e⁻ range | P(e⁻ into gas) |
|---|---|---|---|---|
| PCB_Cu_F (18 µm Cu, 98.3 % coverage) | 2.5×10⁻² | 91 % | 3.1 µm | **2.0×10⁻³** |
| Micromesh (19 µm SS wire) | 7.8×10⁻³ | 87 % | 3.4 µm | **2.0×10⁻³** |
| Drift cathode mylar (2 × 12 µm) | 6.2×10⁻⁴ | 9 % | 20 µm | 2.4×10⁻⁵ |
| Front window mylar (40 µm) | 1.0×10⁻³ | 9 % | 20 µm | 2.4×10⁻⁵ |
| — drift gas, **Ar**, 3 mm | 2.3×10⁻⁴ | | | 2.3×10⁻⁴ |
| — drift gas, **Ne**, 3 mm | 5.4×10⁻⁵ | | | 5.4×10⁻⁵ |

**Wall total ≈ 4.0×10⁻³ per photon — 17× the argon gas rate, 74× neon's.**
At 100 keV the wall total falls to 1.9×10⁻³ and the ratios are 19× (Ar) and
47× (Ne). So:

| | total fake source per photon @60 keV |
|---|---|
| Ar/Iso | 4.0×10⁻³ (wall) + 2.3×10⁻⁴ (gas) = **4.3×10⁻³** |
| Ne/Iso | 4.0×10⁻³ (wall) + 5.4×10⁻⁵ (gas) = **4.1×10⁻³** |

A **4 % difference**, against a ~30 % loss in primary ionization and a
30× worse zero-cluster floor (plan §1: P₀ at 3 mm is 0.021 % for Ar/Iso vs
0.64 % for Ne/Iso 95/5). That is the whole argument in one table.

> **Correction, 2026-08-05 (later).** An earlier version of this table put
> the mesh at 3.5×10⁻³ and the wall total at 5.5×10⁻³, from a
> projected-wire-coverage estimate of the mesh areal density that was 1.7×
> too high. The numbers above use the correct plain-weave value
> (t_eff = 2·(πd²/4)/pitch = 8.46 µm solid-equivalent). The conclusion is
> unchanged — the wall/gas ratio moves from 24× to 17× in argon — but see
> §4.2, which is also corrected: **the Geant4 mesh model was not the thing
> that was wrong.**

The two big contributors are both *unavoidable by gas choice* and both
made of the highest-Z material in the detector. Note the asymmetry in where
their electrons land:

- **Mesh** photoelectrons go into *both* gaps. The ones emitted upward
  cross the full drift gap at full gain — the worst case: track-like
  topology, track-like timing, full gain.
- **Pad-plane Cu** photoelectrons are born in the 150 µm amplification gap
  and get **z-dependent partial gain** (birth at the pad surface ⇒ ~no
  gain, at the mesh ⇒ full gain). Their charge spectrum is therefore
  smeared over four decades and mostly suppressed — plan §2 already flags
  this and it is the one thing working in our favour here.

### 4.2 How much to trust this

The escape model is deliberately crude: uniform interaction depth,
straight-line practical range = 0.5 × CSDA, half the electrons emitted
toward the gas. Realistic uncertainty on the wall numbers is a factor 2–3.
They cannot plausibly be a factor 20 lower, which is what it would take to
restore the gas argument.

**On the Geant4 mesh model — checked 2026-08-05, and it is fine.** The code
builds the mesh as a 2·d-thick slab of steel at fill fraction π·d/(4·pitch)
(`DetectorConstruction::ConstructP2`). That is not an ad-hoc dilution: its
areal density is π·d²/(2·pitch) = 8.46 µm of solid steel, which is
*exactly* the plain-weave value for one x- and one y-wire per pitch² cell.
So the **conversion rate in the mesh is right**, and the earlier claim in
this note that the slab "over-counts interactions" was wrong — the error
was in this note's own arithmetic, not in the code.

What the slab does approximate is the *escape geometry*. An electron born
in a 22 %-density slab has the same range in g/cm² but ~4.5× the geometric
range in µm, and the slab is 2× thicker than a wire is wide, so the two
effects largely cancel: escape fraction ≈ 0.5·(14 µm/38 µm) = 0.18 for the
slab versus 0.25 for a solid 19 µm cylinder. **A ~1.4× underestimate, not
3×.** Two things the slab genuinely cannot represent: the 51 % optical
transparency (irrelevant for a thin absorber where P ≪ 1, but not for
field/transparency modelling) and the correlation between where a photon
converts and where the electron can get out along a wire.

Net: a cylindrical-wire cross-check is **worth doing and no longer
blocking** (plan P0.13, downgraded). Existing plan risk #3 stays a genuine
caveat rather than a decision-stopper.

### 4.3a Measured — Geant4, corrected gas *(2026-08-05, night)*

**These supersede an earlier pass taken before the P0.1 gas fix**, which had
argon 7 % too dense, neon 10 % too dense, and neon carrying 2.6× the
intended isobutane by mass (`GAS_FIX_NOTES.md`). Conclusions are unchanged;
the Ar-vs-Ne gap narrowed slightly.

4×10⁵ γ at 60 keV, normal incidence, 3 mm gap, per incident photon:

| class | Ar/Iso | Ne/Iso | median edep (Ar) | median edep (Ne) |
|---|---|---|---|---|
| **wall photoelectric** | **3.16×10⁻³** | **3.01×10⁻³** | **567 eV** | **370 eV** |
| gas photoelectric | 1.45×10⁻⁴ | 1.25×10⁻⁵ | 6125 eV | 13 102 eV |
| fluorescence re-absorption | 1.15×10⁻⁴ | 1.50×10⁻⁵ | 6350 eV | 7418 eV |
| gas Compton | 7.50×10⁻⁵ | 5.75×10⁻⁵ | 4530 eV | 3142 eV |
| wall Compton | 3.0×10⁻⁵ | 1.25×10⁻⁵ | 3150 eV | 3123 eV |
| **total** | **3.53×10⁻³** | **3.11×10⁻³** | | |

Conversion-layer budget (where the ionizing chain started), per photon:

| layer | Ar/Iso | Ne/Iso | analytic |
|---|---|---|---|
| Micromesh | 2.10×10⁻³ | 2.07×10⁻³ | 2.0×10⁻³ |
| PCB_Cu_F | 1.03×10⁻³ | 9.25×10⁻⁴ | 2.0×10⁻³ |
| DriftGas | 3.10×10⁻⁴ | 8.0×10⁻⁵ | 2.3 / 0.54×10⁻⁴ |
| everything else | ≤2.5×10⁻⁵ | ≤1.5×10⁻⁵ | ~5×10⁻⁵ |

The mesh number is **identical in the two gases to 1.4 %** (inside the ~3.5 %
counting error), which is the internal check that the wall floor really is
gas-independent. The analytic mesh estimate is dead on; the pad-copper
estimate was 2× high, the only place the escape model materially missed.
**Walls are 90 % of all fakes in argon and 97 % in neon.**

**The Ar-vs-Ne comparison, measured.** Total photon-induced hit rate:

| gas | per 4×10⁵ γ | per photon |
|---|---|---|
| Ar/Iso | 1410 | 3.53×10⁻³ |
| Ne/Iso | 1244 | 3.11×10⁻³ |

**A (12 ± 3) % reduction** — against a naive gas-only expectation of **77 %**.
Direct gas conversions do fall by the predicted factor (Ar 2.20×10⁻⁴ vs Ne
7.0×10⁻⁵, ratio 3.1 against 4.3 predicted), and the fluorescence chain by
7.7× as expected — but both are minorities riding on a gas-independent wall
floor. That is the whole argument, measured on the as-built geometry.

**The deposit spectra, and why argon separates slightly *better*.** MIP
reference from 16 000 e⁻ at 119 MeV in the same geometry:

| gas | MIP median | MIP mean | wall-PE median | wall / MIP |
|---|---|---|---|---|
| Ar/Iso | 738 eV | 1067 eV | 567 eV | **0.77** |
| Ne/Iso | 354 eV | 603 eV | 370 eV | **1.05** |

The analytic pass predicted wall photoelectrons would land at 3–5× the MIP
deposit. They do not — **they are MIP-sized or smaller**, because a
photoelectron escaping a metal surface has already lost most of its energy
there and then only clips the gap. Consequences:

- an **upper** charge cut does almost nothing against the dominant fake
  class; the **lower** threshold is the cut that bites it, and it trades
  directly against electron efficiency through the Landau tail;
- **argon is marginally better placed than neon here too**: its wall fakes
  sit at 0.77× the MIP median, neon's at 1.05× — i.e. in neon the dominant
  background is indistinguishable from signal by charge, while in argon it
  is at least slightly below. The wall deposit is set by the gas density
  just as the MIP is, but wall electrons are soft enough that their dE/dx
  does not scale the same way.

The cleanly separable classes (gas photoelectric ~8× MIP, fluorescence
~9× MIP, gas Compton ~6× MIP, all point-like) are only ~10 % of argon fakes
and ~3 % of neon's.

### 4.3 The Cu-K fluorescence channel — the one gas-sensitive wall channel

A photoabsorption in the copper pad plane leaves a K vacancy, which
fluoresces at 8.05 keV with yield 0.44. That photon is *soft*, so its
re-absorption is strongly Z-dependent: in 3 mm it is absorbed with
probability 5.9 % in Ar but only 0.6 % in Ne. Chain probability per
incident 60 keV photon: **1.8×10⁻⁴ in Ar, 1.7×10⁻⁵ in Ne** — comparable to
direct gas conversion in argon, and a genuine 10× argon penalty.

These produce the *classic* discriminable signature: a point-like ~8 keV
blob, ~10× MIP charge, δ-like in arrival time. So they are the population
that both the charge window and the timing handle kill most efficiently.
Their main cost is that they are ~5 % of all argon-side fakes before cuts.

**This channel is currently suppressed in the code by a bug** (§6.1) and
must be switched on before any of this is measured.

## 5. From energy deposited to ADC — what breaks the proportionality

Answering Dylan's second question directly. Ordered by how much damage each
does to "ADC ∝ energy deposited in the gas". The first group are the ones
that can move a conclusion; the last group are percent-level.

### 5.1 Order-unity or worse

1. **Where the deposit happened (binary).** Only `DriftGas` and `AmpGas`
   produce signal at all. Deposits in `FrontGas`, `DriftCathode_Gas` and
   `BackGas` are in field-free or reversed-field regions and yield exactly
   zero — the simulation currently doesn't even record them, so we cannot
   presently tell "no conversion" from "conversion somewhere useless".
2. **Partial gain in the amplification gap.** An electron born at depth z
   in the 150 µm gap is amplified only over the remaining distance:
   G(z) ≈ G_full^((z_mesh − z)/d_amp). Pad-surface birth ⇒ gain ~1;
   mesh-side birth ⇒ full gain. **Four decades of dynamic range within one
   volume**, and it is precisely the population (Cu photoelectrons) whose
   rate matters most. A flat full-gain treatment is not merely imprecise,
   it is qualitatively wrong.
3. **Space-charge gain suppression on dense deposits.** A contained 8 keV
   blob is ~300 primaries inside ~100 µm; at gain 10⁴ that is ~3×10⁶
   electrons in one avalanche region, approaching the Raether limit
   (~10⁷–10⁸) where the local field is screened and gain drops. This
   **compresses exactly the signal the discrimination relies on** — big
   deposits read back smaller than proportional — and at the top end it is
   also the spark mechanism. Micromegas literature documents gain vs
   primary-density curves; we have none for our geometry.
4. **Ballistic deficit / arrival-time spread.** A track's charge arrives
   spread over the full drift time (60–110 ns depending on gas); a point
   deposit arrives all at once. With peaking time t_p ≳ 200 ns the shaper
   integrates both fully and PDO ∝ charge; at t_p = 25 ns the track is
   suppressed ~3–4× and the point deposit is not. This breaks
   proportionality *in our favour* and is the subject of
   `TIMING_PSD_NOTES.md` — but note it is defeated by wall photoelectrons,
   which are track-like in time as well as in charge.
5. **Charge sharing across pads with a per-pad threshold.** A MIP spread
   over two pads can fall below threshold on both while a point deposit on
   one pad fires. Helps a max-pad discriminant, hurts a summed-charge one;
   neighbor logic partially recovers it. At the nominal 10° incidence a
   track averages only ~1.1 pads, so this is smaller than it sounds.
6. **ADC saturation** (10-bit PDO, ~8-bit effective). Rails on the largest
   deposits. Harmless for a *window* cut (saturated ⇒ rejected) but it
   destroys any shape information above the rail, and it interacts with
   gain choice. Less critical than the plan assumed, since §3 shows the
   band deposits 3–30× MIP rather than 100×.

### 5.2 Percent-level, but they bias Ar-vs-Ne comparisons specifically

7. **Penning transfer.** In neon mixtures Ne* (16.6 eV) ionizes every
   candidate quencher, so the *effective* W is lower than the tabulated
   ionization W and the same deposited energy yields more electrons. Any
   edep→primaries conversion that uses a single W per gas silently mis-scales
   Ne relative to Ar. Published r exists only for Ne/CO₂ 90:10 (0.577) —
   this is why N10 was added to the gas list.
8. **W-value treatment in the code.** `SteppingAction` uses
   `n = floor(edep/W)` plus a uniform-random remainder — that is a
   reasonable mean but gives Bernoulli, not Fano, fluctuations, and W is
   not constant below ~1 keV. For a Ne MIP with ~20 primaries in 3 mm this
   distorts the *width* of the signal distribution, which is exactly what
   the ROC integrates over.
9. **Pillars.** 4.8 % of the amplification-gap area is dielectric: no gain
   there, and drifting electrons landing on a pillar are lost. Roughly a
   uniform 5 % loss, plus rate-dependent charging-up.
10. **Mesh transparency** ε_mesh(E_drift/E_amp) — flat near 100 % at the
    field ratios we will run, so mostly a constant, but it is field-point
    dependent and therefore couples to the drift-field choice per gas.
11. **Attachment on O₂/H₂O impurities**: exp(−z/λ) with λ ≫ 3 mm at
    reasonable purity — sub-percent, but z-dependent, so it very slightly
    tilts the drift-gap response and mimics a topology effect.
12. **Ar K-fluorescence escape**: 12 % of argon photoabsorptions lose
    2.96 keV, which then escapes the gap (attenuation length ~40 cm below
    the K-edge). Produces a visible escape structure in argon and none in
    neon. Small in the 50–100 keV band, prominent for the ⁵⁵Fe bench runs.
13. **Gain uniformity across the wedge** — amp-gap variation from the bulk
    process and mesh sag under the overpressure bulge. Unknown; a
    prototype-measurement question, not a simulation one.

### 5.3 The answer to "do these deposits actually reach the pads"

For everything born in the **drift gap**: yes, essentially all of it —
items 9–11 cost ~5–10 % in total and are close to deposit-independent, so
proportionality holds well there. For everything born in the
**amplification gap**: no, and unpredictably so (item 2). For everything
born **outside both** (front gas, cathode gas, back gas): no, ever.

So the sim's job is less "does charge reach the mesh" and more **"classify
which of those three regions each deposit belongs to, and for the amp gap,
record z"**. That is a bookkeeping requirement, and the current output
cannot meet it.

## 6. What the simulation must record — gap analysis

### 6.1 A bug that has to be fixed first

`PhysicsList::SetCuts()` sets a **0.1 mm gamma production cut**. In copper
that corresponds to roughly 5–10 keV, i.e. **Cu-K fluorescence at 8.05 keV
sits right at the production threshold and is likely being suppressed** —
the exact channel §4.3 says is the dominant gas-sensitive wall process, and
the one whose presence Ar-vs-Ne comparisons hinge on. Fix:
`G4EmParameters::Instance()->SetDeexcitationIgnoreCut(true)` (or
`/process/em/deexcitationIgnoreCut true`), with fluorescence, Auger and
PIXE explicitly enabled. Plan item P0.6 already asks for deexcitation to be
on; it does **not** currently say to defeat the cut, which is the part that
actually matters.

Related: the 10 µm e⁻ range cut is ~1 keV in gas (fine) but ~80 keV in
copper (coarse). Photoelectric and Compton secondaries are produced
regardless of cuts, so escape physics survives, but the energy sharing
inside the metal is crude. A `G4Region` over the mesh + pad Cu + gas with a
~1 µm cut is the clean fix.

### 6.2 What is missing from the output

The current `IonizationCluster` (`include/EventData.hh`) carries
x, y, z, edep, nPrimary, trackID, parentID, volumeName, particleName, ke.
That is not enough to separate interaction types. Missing, in priority
order:

1. **Creator process of the depositing track** (`phot`, `compt`, `Rayl`,
   `eIoni`, `eBrem`, `conv`…). Without this, "was it Compton or
   photoelectric" is unanswerable — this is the headline ask.
2. **Inherited origin tag.** A δ-ray from a photoelectron has
   `creatorProcess == eIoni` and `parentID` pointing at the photoelectron,
   not at the conversion. Ancestry in the current output is one level deep,
   so the chain is broken. Needs a `G4VUserTrackInformation` carrying
   {origin process, origin volume, origin energy, ancestor trackID}, copied
   to every secondary at creation — which requires a
   **`TrackingAction` class that does not currently exist in this repo**
   (`ActionInitialization.cc` registers only primary/run/event/stepping).
3. **Per-event first-interaction record for the primary photon**: process,
   volume, position, energy before/after. This is what turns a photon run
   into the per-layer conversion budget of §4, and it is one small
   `G4UserTrackingAction`/stepping hook.
4. **Per-event region classification and the region-resolved edep sums** —
   `edepDrift` and `edepAmp` exist; add the currently-unscored gas volumes
   (`FrontGas`, `DriftCathode_Gas`, `BackGas`) and the solid layers so
   "converted somewhere useless" is distinguishable from "did not convert".
5. **z relative to the mesh for amp-gap clusters** — derivable in analysis
   from absolute z and the geometry, but it should be explicit so Stage B
   cannot get it wrong.
6. **A per-event summary tuple that makes the money plot a one-liner**:
   `(edepDrift, edepAmp, nPrimDrift, interactionClass, conversionVolume,
   nClusters, zSpread)` where `interactionClass` ∈ {none, gas-photoelectric,
   gas-Compton, gas-Rayleigh, wall-photoelectric, wall-Compton,
   fluorescence-reabsorption, MIP-crossing, other}. **This is the branch the
   requested histogram is made from** — "event-by-event energy deposition,
   separated by interaction type" is then a single `TTree::Draw` with a
   class cut, per gas.

### 6.3 Cheap validation that does not need any of the above

`scripts/photon_budget.py` reproduces §2–§4 analytically. The first Geant4
photon run should be compared against it *before* any campaign conclusion:
if simulated P_int in the drift gas disagrees with table §2 by more than
~10 %, something is wrong with the gas density or the mixture builder
(cf. the P0.1 mass-vs-volume-fraction bug, which would show up here as
exactly this kind of discrepancy).

## 7. What this changes about the campaign

1. **Argon is promoted to a co-baseline.** The gas scan still runs — the
   neon mixtures have real advantages in max-gain headroom and discharge
   probability (`GAS_NOTES.md` §7.2), and those are worth having — but
   photon rejection is no longer the reason to prefer them, and the
   recommendation memo should say so explicitly.
2. **The conversion-layer budget is promoted from an observable to *the*
   decision variable.** Phase 1's stacked-bar figure now decides the gas
   question, not the sensitivity curves.
3. **The mesh model becomes blocking**, not a caveat (§4.2).
4. **Get the background spectrum.** Whether the band is soft or hard moves
   the Ar/Ne ratio between 5× and 2× and is currently unknown to us (§2).
5. **The discrimination study is still worth doing** — it is what sets the
   charge window, it kills the fluorescence and contained-Compton
   populations efficiently, and it is measurable on the bench with ⁵⁵Fe vs
   ⁹⁰Sr. It is just not a gas-choice argument.

## 8. What would falsify the §4 conclusion

- Geant4 with a proper woven-mesh geometry giving a wall-conversion rate
  ≳10× below the estimate here (would need the escape model to be wrong by
  much more than the factor 2–3 claimed).
- A background spectrum concentrated at ≤50 keV, where the Ar/Ne ratio is
  5.4× and the wall floor is also higher — worth recomputing, the ratio of
  wall to gas is roughly energy-independent so this probably does *not*
  rescue neon, but it should be checked rather than assumed.
- Amp-gap partial gain suppressing the Cu contribution far more than
  assumed, leaving the mesh alone — that halves the floor but does not
  change the conclusion.
- Wall photoelectrons turning out to be rejectable after all, e.g. if their
  charge spectrum is harder than a MIP's by enough to sit above the upper
  cut. §3 suggests they land at ~3–5× MIP, uncomfortably in the middle;
  this is a real possibility and the sim will answer it directly.
