# The gas-mixture fix (P0.1) — two bugs, and what they invalidated

**Written:** 2026-08-05. Status: **fixed, verified in Geant4, and the
affected results re-measured.** Code: `DetectorConstruction.cc`
(`idealRho` + `makeMixV`). Campaign context: `../SIM_CAMPAIGN_PLAN.md` §3
task P0.1.

Every simulated gas in this repo was wrong in two independent ways until
2026-08-05. Both are now fixed by construction rather than by correcting
the numbers, so neither can silently come back.

---

## 1. Bug A — volume fractions passed as mass fractions

`G4Material::AddMaterial(mat, fraction)` takes a **mass** fraction. The
builder handed it the **volume** fraction, which is how gas mixtures are
specified everywhere else (gas systems, literature, our own plan).

For ideal gases at a common temperature and pressure the conversion is

```
    f_mass,i = f_vol,i · ρ_i / Σ_j (f_vol,j · ρ_j)
```

so the error scales with how much the component densities differ from the
carrier's. Isobutane is 2.9× denser than argon and 5.8× denser than neon:

| mixture | intended (volume) | what was simulated (mass) |
|---|---|---|
| Ar/Iso 95/5 | 95 / 5 | **92.9 / 7.1** |
| Ne/Iso 95/5 | 95 / 5 | **86.8 / 13.2** |
| Ar/CO₂ 70/30 | 70 / 30 | 67.9 / 32.1 |
| He/Ethane 96.5/3.5 | 96.5 / 3.5 | 78.5 / 21.5 |

The neon case is the serious one: **2.6× the intended isobutane content by
mass.** Since the whole campaign exists to compare argon against neon, and
the error is 2.4 % in argon but 13 % in neon, it biased exactly the
comparison the campaign is for. He/Ethane is worse still (6× the intended
ethane) but that mixture is inherited MX17 and not on the campaign list.

## 2. Bug B — densities were STP values in a 20 °C material

Independently, every component density was hard-coded at the **0 °C** value
while the `G4Material` declared 293.15 K:

| gas | hard-coded | ideal gas at 20 °C | ratio |
|---|---|---|---|
| Ar | 1.7820 | 1.6607 | 1.073 |
| Ne | 0.8999 | 0.8389 | 1.073 |
| He | 0.1786 | 0.1664 | 1.073 |
| CO₂ | 1.9770 | 1.8295 | 1.081 |
| CF₄ | 3.7200 | 3.6584 | 1.017 |
| iC₄H₁₀ | 2.6700 | 2.4162 | 1.105 |
| C₂H₆ | 1.3560 | 1.2500 | 1.085 |

(mg/cm³.) The uniform 1.073 factor is exactly 293.15/273.15 — these are
textbook STP densities. **Every gas was ~7 % too dense**, which propagates
directly into dE/dx, primary-ionization counts, and photon interaction
probability. It also explains a discrepancy noticed earlier: the analytic
budget in `PHOTON_DISCRIMINATION_NOTES.md` used 20 °C densities and so
disagreed with the simulation's gas conversion rate by ~7 % for reasons
that had nothing to do with physics.

The isobutane value (2.67) is not even the 0 °C ideal figure (2.593); it
looks like a real-gas density at ~0 °C. The plan's P0.1 already flagged it
as suspicious.

Combined, the mixture densities were off by **−7.0 % (Ar/Iso)** and
**−10.3 % (Ne/Iso)** — neon worse because both bugs pushed the same way.

## 3. The fix

Both are now derived rather than tabulated:

- `idealRho(molarMass)` computes each pure-gas density from the molar mass
  at the declared 20 °C / 1 atm. Changing the operating temperature or
  pressure now changes the densities consistently, and there is no second
  place for a stale number to hide.
- `makeMixV(name, {{mat, volumeFraction}, ...})` takes **volume** fractions,
  derives the mixture density as Σ fᵢρᵢ and the mass fractions as
  fᵢρᵢ/ρ_mix, and throws if the volume fractions do not sum to 1. The
  density is no longer computed separately at each call site (which was the
  structural reason the two numbers could disagree in the first place).

Ideal gas is good to <1 % for the noble gases, CO₂ and CF₄ at 1 atm. The
worst real-gas deviation is isobutane, ~2 % denser than ideal at 20 °C; at
≤20 % quencher that is <0.4 % on any mixture here, well inside the P/T
control of a real gas system. Documented in the code rather than corrected
with an invented compressibility factor.

**Verified in Geant4**, not just on paper — `printMaterial` reports:

```
Material: ArIso  density: 1.698 mg/cm3   ElmMassFraction: Ar 92.89 %
          ElmAbundance: Ar 57.58 %, C 12.12 %, H 30.30 %
Material: NeIso  density: 0.918 mg/cm3   ElmMassFraction: Ne 86.84 %
          ElmAbundance: Ne 57.58 %, C 12.12 %, H 30.30 %
```

The atom fractions are the check that matters: 95 Ar + 5 C₄H₁₀ per 100
volumes is 95 : 20 : 50 atoms = 57.6 : 12.1 : 30.3 %. Exact.

## 4. What this invalidated

Everything simulated before 2026-08-05 night, which is one day's worth of
first runs. Specifically:

- The **first photon run numbers** in `PHOTON_DISCRIMINATION_NOTES.md`
  §4.3a were taken with argon 7 % too dense and neon 10 % too dense, and
  with neon carrying 2.6× the intended isobutane.
- The **Ar-vs-Ne comparison** (347 vs 301 events per 10⁵ γ) was the result
  most at risk, since the two gases were wrong by different amounts (−7 %
  vs −10 % in density, 2.4 % vs 13 % in composition). In the event it
  barely moved — 13 % → (12 ± 3) % — because the gas-dependent part of the
  rate is only ~10 % of the total to begin with.
- The **benchmarks** in `BENCHMARKS.md` are unaffected in any way that
  matters: runtime is dominated by geometry navigation, and a 7 % density
  change moves the interaction rate, not the cost per event.

**Corrected numbers are in `PHOTON_DISCRIMINATION_NOTES.md` §4.3a**
(re-measured at 4×10⁵ γ per gas, 4× the original statistics). What moved:

| quantity | pre-fix | post-fix |
|---|---|---|
| photon-induced hit rate, Ar/Iso | 3.47×10⁻³ | 3.53×10⁻³ |
| photon-induced hit rate, Ne/Iso | 3.01×10⁻³ | 3.11×10⁻³ |
| Ar → Ne reduction | 13 % | **(12 ± 3) %** |
| Micromesh conversion rate | 1.98×10⁻³ | 2.10×10⁻³ |
| DriftGas conversion rate (Ar) | 3.05×10⁻⁴ | 3.10×10⁻⁴ |

Almost nothing, and for a reason worth stating: **the dominant fake source
is the walls, whose conversion rate does not depend on the gas at all.** The
mesh and pad-copper rates should have been identical before and after, and
they are, within counting statistics. Only the small direct-gas component
moved, and the two bugs partly cancelled in it (denser gas ⇒ more
conversions; the density error was the larger effect and pushed both gases
the same way).

The qualitative conclusion — walls dominate, argon and neon are close — was
never at risk from this. But the *number* people will quote (the Ar-vs-Ne
gap) was, and it is now measured on a correctly-constituted gas.

## 5. Still open

- **The campaign gas list is not implemented** (plan P0.2). Only the
  inherited eight mixtures exist; the §7.4 candidates (Ne/CO₂ 90/10 = N10,
  Ne/Iso 90/10, the ternaries) still need adding. With `makeMixV` this is
  now a one-line entry each, specified the way the plan quotes them.
- **W values are still per-mixture constants** in `SteppingAction`, not
  derived from composition, and they do not account for Penning transfer.
  That is a separate bias on Ar-vs-Ne comparisons made on `nPrimary` rather
  than `edep` — see `PHOTON_DISCRIMINATION_NOTES.md` §5.2 and plan §5
  step 3a.
- Pressure is assumed 1 atm everywhere. If the chamber runs at slight
  overpressure (the windows bulge, so it does), the densities scale
  linearly and `idealRho` is where to apply it.
