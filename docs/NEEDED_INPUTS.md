# Needed external inputs — the register of what we are guessing

**Written:** 2026-08-05. One place for everything the simulation currently
*assumes* because nobody has given us the number, ordered by how much the
answer would move a campaign conclusion. Each item says what we assume now,
why it matters quantitatively, and what specifically to ask for.

Companion docs: `P2_MODEL.md` (geometry assumptions — Alexandra's review),
`testbeam/TB_CONDITIONS.md` (test-beam as-run conditions),
`research/PHOTON_DISCRIMINATION_NOTES.md` (photon physics),
`SIM_CAMPAIGN_PLAN.md` (the campaign these feed).

Legend: ⬜ needed · 🟡 partially known · ✅ have it.

---

## 1. ⬜ Spatial and angular distribution of the signal electrons

**Raised by Dylan, 2026-08-05: "I don't know the spatial distribution of the
electrons nor the X-ray background, and I assume that matters."** It does,
and more than the campaign currently accounts for.

### What we assume today

A **pencil beam at one point**, `--gun-x 307.4 --gun-y 177.5`
(r = 355 mm, φ = 30°, mid-active-area), travelling along +z, i.e. normal
incidence, with the incidence angle treated as a free scan parameter
0°–40° (`P2_EXPERIMENT.md` §3 bottom line). Position dependence is
currently "deferred to a small dedicated scan" (plan §4.1).

### Why it matters

1. **Rate per pad is a local quantity, and we only have an average.** The
   pads are area-equalized at ≈ 1.32 cm², and the quoted **15–20 kHz/pad**
   is presumably an area average over the 1689 cm² active area. If the
   electron flux falls like 1/r² over the active span r = 120 → 590 mm,
   inner-to-outer flux differs by **(590/120)² ≈ 24×**, so the hottest pads
   could run at ~5–10× the mean. Everything rate-driven keys off the
   *worst* pad, not the average: streaming bandwidth, VMM occupancy and
   pile-up, baseline shift from the ion tail, and — most importantly — the
   accidental-coincidence arithmetic of plan §5 step 5, where the fake rate
   goes as the *square* (2-fold) or *cube* (3-fold) of the per-layer hit
   probability. A factor 5 in local rate is a factor 25–125 in accidental
   triples. **This alone can change the headline rejection number by more
   than the entire gas choice does.**

2. **Incidence angle is correlated with radius, and we treat it as a free
   global parameter.** With no B field at the wheels (our working
   assumption) and a target at distance L, a straight track landing at
   radius r arrives at θ = atan(r/L) — so θ varies *across the wedge*,
   not run to run. Consequences:
   - path length in the gap goes as 1/cos θ (1.00 → 1.16 over 0°–30°),
     shifting the Landau MPV and the zero-cluster floor with radius;
   - **pad multiplicity is strongly angle-dependent**, and the plan's
     working assumption that a track hits **~1.1 pads at 10°** — which is
     why topology was demoted to a secondary discriminant (plan §5 step 3)
     — does not survive at 30–40°, where a 3 mm gap smears the track over
     ~1.7–2.5 mm transversely. If the outer radius sits at 30°+, topology
     comes back as a real discriminant *there* and the ROC is
     radius-dependent.
   - Also worth checking: with a point target at fixed L, a 30°–40°
     acceptance maps to a **narrow annulus**, not the full r = 120–590 mm
     span (at L = 500 mm, 30°–40° → r = 289–419 mm). The full active area
     is only illuminated once the 60 cm target length, the field, and the
     θ = 140°–150° acceptance are folded in. **We may be simulating an aim
     point that is barely illuminated, or missing that half the wedge sees
     almost nothing.**

3. **Azimuthal uniformity is assumed and probably fine** (the wedge is one
   of six in a wheel, and the physics is azimuthally symmetric) — but the
   wedge has an asymmetric outline (the y ≤ 540 mm top cut) so the corner
   regions are worth a sanity check.

### What to ask for

The clean deliverable is **not** a set of summary numbers but a
**phase-space file at the wheel plane**: one row per particle crossing,
with `(x, y, z, px, py, pz, E, PDG, weight)`, from whatever Geant4/simulation
the collaboration already runs for the P2/BASKET acceptance. That makes our
wedge simulation directly normalizable and removes every assumption above
at once. Failing that, ask for:

- dN/dA(r, φ) for signal electrons at the wheel plane, and the same for the
  photon background (§2);
- the joint distribution of incidence angle with radius, or just the wheel
  z-position and whether there is field at the wheels — from which we can
  compute it;
- whether 15–20 kHz/pad is an average or a peak, and over what area.

**Follow-up task**: this needs a new gun mode in the simulation —
`--phase-space <file>` sampling position + direction + energy from such a
file, alongside the existing pencil beam. Filed as **P0.15** in the campaign
plan. Cheap to write; useless until the file exists.

## 2. ⬜ Spectrum and spatial distribution of the X-ray background

### What we assume today

A **mono-energetic pencil beam** at each of the plan §4.2 scan energies,
normal incidence, at the same single aim point. The band is known only as a
label: "50–150 keV", now refined by colleagues to **50–100 keV with a
higher-energy tail** (no shape).

### Why it matters

1. **The spectrum shape decides the Ar-vs-Ne ratio.** Interaction
   probability ratio Ar/Ne is 5.4× at 50 keV and 2.1× at 150 keV
   (`research/PHOTON_DISCRIMINATION_NOTES.md` §2). A soft spectrum and a
   hard one give materially different answers to the gas question. We are
   currently reporting per-energy and cannot fold.
2. **The photon and electron spatial distributions are probably
   *different*, and their ratio is what matters.** Photons are unaffected
   by any residual field and map straight from the target; electrons at
   30–40° may be field-steered. So the quoted **10³ photon/electron ratio
   is an average that could be badly wrong locally** — if photons
   concentrate at small radius (near the beam axis) while electrons land
   further out, the local signal-to-background varies by an order of
   magnitude across the wedge. The rejection requirement is set by the
   worst region, not the mean.
3. **Incidence angle changes the wall conversion rate**, which is the
   dominant fake source (§4 of the photon note): path length through the
   mesh and the copper pad plane scales as 1/cos θ, so a 40° background
   converts ~30 % more often than a normal-incidence one. It also changes
   *which* gap the escaping photoelectron enters, which sets whether it
   gets full or partial gain.
4. **The high-energy tail needs bounding.** Above ~200 keV the photon
   barely interacts (P ≈ 5×10⁻⁵ in 3 mm of argon) but Compton electrons
   get energetic enough to look exactly like signal. If the tail carries
   appreciable flux it may dominate the irreducible fakes.

### What to ask for

- The generated photon spectrum at the wheel plane (histogram or the
  phase-space file of §1 — same file, both species).
- Its origin breakdown (target bremsstrahlung vs beamline vs collimator),
  because that determines whether a local shield is even geometrically
  possible — the shield trade study is currently sidelined (plan §8).
- Whether "50–150 keV" was a full-spectrum statement or the range of a
  dominant peak.

## 3. 🟡 Geometry — Alexandra's review

Tracked in `P2_MODEL.md` §"Still to confirm"; production runs are gated on
it. Open: mylar thicknesses/aluminisation, back-frame depth, frame
opening/material, overpressure/sag, mesh spec, which pad-mapping revision
the DAQ uses. Not repeated here.

## 4. 🟡 Test-beam as-run conditions

Template at `testbeam/TB_CONDITIONS.md`, handed to Alexandra. Everything in
`TESTBEAM_PLAN.md` blocks on it.

## 5. ⬜ VMM3a noise at the real pad capacitance

No published ENC curve covers 75–170 pF at our gain and peaking-time
settings (`vmm/README.md` §4, §6). This is the make-or-break number for the
short-peaking-time discrimination handle
(`research/TIMING_PSD_NOTES.md` §2). Measurable on the bench: pedestal σ vs
peaking time on real pads.

## 6. ⬜ Gas gain curves for neon mixtures in a 128–150 µm bulk Micromegas

None published for any Ne mixture (`research/GAS_NOTES.md` §7.2), and no
Penning transfer coefficient for Ne/iC₄H₁₀, Ne/CH₄ or any ternary. Mitigated
by bracketing and by reporting in primary-electron equivalents; resolved
properly only by a ⁵⁵Fe/⁹⁰Sr gain scan on the prototype.

---

## Priority

If only one thing can be asked for, ask for **the phase-space file at the
wheel plane covering both electrons and background photons** (§1 + §2). It
closes the two largest open axes at once, converts the campaign from
"scan everything and hope" to "weight by the real distribution", and it
almost certainly already exists inside the collaboration's acceptance
simulation.
