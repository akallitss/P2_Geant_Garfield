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
     radius-dependent. **The ~1.1-pad number now carries more weight than
     when it was written** (2026-08-06): it is also the basis for running
     with VMM neighbor logic off and for treating clustering as a minority
     effect (plan §2, `vmm/README.md` §3.1). If the real angular
     distribution puts a significant fraction of tracks at 30°+, both of
     those conclusions become radius-dependent too — another reason the
     phase-space file below is the highest-value ask.
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

### 🟡 Received 2026-10-06: `BkgHistograms.root` (Matthieu, P2Sim detTest, MMP2_v4)

`../data/Simulations/background_mattieu_data/BkgHistograms.root`. Per MM
layer i = 0, 1, 2: `hEnergy_i`, `hEnergy_response_i` (0–1 MeV, 5 keV bins),
`h2EnergyTheta_i` (x = E [MeV], y = θ [deg], z = rate [MHz]), `h_map_xy_i`
(TH2Poly), plus `h_rate_MM` / `hr_rate_MM` (totals). Matthieu: a photon is
counted **once, when it first enters the drift gap**. Mostly that is the
entrance face, but photons from behind are recorded at the exit face. The
weights are a plain rate in MHz, with no further normalisation.

What is in it (checked 2026-10-06):

| | MM0 | MM1 | MM2 |
|---|---|---|---|
| photons entering the drift gap [MHz, whole layer] | 47 120 | 38 928 | 32 824 |
| response-weighted (`hr_rate_MM`) [MHz] | 101.8 | 42.0 | 35.2 |
| MC entries (each ≈ 9.33 MHz) | 5 052 | 4 177 | 3 522 |
| mean E | 100 keV | 107 keV | 109 keV |
| E < 20 keV | 3.6 % | 0.3 % | 0.2 % |
| θ < 90° (in the 2D underflow) | **5.2 %** | 3.7 % | 3.1 % |
| E > 1 MeV (overflow) | 0.4 % | 0.5 % | 0.5 % |
| incidence (180° − θ): median / within 30° | 17.5° / 78 % | 16.5° / 82 % | 16.5° / 85 % |

- The XY map covers the full 360° (all wedges, flat in φ to ±15 %) over
  r 121–589 mm, roughly flat per unit area. Its integral is 1000 ×
  `hr_rate_MM`, so it is the **response-weighted rate in kHz**, not the
  photon flux (to confirm).
- 5 000 entries spread over 200 × 90 bins is sparse: fine for 1D marginals,
  thin for a correlated E–θ sampling.

**What to ask Matthieu:**

1. **θ convention**: measured from +z (beam axis)? θ in 90–180° means the
   photons travel upstream. Which face of the MM (drift-cathode or PCB
   side) faces +z in detTest? That decides whether a photon crosses the
   cathode first or the PCB first in our model.
2. **Extend θ to 0–180°**: 3–5 % of the flux (the "from behind" photons)
   is in the underflow and has lost its angle.
3. **Record the flux at the outer window instead of the drift gap**, if
   possible. The current spectrum already includes his windows, cathode,
   copper and FR4 (which differ from ours, `figures/p2_stack_vs_p2sim.png`).
   Feeding it through our geometry would attenuate it twice.
4. **Log or finer energy binning below 50 keV, and range up to ~10 MeV.**
   Fluorescence lines and the soft part dominate the conversion response.
5. **More statistics**, or the raw per-photon list (E, θ, φ, x, y, entry
   face) instead of histograms. A list is what our generator can sample
   from directly.
6. Beam current / conditions the MHz refer to, and whether the map is
   kHz response-weighted.

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

---

## Addendum (2026-08-06, from the MX17 model work): bulk pillars not modelled

The P2 wedge model has **no bulk pillars** — the amplification gap is pure
gas. MX17's bulk gerber (`3498A_bulk.gbr`) shows the real pattern for that
detector: Ø 0.6 mm pillars on a regular 4.68 mm grid, i.e. ~1.3 % of the amp
gap is polyimide, not gas, with locally dead amplification spots. The MX17
model now places them (`MX17_Geant/shared/MX17ModuleGeometry.hh`); if the P2
production uses a similar bulk process, the same pattern likely applies.

**Ask for:** the P2 bulk pillar mask (diameter + pitch), or confirm the
MX17-like 0.6 mm / 4.68 mm pattern, then add pillars as daughters of the
`AmpGas` volume as MX17 does.

### ⬜ TO-DO before any production run: get the Dynamask pillar map

*(Dylan, 2026-08-06.)* **Get the Dynamask pillar map from
`P2_Basket_Analysis` (or wherever it lives) before running.** This gates
production runs — the pillar layout has to be in the geometry, not
retrofitted afterwards, because it changes both the active-gas fraction and
where the dead spots are.

What we already have in-repo, and what is still missing:

| | status |
|---|---|
| Pillar **pattern** from the CERN bulk mask `design/gerbers/bulk_masks_CERN/P2_Mask2.gbr` | ✅ Ø **0.5 mm**, pitch **2.000 mm** exact, **41 366** pillars, **4.8 %** of the amp gap (`P2_GEOMETRY.md` §2) |
| Pillar **material** | ✅ **Dynamask** photoimageable dry film, ρ ≈ 1.2–1.4 g/cm³, ε_r ≈ 3.9 — *not* kapton/FR4 (`HANDOFF_MX17_RESPONSE.md` §2.4) |
| Pillar **map actually used by the analysis** (as-fabricated positions / dead-channel-adjacent pillars, in whatever form `P2_Basket_Analysis` consumes) | ⬜ **NEEDED — this item** |
| Pillars in the Geant4 geometry | ✅ 2026-10-06: `--pillars saclay` (V1 mask, default) / `cern` / `none`; see `P2_MODEL.md`. Stage B dead spots still open |

Why the mask alone may not be enough, i.e. why to go get the analysis's map:
the gerber gives the *design* pattern, while the analysis presumably carries
the map it actually uses to mask pad regions — and if the two disagree
(revision, origin offset, masked/edge pillars), the discrepancy is exactly
the kind of thing that silently biases per-pad efficiency. Getting the same
file the data analysis uses also keeps sim and data on one definition, the
same argument as the pad-mapping revision (`SIM_CAMPAIGN_PLAN.md` §10.6).

What it feeds: pillar dead spots in Stage B (a drifting electron landing on a
pillar is lost, and the amplification field is locally absent), the 4.8 %
effective-gas correction in the amp gap, and the Dynamask material budget.
Filed as **P0.16** in the campaign plan.

---

## ⬜ Readout-copper zoning, real pad/strip pattern, resistive layer

*(Dylan / MX17_Geant, 2026-08-06.)* Three accuracy upgrades landed in
MX17_Geant that apply just as well to the P2 wedge. Ordered by effort/benefit.

**Ticketed 2026-08-07** as **P0.18** (Cu zoning), **P0.19** (real pad
artwork) and **P0.20** (resistive layer) in `SIM_CAMPAIGN_PLAN.md` §3 — until
then these were written down here but were not on anyone's task list, which
is how item 1 (an actual error in the current model, not an input request)
sat unaddressed. This section stays as the detailed write-up; the plan
carries the status.

### 1. Zone the Cu coverage — cheap, and P2 currently has the error

`SimConfig.hh` documents `p2_fcu_coverage = 0.983` / `p2_bcu_coverage = 0.174`
as measured **over the active area**, but `EffCu` then applies them to a
copper slab spanning the whole wedge board. So the model puts near-solid
copper out in the fan-out/periphery, where the real coverage is far lower.

MX17 had the mirror-image of this bug (one *board-wide* average applied
everywhere, giving 26 % too little copper in the active area and ~10× too
much outside). The fix costs two volumes and no CPU: build each Cu layer as
an active-area zone plus a homogenized remainder, each at its own measured
coverage. See `BuildReadoutZone` in `MX17_Geant/shared/MX17ModuleGeometry.hh`.
**Needed input:** re-run the P2 coverage script reporting active *and*
outside-active fractions separately (it currently reports one number).

### 2. Real pad/strip pattern — ~3 % CPU if the artwork is regular

If the P2 readout artwork is a regular grid (MX17's is exactly 512 × 512 on a
0.78 mm pitch), it can be built as real geometry with two nested
`G4PVReplica` levels: the whole 786432-feature MX17 pattern costs **15 extra
volumes**, not 786432 placements, ~+3 % CPU and no extra memory. Only worth it
for backscatter off the board — the signal copper is downstream of both gas
gaps. **Needed input:** confirm the P2 pad/strip artwork is periodic.

### 3. Resistive layer — P2 models it as a solid slab

`DetectorConstruction.cc` places `ResistivePaste` as a **100 µm slab of
full-density (1.4 g/cm³) paste** covering the whole wedge, with no coverage
scaling and no structure. MX17's is now 515 discrete 550 µm ESL strips on a
0.8 mm pitch inside a gas-filled envelope, so the inter-strip grooves are real
chamber gas. This matters more than the copper pattern: the resist is the
first solid the avalanche region sees and it is 100 µm thick, ~5× a P2 Cu
layer. **Needed input:** is the P2 resistive layer strips, a uniform DLC/paste
sheet, or pads? If it is uniform, the current solid slab is right and only the
thickness/density need confirming — but that should be stated, not assumed.

### Do NOT port the MX17 rasterizer fix — but P2's script has its own bug

MX17's `analyze_cu_coverage.py` had an endpoint-inclusive PIL bug that
inflated every feature by one pixel per dimension (+16 % on 0.68 mm pads at
0.05 mm/px), so its published coverages were ~13 % high until 2026-08-06.
P2's equivalent script uses shapely (exact vector geometry) and is **not**
affected by *that* bug. Flagged so nobody "fixes" a bug that is not there.

**Corrected 2026-08-07 — the original wording ("is not affected", full stop)
was too strong and told the next reader not to look.** P2's script has a
different defect with a similar effect. `ApertureDef.size` in
`scripts/gerber/gerber_outline.py` returns `max(self.params)`, which is right
for `C`/`R`/`O` but wrong for KiCad's **`RotRect` macro apertures**, whose
third parameter is a **rotation angle in degrees**, not a dimension:

```
%ADD19RotRect,0.300000X1.800000X302.763000*%     <- 0.3 x 1.8 mm, rotated 302.763 deg
```

`copper_union` sends any non-`R` flash down the `Point(...).buffer(ap.size/2)`
branch, so each of these paints a disc of radius ~151 mm. `P2_BASKET-F_Cu.gbr`
declares 21 such apertures (verified 2026-08-07), flashed at the ten connector
footprints. The script also counts apertures tagged `NonConductor` and
`Profile` as copper.

⇒ **The 0.937 board-level figure was inflated by 19 %.** ✅ **Fixed and
cross-checked 2026-08-07 by two independent methods**: the rasterized
extractor in `extract_readout_pattern.py`, and the shapely (exact vector
geometry) path in `analyze_cu_coverage.py` after its aperture handling was
repaired. They agree to 0.05–0.3 %:

| | shapely (exact) | rasterized | old (buggy) |
|---|---|---|---|
| F.Cu board-wide | **0.7866** | 0.787 | 0.937 (+19.1 %) |
| B.Cu board-wide | **0.1803** | 0.180 | — |
| F.Cu active area | **0.9749** | 0.978 | 0.983 (+0.8 %) |
| B.Cu active area | **0.1743** | 0.171 | 0.174 (−0.2 %) |

**The important nuance, and why nothing published needs retracting:** the
`RotRect` bug painted ~151 mm discs at the ten connector footprints, which
sit *outside* the active area. So it wrecked the board-wide number (+19 %)
while leaving the active-area number nearly right (+0.8 %). Every result so
far came from a beam inside the active area, which is why the conversion
budget survived the correction unchanged. And `p2_bcu_coverage = 0.174` was
correct all along, to 0.2 % — only `P2_GEOMETRY.md`'s prose describing B.Cu
as a solid ground plane was wrong.

Also correct in passing: `P2_GEOMETRY.md` §3 describes B_Cu as "a single solid
ground plane". It is not — it is stroked signal traces plus via pads, which is
what the coverage script's own docstring says ("signal lines, not a plane") and
what `p2_bcu_coverage = 0.174` already encodes. The model input looks right;
only the prose is wrong.
