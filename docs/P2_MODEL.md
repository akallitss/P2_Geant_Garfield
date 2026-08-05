# P2 wedge — as-built Geant4 model

**Written:** 2026-08-05. Companion to [`P2_GEOMETRY.md`](P2_GEOMETRY.md) (what
the fab files say) — this file says what the **simulation actually builds**
(`-m p2`, the default mode) and, critically, **which numbers are guesses**.

Implementation: `DetectorConstruction::ConstructP2()` in
`src/DetectorConstruction.cc`, profile + constants in `include/P2Wedge.hh`,
run-time knobs in `include/SimConfig.hh` (`p2_*`) / `mm_sim` CLI.
Python mirror for plotting (kept in sync by hand): `scripts/model/p2_model.py`.

```bash
python scripts/model/plot_p2_model.py        # regenerates the figures below
```

## Coordinates

World x/y = **gerber coordinates**: the wedge apex = MESA beam axis is at
(0,0), so pad-mapping files and simulation output share one frame with no
offsets. z = 0 at the **front window plane** (top of the gas frame), +z points
downstream (beam direction). The default pencil beam aims at
(x, y) = (307.4, 177.5) mm — mid-active-area, r = 355 mm, φ = 30° — from 5 cm
upstream of the front window; move it with `--gun-x/--gun-y`.

## Stack (front → back)

| z [mm] | volume | material / thickness | source |
|---|---|---|---|
| −10.04 … 0 | `FrontWindow_*` | 40 µm mylar, bulged, sag 10 mm | sheets: MX17-like (user); **sag estimated, see below** |
| 0 … 4.0 | `FrontGas` | chamber gas, 4 mm | window→mesh = 8 mm from frame STEP, minus the 1 mm foil gap below |
| 4.0 | `DriftCathode_Mylar1` | 12 µm mylar | two-foil cathode **confirmed** (2026-08-05); thickness **GUESS** |
| 4.012 … 5.012 | `DriftCathode_Gas` | chamber gas, 1 mm | foils "maybe 1 mm" apart (Alexandra 2026-08-05) |
| 5.012 | `DriftCathode_Mylar2` + `_Al` | 12 µm mylar + 0.1 µm Al on drift side | Al on gas side **confirmed**; thicknesses **GUESS** |
| 5.024 … 8.024 | `DriftGas` *(sensitive)* | 3.0 mm | **confirmed baseline** (2026-08-05); campaign scans 1–4 mm |
| 8.024 | `Micromesh` | woven SS 48 µm opening / 19 µm wire → 38 µm slab, ρ×0.223 | spec **confirmed** (2026-08-05, "48×19", reading re-confirmed by Dylan); effective-density model |
| 8.062 … 8.212 | `AmpGas` *(sensitive)* | 150 µm | **confirmed** (2026-08-05) |
| 8.212 … 8.448 | `PCB_Cu_F` / `PCB_FR4` / `PCB_Cu_B` | 18 µm Cu·**0.983** / 200 µm FR4 / 18 µm Cu·**0.174** | Stack_Up_P2.txt (exact); Cu density scaled by the **gerber-measured area coverage** over the active area (B.Cu = signal lines; radial 0.03→0.26) |
| 8.448 … 9.448 | `BackGas` | 1 mm | = carbon back-frame depth: 1 mm assumed for this wedge (Dylan 2026-08-05; a few are 3 mm, `--back-gap`) |
| 9.448 … 14.49 | `BackWindow_*` | 40 µm mylar, bulged, sag 5 mm | glued to rear face of the carbon frame (**confirmed** 2026-08-05) |
| 0 … 8.21 | `GasFrame` (ring) | plastic (**confirmed**); polycarbonate assumed for the type | footprint from STEP |
| 8.448 … 9.448 | `GasFrameBack` (ring) | carbon fiber, 1.6 g/cm³ | **confirmed carbon**, glued to PCB periphery (photos 2026-08-05) |

Transverse profiles (all in `P2Wedge.hh`):

- **board** (PCB layers): r 95→650 mm, edges +10 mm, top cut y ≤ 540.01 — exact from Edge_Cuts;
- **frame footprint**: r 104.5→605.5, edges +7.5, top 537.5 — from STEP point clusters;
- **frame opening** (= all gas layers, windows, mesh, cathode): r 107→603,
  edges +0, top 535 — **approximate** read of the STEP inner walls.

## Window bulge model

The few-mbar overpressure bloats the outer mylar sheets. Each window is a
**terraced dome**: N=6 stacked gas prisms whose profile shrinks about the
window centroid following a spherical-cap law (σ_k = cos(kπ/2N), height
h_k = H·sin(kπ/2N)), each step capped by a flat mylar annulus so **every
z-parallel path crosses exactly one 40 µm mylar thickness**. Oblique tracks
through terrace risers are the residual approximation error. H = 10 mm front,
5 mm back (front bloats more, per the filling described), settable with
`--bulge-front/--bulge-back`.

The sag magnitude is physically motivated rather than measured: the
overpressure is known only as **~1–10 mbar** (2026-08-05). Hencky's
large-deflection membrane solution, w₀ ≈ 0.662·a·(p·a/(E·t))^⅓, with
half-span a ≈ 220 mm, E(mylar) ≈ 4.5 GPa, t = 40 µm, gives w₀ ≈ 10 mm at
3 mbar; the p^⅓ scaling means the full 1–10 mbar range only spans
w₀ ≈ 7–15 mm, so the default is representative. Revisit if the pressure or
an actual sag measurement turns up.

![cross-section](figures/p2_stack_xsec.png)
![3D views](figures/p2_3d_overview.png)
![exploded](figures/p2_3d_exploded.png)
![board copper](figures/p2_board_copper.png)

## Settled 2026-08-05 (Dylan)

- **Drift gap: 3 mm baseline.** The simulation campaign will scan drift gaps
  from 4 mm (trending) down to 1 mm (`--drift-gap`) — mechanics prefers large,
  signal/noise may prefer small.
- **Amplification gap: 150 µm** (not the 128 µm bulk standard).
- **Gas: keep the current set**, `ArIso` default for now; many gases will be
  simulated later.
- **Overpressure: ~1–10 mbar**, exact value unknown → sag defaults justified
  by the Hencky estimate above.

## Settled 2026-08-05 evening (Alexandra, WhatsApp + photos)

- **Drift cathode = TWO mylar foils** "at maybe 1 mm distance", the one on
  the gas side aluminised facing the drift gas. Modelled as
  mylar / 1 mm gas / mylar+Al; window→mesh kept at the STEP-derived 8 mm, so
  the front gas gap dropped 5 → 4 mm. Foil thicknesses still guesses (12 µm
  each). The outer windows are a *different*, non-aluminised foil.
- **Mesh: woven stainless steel, "48×19 opening and pitch"** — read as
  48 µm opening / 19 µm wire (pitch 67 µm, the standard 45/18-style weave).
  Modelled as a 38 µm (2 wire Ø) slab of steel at fill fraction
  π·d/(4·pitch) ≈ 0.223, matching the woven areal mass. Optical transparency
  (48/67)² ≈ 51 % is *not* geometrically modelled.
- **Back frame = carbon**, a wedge ring glued to the PCB periphery ("no
  carbon in the active area"), back mylar glued to its rear face, then the
  whole thing integrates with the drift side. Normally **1 mm** thick, up to
  3 mm in this production ("nothing was the same") → default 1 mm,
  `--back-gap` to change. Photos in the repo owner's archive (2026-08-05).
- **Top frame:** aluminum was discussed for the top frame at one point
  (never carbon on top), but the **top gas frame is plastic** (Dylan
  follow-up) — polycarbonate assumed for the exact type.
- **B.Cu is signal lines only, not a full ground plane.**

## Settled 2026-08-05 (Dylan follow-up)

- **"Front should be 3 mm" reconciled:** the 3 mm is the frame ledge → drift
  foil distance, i.e. exactly the STEP geometry already built. No change.
- **Mesh reading confirmed:** wire Ø 19 µm, opening 48 µm.
- **Gas between the two cathode foils:** chamber gas (as modelled).
- **Back frame: assume 1 mm** for this wedge (`--back-gap` for the 3 mm ones).
- **Cu coverage measured from the gerbers**
  (`scripts/gerber/analyze_cu_coverage.py`, exact shapely union): over the
  active area **F.Cu 0.983**, **B.Cu 0.174** (radial profile 0.03 inner →
  0.26 outer as the signal lines accumulate toward the outer-arc
  connectors). Both Cu layers are now full-thickness slabs of
  density-scaled copper. A radially-dependent B.Cu density would be the next
  refinement if it ever matters.

## Still to confirm with the collaboration

One-figure summary (send this): `figures/p2_stack_questions.png`, regenerate
with `python scripts/model/plot_p2_model.py --only questions`.

![stack questions](figures/p2_stack_questions.png)

Ordered roughly by impact on the physics:

1. **Foil thicknesses** — outer windows (40 µm assumed, MX17-like) and the
   two drift-cathode foils (12 µm each assumed); the foil separation is
   "maybe 1 mm", modelled as exactly 1 mm.
2. **Overpressure / sag** — pin down the mbar value or measure the sag.
3. **Gas mixture & operating point** — `ArIso` (Ar/iC₄H₁₀ 95/5) default;
   campaign will scan gases anyway.
4. **Pillars** — 4.8 % of the amp gap volume; currently ignored (not folded
   into gas density).

Items already settled by data or by Alexandra/Dylan (2026-08-05): board
outline, active area, pillar pitch/count, PCB stack + measured Cu coverage,
drift gap, amp gap, mesh spec, front gap (STEP, "3 mm" = ledge→drift foil),
two-foil drift cathode with chamber gas between, back-frame anatomy (carbon,
1 mm, glued to PCB periphery), top frame plastic, pad map (modulo the
two-revision discrepancy, still open — handoff question 6).
