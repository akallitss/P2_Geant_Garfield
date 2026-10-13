# Mesh field map and transparency (P0.13)

COMSOL unit cell of the woven mesh (scripts/comsol/, 18 µm wire, 68 µm pitch,
pads 0 V, mesh −450 V, drift −700 V over 4 mm → 632 V/cm drift, 31.5 kV/cm
gap) read into Garfield++ (pinned build 927e5c21) by
`scripts/garfield/mesh_transparency.cc`, run on HTCondor with
`scripts/garfield/transparency.sub` (cluster 12522802). Raw outputs:
`/eos/user/a/akallits/p2_sim/transparency/`. Figures:
`python3 scripts/analysis/mesh_field_maps.py <copy of that dir>`.

Transparency = share of electrons started 150 µm above the mesh that pass
below the wires, microscopic tracking (diffusion included), 5000 per point:

| gas | 4 mm cell | 2 mm cell |
|---|---|---|
| Ar/iso 95/5 | 90.6 ± 0.4 % | 90.6 ± 0.4 % |
| Ar/iso 90/10 | 98.7 ± 0.2 % | 98.6 ± 0.2 % |
| Ne/iso 95/5 | 93.2 ± 0.4 % | 92.5 ± 0.4 % |
| Ne/iso 85/15 | 99.9 % | 99.9 % |
| Ar/CO2/iso 93/5/2 (SPS) | 95.9 ± 0.3 % | 95.4 ± 0.3 % |
| Ar/CF4/iso 88/10/2 (SPS) | 94.6 ± 0.3 % | 94.9 ± 0.3 % |

Field lines alone (RKF, no diffusion) give 100 % in every case. The two cells
agree, so the 2 mm cell is enough for a drift-voltage sweep.

In Ar/CF4/iso the loss is almost all CF4 attachment (5.3 % of electrons,
status −7), not the mesh (0.1 %).

## Drift-field sweep and literature check

`scripts/comsol/make_sweep.sh` solves the 2 mm cell once per drift field
(150–6000 V/cm, mesh −450 V; COMSOL re-meshes per solve, so fields cannot be
superposed node by node), then `scripts/garfield/transparency_sweep.sub` runs
3000 electrons per field and gas. Numbers in `transparency_sweep.csv`,
figure `fig_transp_sweep.png`.

Comparison with Nikolopoulos et al., JINST 6 (2011) P06011 (arXiv:2005.01065;
Ar/CO2 85/15, 18 µm calendered mesh, 63.5 µm pitch, 128 µm gap,
E_amp 41.4 kV/cm; measured points read off fig. 5, ±0.02):

- Measurements are relative (normalised below 600 V/cm), so only the shape vs
  E_drift/E_amp is tested. Our Ar/CO2 85/15 is flat to ratio 0.02 and then
  falls like the data, 3–9 % above it (our mesh is more open, 54 % vs 51 %
  optical; pitch, gap and E_amp also differ).
- Field lines alone start to fall at ratio 0.04 in both (their 1.7 kV/cm).
- Absolute plateau in Ar/CO2 85/15: theirs 98.6 % (calculation scaled by
  1.4 %), ours 98.7 % (1.1–1.4 % is CO2 attachment).
- Bulk Micromegas measurement (arXiv:2602.21238, 18/63 µm mesh) puts the
  optimum ratio for Ar/iso 95/5 at 0.0094. Ours starts falling above ~0.01.

So P2 at SPS (ratio 0.020) is at the edge of the plateau for cold mixtures
but past it for Ar/iso 95/5 and Ne/iso 95/5 (hot electrons, D_T ≈ 430–490
µm/√cm at 632 V/cm), which lose 7–10 %. A drift field of ~300 V/cm would
recover most of that (Ar/iso 95/5: 98.4 %), at the cost of a slower drift.

## Review (2026-10-14) — supersedes the numbers above

A review of the method raised five points; what each turned out to be:

1. *Every non-passing electron counted as a wire hit.* Checked from the
   saved status codes: every loss is a wire hit (−5) or attachment (−7), none
   left through the top. Attachment is now shown separately (it is the whole
   Ar/CF4/iso loss). The top of the volume is now 1 mm above the start: no
   change.
2. *Steps clipping a wire.* The pinned Garfield has the check but it is
   opt-in for field maps (`ComponentFieldMap::EnableInterfaceCrossings`),
   which the first runs did not call. Switching it on exposed a Garfield bug
   on periodic maps (steps across the cell face were tested along a chord
   through the whole cell: 45 % instead of 90 %), fixed locally in
   `scripts/garfield/patches/`. With the fix: < 0.1 % change.
   **Null-collision steps** (`EnableNullCollisionSteps`), the other test,
   matter: by default a free flight between two real collisions uses the field
   at its start; near the Ramsauer minimum (Ar or Ne with 5 % isobutane) the
   free paths exceed the 18 µm wire, so near the mesh this is wrong. With the
   field re-read at every null collision the wire losses halve in those gases
   (Ar/iso 95/5 90.6 → 96.1 %, Ne/iso 95/5 93.0 → 96.0 %); cold mixtures
   move < 2 %. All runs below use null-collision steps.
3. *Geometry.* The colleague's `meshHeight` is the wire centre plane, so the
   mesh bottom was at 135.5 µm; the mesh is calendered in the model (cut flat
   at ±1.5 rCut = ±14.85 µm, 29.7 µm thick, as the 30 µm of Nikolopoulos et
   al.). Re-solved with meshHeight 164.85 µm: mesh bottom on the 150 µm
   pillars, E_amp 31.5 → 28.5 kV/cm. Whether the P2 mesh is calendered, and
   its thickness, is not measured. `docs/P2_MODEL.md` updated for the
   measured 50/68 µm mesh (Geant4 default still 45/63).
4. *Voltages.* Mesh −450 V / drift 625 V/cm are assumed (the beam-test
   voltages are still open in `docs/testbeam/TB_CONDITIONS.md`), and all gases
   are compared at the same voltages, not at equal gain. Read the sweep at the
   field ratio of the real operating point.
5. *Wording.* Stage B still uses transparency 1.0; the 2 mm/4 mm agreement is
   tested at 625 V/cm only.

Also found: the Ar/CO2/iso Townsend coefficients in the Stage B JSON reached
1e33 /cm between 3 and 15 kV/cm (export sampled a regular grid through the
table's gap). `export_gas_json.C` now writes the table's own field points;
JSON regenerated; nothing used alpha yet.

### Final SPS-point numbers (mesh bottom on the pillars, null-collision steps)

`transparency_review.csv`, `fig_transp_review.png`; 2 mm cell, 3000 electrons:

| gas | T | of the loss: attachment |
|---|---|---|
| Ar/iso 95/5 | 94.3 ± 0.4 % | 0 |
| Ar/iso 90/10 | 99.2 ± 0.2 % | 0 |
| Ne/iso 95/5 | 93.8 ± 0.4 % | 0 |
| Ne/iso 85/15 | 99.9 % | 0 |
| Ar/CO2/iso 93/5/2 | 97.5 ± 0.3 % | 0.4 % |
| Ar/CF4/iso 88/10/2 | 95.5 ± 0.4 % | 4.5 % (wires ~0) |

Drift sweep: `fig_transp_sweep_gap150.png`, `transparency_sweep_gap150.csv`.
All gases are flat to E_drift/E_amp ≈ 0.016; at 450 V/cm Ar/iso 95/5 reaches
99.3 %.

### Validation on Nikolopoulos et al. 2011

Their geometry rebuilt in the same COMSOL cell (`PITCH=63.5 MH=142.85 VM=563
TAG=_nik make_sweep.sh`: 18 µm wire, 63.5 µm spacing, calendered ~30 µm,
mesh bottom on 128 µm pillars, E_amp 41.6 kV/cm vs their 41.4) and their gas
(Ar/CO2 85/15), 4000 electrons per point, both stepping modes
(`transparency_sweep_nik.sub`, `fig_transp_nikolopoulos.png`,
`transparency_nikolopoulos.csv`). Normalised like the data (< 600 V/cm), the
rms deviation from the measured points over 600–6000 V/cm is 2.4 % with
null-collision steps (max 3.8 %) and 3.4 % with default steps (max 6.6 %);
their own calculation agreed to 5 %. The earlier 3–9 % offset came from
comparing different geometries. Limitation: Ar/CO2 85/15 is a cold gas and
barely tests the stepping issue; the hot mixtures are not validated on data.

### Why the gases differ (`fig_transp_ek.png`)

Wire losses at the SPS point follow the electron characteristic energy
ε_k = e·D_T/μ = D_T²E/2 (Magboltz, 638 V/cm): Ar/CF4/iso 0.13 eV → 0.0 %,
Ne/iso 85/15 0.23 eV → 0.1 %, Ar/iso 90/10 0.40 eV → 0.8 %, Ar/CO2/iso
0.46 eV → 2.0 %, Ne/iso 95/5 0.60 eV → 6.2 %, Ar/iso 95/5 0.76 eV → 5.7 %.
Noble gases heat the electrons (elastic collisions only below 11.5/16.6 eV),
isobutane cools them (vibrational levels at 0.1–0.4 eV); hot electrons diffuse
off the converging field lines onto the top of the wires.
