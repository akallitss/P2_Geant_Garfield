# P2 detector design data

Copied 2026-08-05 from
`/media/dylan/LaCie/Extras/Physics/Post-Doc-Saclay/P2-detector-design/`.
Dates below are the ones written *inside* the files, which do not always match
the directory names on the drive (`Version-Nov-25` holds March 2024 output;
`Version_Apr26` holds both January and May 2026 output).

Measurements taken from these files are in [`../docs/P2_GEOMETRY.md`](../docs/P2_GEOMETRY.md).

## `gerbers/P2_BASKET_Apr26/` — current production readout ★

KiCad 9.0.2, generated 2026-01-20, project `P2_BASKET`.
From `wedge_det_design/Version_Apr26/FABRICATION_P2/`.

```
Gerber/     F_Cu, B_Cu, Edge_Cuts, F_ENIG_Finition, User_Drawings
NCDrill/    P2_BASKET.drl + drl_map (dxf, pdf)
Artworks/   per-layer PDFs
Stack_Up_P2.txt    2-layer: 18 µm Cu / 200 µm FR4 core / 18 µm Cu
TOP.jpg BOT.jpg    photos of the finished board
```

This is the primary geometry source. `Edge_Cuts.gbr` also contains the 610 × 610
mm fabrication panel and the 2 mm routing channel — neither is part of the
detector.

## `gerbers/bulk_masks_CERN/` — bulk Micromegas masks ★

Ucamco UcamX, 2026-05-13. From `wedge_det_design/Version_Apr26/Bulk_CERN/`.

`P2_Mask1.gbr` / `P2_Mask2.gbr` (+ PDFs) are the CERN bulk soldermask / pillar
masks. `P2_Mask2` carries the 41 366 pillars and therefore defines the active
area. `3717-A-Alexandra Kallitsopoulou - P2 Micromegas bulk.pdf` is the bulk
run request.

`older_revisions/` holds the V1 and V2 revisions of the same two masks
(`P2_BASKET-Mask_M1/M2_V1|V2.gbr`), from `Version_Apr26/V1/` and `V2/`.

## `gerbers/DummySector0_Nov25/` — earlier wedge

KiCad 7.0.7, dated 2024-03-26, rev 0.30. From
`wedge_det_design/Version-Nov-25/GERBERS/`. F_Cu, B_Cu, Edge_Cuts,
User_Drawings, drill and a `.gbrjob`. Kept for diffing against the current
design; not the production geometry.

## `gerbers/RKP2_10x10/` — prototype readout

KiCad 6.0.1, 2022-10-21. From `10x10-design/RKP2/`. Gerbers, both drill files
and the `.kicad_pcb`. Relevant only if a small-prototype simulation is wanted.

Despite the `10x10-design` directory name, `Edge_Cuts` is a 300 × 300 mm square
centred on the origin; the copper only covers ~256 × 256 mm of it. The "10×10"
presumably refers to the active area, which was not separated out here.

## `gerbers/K59V_Fx2Mec8/` — interconnect board

KiCad 8.0.0, 2024-05-02. From `wedge_det_design/K59V/_ EXPORT/`. 14 gerber
layers (including 4 inner signal/GND layers), NCDrill and the BOM. Sits off the
active area — material budget only.

Left on the LaCie drive: the KiCad project and schematic, the DXF exports, the
12.7 MB `Fx2Mec8.step`, six dated backup zips, the Würth manufacturer copy, and
the Hirose FX20-140 / Samtec MEC8 vendor 3D models.

## `mapping/` — channel mapping

From `wedge_det_design/Maarten-Irakli-mapping/`. Two revisions of the same
1280-pad map:

- `connector_0..9.txt` — 7 columns, has `AbsIndex`, ordered outer radius inward
- `Mapping/connector_0..9.txt` — 9 columns, adds `X`, `Y`, `PadName`, ordered
  inner radius outward, plus `MappingFx2Mec.xlsx`

They agree on 1201 of 1280 pad positions. **Which one the DAQ uses is an open
question** — see the handoff.

`README` (French) documents the index conventions. `module_equalized.gbr` is a
non-CAD gerber of raw draw commands in mm that goes with the mapping, with a
`.png` preview.

> `Phi` and `DeltaPhi` in these files are in **radians**, not degrees.

## `mechanical/` — drift frame

From `DRIFT_P2_261124/DRIFT_P2_261124/`.

- `P2_Frame_V1_3mm.stp` — 3 mm drift variant
- `P2_Frame_V2_4mm.stp` — 4 mm drift variant

The matching `.stl` meshes (~280 MB each, very finely tessellated) were left on
the drive.

## `docs/`

`FX20-60P-0.5SV15...pdf` (Hirose connector catalog) and `MariamInventaire_P2.pptx`.

## Not copied

`detectors-labo-bulk-info/` on the drive holds no design files — it is bulk
fabrication travelers (Excel *Fiche de préparation* sheets for mesh and bulk
runs, 2024–2025).
