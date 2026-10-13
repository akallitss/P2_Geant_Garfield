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
