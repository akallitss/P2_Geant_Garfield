# COMSOL mesh unit cell (P2)

Base model: `MeshUnitCellIon_PressureStudy.mph` (a colleague's woven-mesh unit
cell, COMSOL 6.3, used with permission). Not in the repo (850 MB); a copy is
on EOS: `/eos/user/a/akallits/p2_sim/comsol/`.

Runs on lxplus with COMSOL 6.4 from CVMFS
(`/cvmfs/projects.cern.ch/engtools/comsol/comsol64/multiphysics/bin/comsol`)
and the CERN licence server (`licencomsol.cern.ch:1718`; needs COMSOL +
AC/DC, both available as batch licences).

| file | what |
|---|---|
| `InspectMeshCell.java` | prints parameters, studies, physics, datasets, export nodes |
| `SolveP2MeshCell.java` | sets P2 parameters, solves Study 1, re-runs the Garfield exports |
| `run_mesh_cell.sh`, `mesh_cell.sub` | HTCondor wrapper (8 cores, 24 GB) |

    comsol compile SolveP2MeshCell.java
    condor_submit mesh_cell.sub

P2 settings (SPS July 2026; mesh and drift polarised, pads grounded):
pitch 2×68 µm (measured 18 µm wire / 50 µm opening; the weave repeats every
two wires), rWire 9 µm, meshHeight 150 µm, VAnode 0, VMesh −450 V, drift
field 625 V/cm: h4mm = height 4 mm, VCathode −700 V; h2mm = 2 mm, −575 V.

Output: `/eos/user/a/akallits/p2_sim/comsol/p2_sps/<case>/mesh.mphtxt` +
`potential.txt` (Garfield `ComponentComsol`). Materials: gas ε=1 (domains
1–3), mesh/insulator ε=4 (domains 4–11).

Gotchas
- COMSOL's default security settings forbid Java file I/O and reading
  system properties in a batch class: make folders in the wrapper, print to
  the log, hard-code paths.
- Never let COMSOL save the model to AFS (quota); run with `-nosave`,
  `-tmpdir`/`-recoverydir` on the job scratch.
- In the base model the weighting-potential exports (data2–4, V2/V4/V5) are
  disabled and Study 1 solves only `V`; enable them if the cell's own
  weighting fields are ever needed.
