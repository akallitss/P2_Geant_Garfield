"""frame.py — the drift frame, read from Stage A's RunMeta.

WHY THIS EXISTS AS A MODULE RATHER THAN THREE CONSTANTS

World z is not drift depth. In the P2 stack the drift gas sits at positive
world z, several mm downstream of the front window, and electrons drift
toward **increasing** z (the mesh is downstream of the gas). A digitizer
that treated the cluster's z as "distance drifted" would be wrong by both an
offset and a sign, and would be wrong *quietly* — the numbers stay finite and
plausible. MX17 hit precisely this and had half their gap sitting at negative
depth (`MX17_Geant/response/digitizer/clusters.py`).

So Stage A records the frame (`RunMeta.meshZ_mm`, `driftEntryZ_mm`,
`ampEntryZ_mm`, `padPlaneZ_mm`, `driftSign`) and Stage B looks it up. Never
hardcode a z here: the drift gap is a campaign scan axis, so these move from
run to run.

    depth_mm = driftSign * (meshZ_mm - z_world_mm)

with depth 0 at the mesh and +gap at the cathode. Clusters in `AmpGas` come
back with **negative** depth, which is correct and is how they are told apart
from drift clusters — they are already past the mesh and are amplified over
only the remaining part of the gap.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import uproot


@dataclass
class DriftFrame:
    """Geometry and provenance of one Stage A run point."""

    mesh_z_mm: float
    drift_entry_z_mm: float
    amp_entry_z_mm: float
    pad_plane_z_mm: float
    drift_sign: float

    drift_gap_mm: float
    amp_gap_um: float
    gas: str
    gas_label: str
    w_value_eV: float

    particle: str
    energy_MeV: float
    thrown: int
    geometry_hash: str
    git_hash: str
    git_dirty: bool
    beam_spread_mm: float

    # ------------------------------------------------------------------
    def depth_mm(self, z_world_mm: np.ndarray) -> np.ndarray:
        """Drift depth: 0 at the mesh, +gap at the cathode, <0 past the mesh."""
        return self.drift_sign * (self.mesh_z_mm - np.asarray(z_world_mm))

    def amp_depth_mm(self, z_world_mm: np.ndarray) -> np.ndarray:
        """How far INTO the amplification gap a cluster sits, from the mesh.

        An electron released at depth d inside the gap is amplified over only
        the remaining (gap - d), so the gain spans orders of magnitude across
        this one volume. OUTPUT_FORMAT.md §4 warns about exactly this: do not
        sum edepAmp and treat it like edepDrift.
        """
        return self.drift_sign * (np.asarray(z_world_mm) - self.amp_entry_z_mm)

    @property
    def amp_gap_mm(self) -> float:
        return self.amp_gap_um * 1e-3


def _s(arr, i=0):
    v = arr[i]
    return v.decode() if isinstance(v, bytes) else str(v)


def read_frame(root_file: Path) -> DriftFrame:
    """Read the frame from a Stage A file's RunMeta tree.

    Raises rather than guessing if RunMeta is absent: a file old enough to
    lack it also predates the drift-frame branches, so there is nothing to
    fall back to that would not be an invented number.
    """
    root_file = Path(root_file)
    with uproot.open(root_file) as f:
        if "RunMeta" not in f:
            raise RuntimeError(
                f"{root_file.name} has no RunMeta tree, so its drift frame is "
                f"unknown. It was produced before 2026-08-07; re-run Stage A. "
                f"Do not assume z is the drift depth — see this module's "
                f"docstring for why that fails silently.")
        m = f["RunMeta"].arrays(library="np")

    if "meshZ_mm" not in m:
        raise RuntimeError(
            f"{root_file.name} has a RunMeta tree but no drift-frame branches "
            f"(meshZ_mm). Produced between the RunMeta and drift-frame "
            f"commits; re-run Stage A.")

    for col, label in (("geometryHash", "geometry"), ("gitHash", "code")):
        vals = {_s(m[col], i) for i in range(len(m[col]))}
        if len(vals) > 1:
            raise RuntimeError(
                f"{root_file.name} merges {len(vals)} different {label} "
                f"versions ({', '.join(sorted(vals))}). These are not one run "
                f"point and must not be digitized together.")

    return DriftFrame(
        mesh_z_mm=float(m["meshZ_mm"][0]),
        drift_entry_z_mm=float(m["driftEntryZ_mm"][0]),
        amp_entry_z_mm=float(m["ampEntryZ_mm"][0]),
        pad_plane_z_mm=float(m["padPlaneZ_mm"][0]),
        drift_sign=float(m["driftSign"][0]),
        drift_gap_mm=float(m["driftGap_mm"][0]),
        amp_gap_um=float(m["ampGap_um"][0]),
        gas=_s(m["gas"]),
        gas_label=_s(m["gasLabel"]),
        w_value_eV=float(m["wValue_eV"][0]),
        particle=_s(m["particle"]),
        energy_MeV=float(m["energy_MeV"][0]),
        thrown=int(m["thrown"].sum()),
        geometry_hash=_s(m["geometryHash"]),
        git_hash=_s(m["gitHash"]),
        git_dirty=bool(m["gitDirty"][0]),
        beam_spread_mm=float(m["beamSpread_mm"][0]),
    )
