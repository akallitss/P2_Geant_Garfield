"""clusters.py — read a Stage A ClusterTree into the arrays the digitizer wants.

Only `DriftGas` and `AmpGas` clusters are ionization the response chain should
see. Everything else in the tree (there is none today, but the schema allows
it) is material budget, not signal.

`nPrimary` is the electron count to use, NOT edep/W. SteppingAction already
applied W and carried the sub-W remainder probabilistically; re-deriving here
would double-count the conversion and throw the remainder away. Same note as
MX17's clusters.py, and it applies verbatim.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import uproot

GAS_VOLUMES = ("DriftGas", "AmpGas")

_BRANCHES = ["eventID", "x", "y", "z", "nPrimary", "time", "volume"]


def read_clusters(root_file: Path, max_events: int | None = None) -> dict:
    """Read the gas clusters from one Stage A file."""
    root_file = Path(root_file)
    with uproot.open(root_file) as f:
        if "ClusterTree" not in f:
            raise RuntimeError(f"{root_file.name} has no ClusterTree")
        a = f["ClusterTree"].arrays(_BRANCHES, library="np")

    vol = np.array([v.decode() if isinstance(v, bytes) else str(v)
                    for v in a["volume"]])
    keep = np.isin(vol, GAS_VOLUMES)

    out = {k: a[k][keep] for k in _BRANCHES if k != "volume"}
    out["volume"] = vol[keep]

    if max_events is not None:
        sel = out["eventID"] < max_events
        out = {k: v[sel] for k, v in out.items()}
    return out
