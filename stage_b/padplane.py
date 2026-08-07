"""padplane.py — where does an electron land?

`vmm/nl_map.py` already loads the pad geometry from `design/mapping/` and
gives pad centres plus the channel assignment. What it does not give is the
inverse: point (x, y) -> pad. That is what Stage B needs, so it is built here
from the same source rather than from a second copy of the geometry.

The P2 pad plane is polar: 42 concentric rings on an 11.4290 mm pitch, each
divided into annular-sector pads, 1280 in total. So the lookup is two 1-D
searches — ring from radius, sector from azimuth — not a 2-D one.

A NOTE ON SCALE, because it decides how much of this matters. Pads are about
11.4 mm radially and ~11-12 mm azimuthally at mid-radius. Transverse
diffusion over the full 3 mm gap is ~140 um. So the charge cloud is roughly
1 % of a pad: sharing happens only when a track passes within a few hundred
microns of an edge, which is exactly why the measured figure is ~1.1 pads per
hit and why neighbour logic (designed for multi-strip clusters) does not
transfer. Do not "improve" this into a 2-D charge-integration model without
first checking it changes anything.

MAPPING REVISION. Two revisions ship in the repo and assign *channels* to the
same pads in almost completely different orders (11/1280 agree). Pad
positions are identical to 2.5 um, so the geometry here is revision-
independent, but anything channel-level is not. The revision is recorded in
the output; see SIM_CAMPAIGN_PLAN §10.6 — which revision the DAQ uses is
still an open question.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np

# Reuse the loader rather than re-reading design/mapping.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
if os.path.join(_REPO, "vmm") not in sys.path:
    sys.path.insert(0, os.path.join(_REPO, "vmm"))

from nl_map import PadMap  # noqa: E402


@dataclass
class PadPlane:
    """Point -> pad lookup over the polar pad map."""

    revision: str
    pad_r: np.ndarray          # pad centre radius [mm]
    pad_phi: np.ndarray        # pad centre azimuth [rad]
    ring_of_pad: np.ndarray    # ring index per pad
    conn: np.ndarray
    chan: np.ndarray

    ring_centres: np.ndarray   # [n_ring]
    ring_edges: np.ndarray     # [n_ring+1]
    _phi_by_ring: list         # per ring: sorted pad phi centres
    _phi_edges: list           # per ring: sector edges
    _pad_by_ring: list         # per ring: pad indices, ordered by phi

    # ------------------------------------------------------------------
    @classmethod
    def load(cls, revision: str = "B") -> "PadPlane":
        pm = PadMap.load(revision)
        r, phi = np.asarray(pm.r, float), np.asarray(pm.phi, float)

        centres = np.unique(np.round(r, 4))
        # Edges midway between ring centres; the outermost half-pitches are
        # extrapolated from the neighbouring gap.
        mid = 0.5 * (centres[1:] + centres[:-1])
        edges = np.concatenate(([centres[0] - (mid[0] - centres[0])],
                                mid,
                                [centres[-1] + (centres[-1] - mid[-1])]))

        ring_of_pad = np.searchsorted(edges, r, side="right") - 1
        ring_of_pad = np.clip(ring_of_pad, 0, len(centres) - 1)

        phi_by_ring, phi_edges, pad_by_ring = [], [], []
        for k in range(len(centres)):
            idx = np.flatnonzero(ring_of_pad == k)
            order = np.argsort(phi[idx])
            idx = idx[order]
            p = phi[idx]
            phi_by_ring.append(p)
            pad_by_ring.append(idx)
            if len(p) > 1:
                m = 0.5 * (p[1:] + p[:-1])
                e = np.concatenate(([p[0] - (m[0] - p[0])], m,
                                    [p[-1] + (p[-1] - m[-1])]))
            else:
                e = np.array([p[0] - 1e-3, p[0] + 1e-3]) if len(p) else np.array([0.0, 0.0])
            phi_edges.append(e)

        return cls(revision=revision, pad_r=r, pad_phi=phi,
                   ring_of_pad=ring_of_pad,
                   conn=np.asarray(pm.conn), chan=np.asarray(pm.chan),
                   ring_centres=centres, ring_edges=edges,
                   _phi_by_ring=phi_by_ring, _phi_edges=phi_edges,
                   _pad_by_ring=pad_by_ring)

    # ------------------------------------------------------------------
    @property
    def n_pads(self) -> int:
        return len(self.pad_r)

    @property
    def ring_pitch_mm(self) -> float:
        return float(np.median(np.diff(self.ring_centres)))

    def locate(self, x, y) -> np.ndarray:
        """Pad index per point, or -1 for points off the pad plane.

        -1 is returned rather than a nearest-pad guess: electrons landing
        outside the instrumented area are a real category (they hit the
        fan-out or the frame) and silently snapping them onto the edge pad
        would inflate the edge rings.
        """
        x = np.atleast_1d(np.asarray(x, float))
        y = np.atleast_1d(np.asarray(y, float))
        r = np.hypot(x, y)
        phi = np.arctan2(y, x)

        out = np.full(r.shape, -1, dtype=np.int64)
        ring = np.searchsorted(self.ring_edges, r, side="right") - 1
        inside = (ring >= 0) & (ring < len(self.ring_centres))

        for k in np.unique(ring[inside]):
            sel = inside & (ring == k)
            if not np.any(sel):
                continue
            e = self._phi_edges[k]
            pads = self._pad_by_ring[k]
            j = np.searchsorted(e, phi[sel], side="right") - 1
            ok = (j >= 0) & (j < len(pads))
            idx = np.flatnonzero(sel)
            out[idx[ok]] = pads[j[ok]]
        return out
