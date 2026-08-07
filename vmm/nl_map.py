"""Pad adjacency for the P2 wedge: *readout-channel* vs *geometric* neighbors.

The VMM3a neighbor logic (NL) forces readout of the channels adjacent **in
chip-channel number** to a channel that fired. On a strip detector those are
also the physical neighbors; on the P2 pad wedge they are not, because the
pad → connector-pin assignment snakes through the pad plane. This module
builds both neighbor relations from `design/mapping/` and quantifies the
mismatch, so the emulator can model NL as the chip really does it instead of
as "geometric neighbors".

Two mapping revisions ship in the repo:

  A  design/mapping/connector_<i>.txt          7 columns (Radius, Phi, AbsIndex)
  B  design/mapping/Mapping/connector_<i>.txt  9 columns (adds X, Y, PadName)

They describe the **same 1280 pads at the same positions** (verified here:
max nearest-pad mismatch 2.5 um, i.e. print precision) but assign *channels*
to those pads in almost completely different orders (11/1280 (connector,
channel) pairs agree). Everything channel-level — neighbor logic, dead-channel
masks, per-channel comparisons with test-beam data — depends on which one the
DAQ uses. See docs/SIM_CAMPAIGN_PLAN.md §10.6.

Usage:
    python vmm/nl_map.py                # print the adjacency report
    from nl_map import PadMap
    pm = PadMap.load("B")
    pm.channel_neighbors(pad_id)        # VMM NL set (chan +-1)
    pm.geometric_neighbors(pad_id)      # physically touching pads
"""

from __future__ import annotations

import glob
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Sequence

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
MAPDIR = os.path.join(REPO, "design", "mapping")

CHANNELS_PER_CHIP = 64      # VMM3a
CHANNELS_PER_CONNECTOR = 128


@dataclass
class PadMap:
    """Pad geometry + channel assignment for one mapping revision."""

    revision: str
    conn: np.ndarray          # connector index 0..9
    chan: np.ndarray          # channel index within connector, 0..127
    r: np.ndarray             # mm
    phi: np.ndarray           # rad
    ring: np.ndarray          # radial ring index, 0 = innermost

    # ------------------------------------------------------------------
    @classmethod
    def load(cls, revision: str = "B") -> "PadMap":
        """revision 'A' = 7-column files, 'B' = 9-column Mapping/ files."""
        sub = "Mapping" if revision.upper() == "B" else ""
        rows = []
        for path in sorted(glob.glob(os.path.join(MAPDIR, sub, "connector_*.txt"))):
            for line in open(path).read().splitlines()[1:]:
                p = line.split()
                if revision.upper() == "B" and len(p) >= 9:
                    rows.append((int(p[0]), int(p[2]), float(p[7]), float(p[8])))
                elif revision.upper() == "A" and len(p) >= 6:
                    rows.append((int(p[0]), int(p[1]), float(p[4]), float(p[5])))
        if not rows:
            raise FileNotFoundError(f"no mapping files for revision {revision}")
        a = np.array(rows, dtype=float)
        r = a[:, 2]
        rings = np.unique(np.round(r, 3))
        ridx = {round(v, 3): i for i, v in enumerate(rings)}
        return cls(revision.upper(), a[:, 0].astype(int), a[:, 1].astype(int),
                   r, a[:, 3], np.array([ridx[round(v, 3)] for v in r]))

    def __len__(self) -> int:
        return len(self.r)

    @property
    def x(self) -> np.ndarray:
        return self.r * np.cos(self.phi)

    @property
    def y(self) -> np.ndarray:
        return self.r * np.sin(self.phi)

    # ------------------------------------------------------------------
    def _index(self) -> Dict[tuple, int]:
        if not hasattr(self, "_idx_cache"):
            self._idx_cache = {(int(c), int(ch)): i
                               for i, (c, ch) in enumerate(zip(self.conn, self.chan))}
        return self._idx_cache

    def channel_neighbors(self, pad: int, cross_chip: bool = True) -> List[int]:
        """Pads read out by VMM neighbor logic when `pad` fires.

        NL acts on chan +-1 within a chip. `cross_chip=False` breaks the
        relation at the 64-channel chip boundary (whether the P2 hybrid wires
        the inter-chip NL lines is NEEDS-DATA).
        """
        idx = self._index()
        c, ch = int(self.conn[pad]), int(self.chan[pad])
        out = []
        for s in (-1, 1):
            n = ch + s
            if not (0 <= n < CHANNELS_PER_CONNECTOR):
                continue
            if not cross_chip and (n // CHANNELS_PER_CHIP) != (ch // CHANNELS_PER_CHIP):
                continue
            j = idx.get((c, n))
            if j is not None:
                out.append(j)
        return out

    def geometric_neighbors(self, pad: int) -> List[int]:
        """Physically touching pads: same ring +-1 in phi, or adjacent ring
        with angular overlap. Corner-only contacts are excluded."""
        if not hasattr(self, "_geo_cache"):
            self._geo_cache = self._build_geometric()
        return self._geo_cache[pad]

    def _half_width(self) -> Dict[int, float]:
        hw = {}
        for i in np.unique(self.ring):
            phis = np.sort(self.phi[self.ring == i])
            hw[int(i)] = float(np.median(np.diff(phis)) / 2) if len(phis) > 1 else 0.05
        return hw

    def _build_geometric(self) -> List[List[int]]:
        hw = self._half_width()
        nmax = int(self.ring.max())
        out: List[List[int]] = []
        for i in range(len(self)):
            nb: set = set()
            for dr in (-1, 0, 1):
                jr = int(self.ring[i]) + dr
                if jr < 0 or jr > nmax:
                    continue
                sel = np.where(self.ring == jr)[0]
                dphi = np.abs(self.phi[sel] - self.phi[i])
                if dr == 0:                       # azimuthal neighbors
                    cand = sel[(dphi > 1e-9) & (dphi < 3.0 * hw[jr])]
                else:                             # radial: require angular overlap
                    cand = sel[dphi < hw[jr] + hw[int(self.ring[i])]]
                nb.update(cand.tolist())
            nb.discard(i)
            out.append(sorted(nb))
        return out

    def distance(self, i: int, j: int) -> float:
        return math.hypot(self.x[i] - self.x[j], self.y[i] - self.y[j])


# ----------------------------------------------------------------------
def report(pm: PadMap, cross_chip: bool = True) -> dict:
    """Quantify NL-vs-geometry mismatch. Returns the numbers quoted in docs."""
    n = len(pm)
    geo = [set(pm.geometric_neighbors(i)) for i in range(n)]
    d_all, is_phys, is_radial = [], [], []
    covered, total, wasted = [], [], []
    for i in range(n):
        nl = pm.channel_neighbors(i, cross_chip=cross_chip)
        for j in nl:
            d_all.append(pm.distance(i, j))
            is_phys.append(j in geo[i])
            is_radial.append(abs(int(pm.ring[i]) - int(pm.ring[j])) == 1)
        covered.append(len(set(nl) & geo[i]))
        total.append(len(geo[i]))
        wasted.append(len(set(nl) - geo[i]))
    d_all = np.array(d_all); is_phys = np.array(is_phys)
    covered = np.array(covered); total = np.array(total); wasted = np.array(wasted)

    res = dict(
        revision=pm.revision, n_pads=n,
        n_nl_pairs=int(len(d_all)),
        frac_nl_physical=float(is_phys.mean()),
        median_nl_distance_mm=float(np.median(d_all)),
        p90_nl_distance_mm=float(np.percentile(d_all, 90)),
        max_nl_distance_mm=float(d_all.max()),
        median_nonadjacent_distance_mm=float(np.median(d_all[~is_phys])) if (~is_phys).any() else 0.0,
        frac_nl_radial_step=float(np.mean(is_radial)),
        mean_geometric_neighbors=float(total.mean()),
        mean_neighbors_covered=float(covered.mean()),
        neighborhood_coverage=float(covered.sum() / total.sum()),
        frac_pads_no_real_neighbor_in_nl=float((covered == 0).mean()),
        wasted_channels_per_pad=float(wasted.mean()),
    )
    return res


def sharing_catch_probability(pm: PadMap, p_radial: float = 0.49,
                              cross_chip: bool = True) -> dict:
    """P(the charge-sharing partner pad is in the NL set).

    The pad cell is 11.43 mm (radial) x 11.86 mm (azimuthal), so a cloud that
    spills over an edge does so across a radial edge with probability
    ~11.86/(11.43+11.86) = 0.51 ... near enough to a coin flip; `p_radial` is
    the split (sharing probability scales inversely with the cell dimension in
    that direction). Each of the four edge directions gets weight
    p_radial/2 or (1-p_radial)/2, and the direction counts as "caught" in
    proportion to how many of the neighbors on that edge the NL set contains.
    """
    n = len(pm)
    w = {"in": p_radial / 2, "out": p_radial / 2,
         "phi-": (1 - p_radial) / 2, "phi+": (1 - p_radial) / 2}
    caught = np.zeros(n)
    for i in range(n):
        nl = set(pm.channel_neighbors(i, cross_chip=cross_chip))
        by_dir: Dict[str, List[int]] = {k: [] for k in w}
        for j in pm.geometric_neighbors(i):
            dring = int(pm.ring[j]) - int(pm.ring[i])
            if dring == 0:
                by_dir["phi+" if pm.phi[j] > pm.phi[i] else "phi-"].append(j)
            else:
                by_dir["out" if dring > 0 else "in"].append(j)
        tot = 0.0
        for d, ws in w.items():
            cand = by_dir[d]
            if not cand:          # detector edge: no partner to catch
                continue
            tot += ws * len(set(cand) & nl) / len(cand)
        caught[i] = tot
    return {"p_catch_mean": float(caught.mean()),
            "p_catch_min": float(caught.min()),
            "p_catch_max": float(caught.max()),
            "frac_pads_zero_catch": float((caught == 0).mean())}


def _print(res: dict) -> None:
    print(f"\n--- revision {res['revision']}  ({res['n_pads']} pads, "
          f"{res['n_nl_pairs']} channel+-1 pairs) ---")
    print(f"  NL partner physically adjacent : {100*res['frac_nl_physical']:.1f} %")
    print(f"  NL partner distance            : median {res['median_nl_distance_mm']:.1f} mm, "
          f"90th pct {res['p90_nl_distance_mm']:.1f} mm, max {res['max_nl_distance_mm']:.0f} mm")
    print(f"    (of the non-adjacent ones    : median {res['median_nonadjacent_distance_mm']:.0f} mm)")
    print(f"  step direction                 : {100*res['frac_nl_radial_step']:.0f} % radial, "
          f"{100*(1-res['frac_nl_radial_step']):.0f} % azimuthal")
    print(f"  physical neighbors per pad     : {res['mean_geometric_neighbors']:.2f}")
    print(f"  of those reached by NL         : {res['mean_neighbors_covered']:.2f} "
          f"=> {100*res['neighborhood_coverage']:.0f} % neighborhood coverage")
    print(f"  pads whose NL set contains no  \n"
          f"    true neighbor at all         : {100*res['frac_pads_no_real_neighbor_in_nl']:.1f} %")
    print(f"  non-neighbor channels forced   : {res['wasted_channels_per_pad']:.2f} per fired pad")


def compare_revisions() -> None:
    """Do A and B describe the same pads? the same channel assignment?"""
    A, B = PadMap.load("A"), PadMap.load("B")
    PA = np.c_[A.x, A.y]; PB = np.c_[B.x, B.y]
    d = np.sqrt(((PA[:, None, :] - PB[None, :, :]) ** 2).sum(-1)).min(axis=1)
    print(f"\npad POSITIONS: max nearest-pad mismatch A->B = {d.max()*1000:.1f} um "
          f"({(d > 0.01).sum()} pads off by >10 um) => same pad plane")
    idxA = {(int(c), int(ch)): i for i, (c, ch) in enumerate(zip(A.conn, A.chan))}
    same = 0
    for j, (c, ch) in enumerate(zip(B.conn, B.chan)):
        i = idxA.get((int(c), int(ch)))
        if i is not None and abs(A.r[i] - B.r[j]) < 0.01 and abs(A.phi[i] - B.phi[j]) < 1e-4:
            same += 1
    print(f"pad CHANNEL ASSIGNMENT: (connector,channel) -> same pad for "
          f"{same}/{len(B)} channels ({100*same/len(B):.1f} %) => the revisions "
          f"disagree on the readout order, not the geometry")


if __name__ == "__main__":
    print(__doc__.split("Usage:")[0].strip()[:0] or "", end="")
    print("P2 pad map: VMM neighbor-logic channels vs physical pad neighbors")
    compare_revisions()
    for rev in ("A", "B"):
        pm = PadMap.load(rev)
        _print(report(pm))
        s = sharing_catch_probability(pm)
        print(f"  P(charge-sharing partner is in the NL set)"
              f": {100*s['p_catch_mean']:.0f} % "
              f"(range {100*s['p_catch_min']:.0f}-{100*s['p_catch_max']:.0f} %)")
    pm = PadMap.load("B")
    r_nochip = report(pm, cross_chip=False)
    print(f"\nrevision B with NL broken at the 64-channel chip boundary: "
          f"{r_nochip['n_nl_pairs']} pairs "
          f"({100*r_nochip['frac_nl_physical']:.1f} % physical, coverage "
          f"{100*r_nochip['neighborhood_coverage']:.0f} %)")
