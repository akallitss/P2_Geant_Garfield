# Local patches to the pinned Garfield++ (927e5c21)

Applied by `scripts/setup_garfield.sh --build` after the checkout.

- `interfacepanels_periodic_leaf.patch` (2026-10-14). With
  `ComponentFieldMap::EnableInterfaceCrossings()` on a periodic map,
  `InterfacePanels::CrossedFolded` bisects a step that crosses the periodic
  boundary down to its maximum depth and then treats the last leaf as affine:
  its folded end points sit on opposite faces of the cell, so the chord spans
  the whole cell and reports a crossing with any panel on the way. In the P2
  mesh cell this stopped about half of the electrons on the cell faces
  (x or y = ±68 µm; transparency 45 % instead of ~90 %). The patch skips a
  leaf whose mapped length exceeds twice its real length (it straddles the
  boundary and is ~1e-6 of the step). Worth reporting upstream.
