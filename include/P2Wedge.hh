#pragma once
// P2Wedge.hh
//
// Parametric P2 wedge profile and gas-envelope dimensions.
// Single C++ source of truth for the P2 geometry constants; the Python
// plotting mirror is scripts/model/p2_model.py and MUST be kept in sync
// (both files carry the same table, see docs/P2_MODEL.md).
//
// Coordinates: gerber frame — the wedge apex (= MESA beam axis) is (0,0),
// all dimensions in mm unless suffixed. Sources:
//   docs/P2_GEOMETRY.md            gerber-derived board / active dimensions
//   design/mechanical/*.stp        gas-frame envelope (bbox analysis)
//   docs/P2_MODEL.md               assumptions table for everything else

#include "G4TwoVector.hh"
#include <vector>

namespace P2 {

// ── Board outline (Edge_Cuts gerber, exact) ──────────────────────────────
constexpr double kBoardRIn        = 95.0;
constexpr double kBoardROut       = 650.0;
constexpr double kBoardEdgeOffset = 10.0;    // radial edges pushed outward
constexpr double kBoardTopCut     = 540.01;  // board exists at y <= this
constexpr double kSectorDeg       = 60.0;

// ── Readout PCB stack (Stack_Up_P2.txt, exact) ───────────────────────────
constexpr double kTCuF = 0.018;   // F.Cu readout pads
constexpr double kTFR4 = 0.200;   // FR4 core
constexpr double kTCuB = 0.018;   // B.Cu ground plane

// ── Gas frame (4 mm drift variant, P2_Frame_V2.0.stp — exact, sectioned in
//    FreeCAD 2026-10-06, docs/P2_GEOMETRY.md §5). One 8 mm body on the
//    board; its opening steps out by 3 mm at a ledge 4.0 mm above the board,
//    and the drift-cathode foil sits on that ledge. Same section at every phi.
//    No profile reaches the board's y = 540 chord, so the top cut is inert.
//    Not modelled: the ~2.4 mm glue groove in the board face, the Ø2 mm gas
//    channel at z = 3.5, the mounting ears outside r = 615. ─────────────────
constexpr double kFrameH          = 8.0;     // board face -> window face
constexpr double kFrameLedgeZ     = 4.0;     // board face -> drift-foil ledge
constexpr double kFrameRIn        = 95.0;    // footprint
constexpr double kFrameROut       = 615.0;
constexpr double kFrameEdgeOffset = 10.0;
constexpr double kFrameTopCut     = kBoardTopCut;
constexpr double kOpenRIn         = 107.0;   // window side (ledge -> window):
constexpr double kOpenROut        = 603.0;   //   window, front gas, drift foil
constexpr double kOpenEdgeOffset  = 3.0;
constexpr double kOpenTopCut      = kBoardTopCut;
constexpr double kDriftOpenRIn    = 110.0;   // board side (board -> ledge):
constexpr double kDriftOpenROut   = 600.0;   //   drift gas, mesh, amp gap
constexpr double kDriftOpenEdgeOffset = 0.0;

// ── Active area (bulk mask, exact; informational) ────────────────────────
constexpr double kActiveRIn  = 119.87;
constexpr double kActiveROut = 589.80;

// ── Pad zone: the copper pad field (P2PadMap.hh) with straight edges 5.0 mm
//    inside the 0 / 60 deg lines (5.00 / 5.12 mm measured from the pad
//    table). The Saclay pillar field fills the same zone. The carbon back
//    frame runs from the board outline in to this zone (Alexandra
//    2026-10-06: "exactly the perimeter of the bare PCB without the active
//    zone"), so the back gas and back window use it as their outline. ───────
constexpr double kZoneRIn         = 114.998;
constexpr double kZoneROut        = 594.872;
constexpr double kZoneEdgeOffset  = -5.0;    // negative = inward

// Closed CCW polygon of a wedge profile: annular sector r in [rIn, rOut],
// phi in [0, 60] deg, both radial edges offset outward perpendicular by
// edgeOffset, truncated by the chord y <= topCut. Passing the board
// parameters reproduces the Edge_Cuts profile (docs/P2_GEOMETRY.md §1).
std::vector<G4TwoVector> WedgeOutline(double rIn, double rOut,
                                      double edgeOffset, double topCut,
                                      int nArc = 96);

// Area centroid of a simple polygon.
G4TwoVector Centroid(const std::vector<G4TwoVector>& poly);

// Polygon scaled by s about point c: p' = c + s*(p - c).
std::vector<G4TwoVector> ScaleAbout(const std::vector<G4TwoVector>& poly,
                                    const G4TwoVector& c, double s);

}  // namespace P2
