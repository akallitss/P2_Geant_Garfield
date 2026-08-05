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

// ── Gas frame (from P2_Frame_V1_3mm.stp: body z = 0..8 mm on the board,
//    drift-mylar ledge at z = 3 mm, inner walls ~2.5 mm, edge walls 7.5 mm,
//    footprint r ~ 104.5..605.5). Opening dims are read off dominant point
//    clusters and are approximate — see docs/P2_MODEL.md. ─────────────────
constexpr double kFrameRIn        = 104.5;
constexpr double kFrameROut       = 605.5;
constexpr double kFrameEdgeOffset = 7.5;
constexpr double kFrameTopCut     = 537.5;
constexpr double kOpenRIn         = 107.0;   // frame inner opening = window
constexpr double kOpenROut        = 603.0;
constexpr double kOpenEdgeOffset  = 0.0;
constexpr double kOpenTopCut      = 535.0;

// ── Active area (bulk mask, exact; informational) ────────────────────────
constexpr double kActiveRIn  = 119.87;
constexpr double kActiveROut = 589.80;

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
