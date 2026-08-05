// P2Wedge.cc
// See P2Wedge.hh. The outline traversal is a direct port of
// scripts/gerber/p2_wedge_model.py::outline(), which was validated against
// the Edge_Cuts gerber (docs/figures/p2_wedge_outline.png).

#include "P2Wedge.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace P2 {

std::vector<G4TwoVector> WedgeOutline(double rIn, double rOut,
                                      double edgeOffset, double topCut,
                                      int nArc) {
    const double d2r  = M_PI / 180.0;
    const double phi1 = kSectorDeg * d2r;

    // Angles where the offset edges / top cut meet the arcs.
    //
    // The two radial edges are the rays at 0 deg and 60 deg pushed outward
    // perpendicular by edgeOffset, i.e. the lines p.n = edgeOffset with n the
    // outward normal.  A circle of radius r meets the 60 deg edge at
    // 60 deg + asin(off/r), and the bottom edge at -asin(off/r).
    const double aInLo   = -std::asin(edgeOffset / rIn);
    const double aInHi   =  phi1 - aInLo;                  // 60deg + asin(off/rIn)
    const double aOutLo  = -std::asin(edgeOffset / rOut);
    const double aOutEdge = phi1 + std::asin(edgeOffset / rOut);
    if (topCut >= rOut)
        throw std::runtime_error("P2::WedgeOutline: topCut >= rOut not supported");
    const double aOutCut = std::asin(topCut / rOut);       // chord meets outer arc

    // Does the y = topCut chord actually clip the outer arc, or does the
    // 60 deg edge get there first?  For the board profile the chord clips
    // (aOutCut = 56.2deg < aOutEdge = 60.9deg); for the gas frame and the
    // frame opening it does NOT (the whole profile sits below topCut), and
    // forcing a chord corner in anyway used to send the outer arc past the
    // 60 deg edge and fold the polygon over itself — G4ExtrudedSolid then
    // died in triangularisation.  Fixed 2026-08-05.
    const bool clipped  = (aOutCut < aOutEdge);
    const double aOutHi = clipped ? aOutCut : aOutEdge;

    std::vector<G4TwoVector> pts;
    pts.reserve(2 * nArc + 6);

    auto arc = [&](double r, double a0, double a1) {
        for (int i = 0; i <= nArc; ++i) {
            double a = a0 + (a1 - a0) * i / nArc;
            pts.emplace_back(r * std::cos(a), r * std::sin(a));
        }
    };

    // bottom edge (y = -edgeOffset), inner corner -> outer corner
    pts.emplace_back(rIn  * std::cos(aInLo),  -edgeOffset);
    pts.emplace_back(rOut * std::cos(aOutLo), -edgeOffset);
    // outer arc CCW to the top chord, or straight to the 60deg edge
    arc(rOut, aOutLo, aOutHi);
    // top chord across to the 60deg offset edge (only when it clips)
    if (clipped) {
        const double xTop = (topCut * std::cos(phi1) - edgeOffset) / std::sin(phi1);
        pts.emplace_back(xTop, topCut);
    }
    // down the 60deg edge to the inner arc, then inner arc CW back to start
    pts.emplace_back(rIn * std::cos(aInHi), rIn * std::sin(aInHi));
    arc(rIn, aInHi, aInLo);

    // G4ExtrudedSolid rejects duplicated vertices; they appear where an
    // arc endpoint coincides with an edge corner (e.g. edgeOffset = 0).
    // Drop consecutive duplicates, including across the wrap-around.
    {
        const double tol2 = 1e-12;  // (1e-6 mm)^2
        std::vector<G4TwoVector> ded;
        ded.reserve(pts.size());
        for (const auto& p : pts) {
            if (!ded.empty() && (p - ded.back()).mag2() < tol2) continue;
            ded.push_back(p);
        }
        if (ded.size() > 1 && (ded.front() - ded.back()).mag2() < tol2)
            ded.pop_back();
        pts.swap(ded);
    }

    // Guard: G4ExtrudedSolid wants a counter-clockwise, non-degenerate polygon.
    double area2 = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        const auto& q = pts[(i + 1) % pts.size()];
        area2 += p.x() * q.y() - q.x() * p.y();
    }
    if (area2 < 0.0) std::reverse(pts.begin(), pts.end());

    // Guard: a self-intersecting profile makes G4ExtrudedSolid abort inside
    // AddGeneralPolygonFacets ("Triangularisation has failed") with a core
    // dump and no indication of which solid was at fault. O(n^2) on ~200
    // points, once per profile at construction — free, and it turns a crash
    // into a message that names the parameters.
    {
        const size_t n = pts.size();
        auto cross = [](const G4TwoVector& o, const G4TwoVector& p,
                        const G4TwoVector& q) {
            return (p.x()-o.x())*(q.y()-o.y()) - (p.y()-o.y())*(q.x()-o.x());
        };
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 2; j < n; ++j) {
                if (i == 0 && j == n - 1) continue;      // shared vertex
                const auto &a = pts[i], &b = pts[(i+1)%n];
                const auto &c = pts[j], &d = pts[(j+1)%n];
                const bool hit = ((cross(c,d,a) > 0) != (cross(c,d,b) > 0)) &&
                                 ((cross(a,b,c) > 0) != (cross(a,b,d) > 0));
                if (hit) {
                    throw std::runtime_error(
                        "P2::WedgeOutline: self-intersecting profile (edges " +
                        std::to_string(i) + "/" + std::to_string(j) +
                        ") for rIn=" + std::to_string(rIn) +
                        " rOut=" + std::to_string(rOut) +
                        " edgeOffset=" + std::to_string(edgeOffset) +
                        " topCut=" + std::to_string(topCut));
                }
            }
        }
    }

    return pts;
}

G4TwoVector Centroid(const std::vector<G4TwoVector>& poly) {
    double a2 = 0.0, cx = 0.0, cy = 0.0;
    for (size_t i = 0; i < poly.size(); ++i) {
        const auto& p = poly[i];
        const auto& q = poly[(i + 1) % poly.size()];
        const double w = p.x() * q.y() - q.x() * p.y();
        a2 += w;
        cx += (p.x() + q.x()) * w;
        cy += (p.y() + q.y()) * w;
    }
    return { cx / (3.0 * a2), cy / (3.0 * a2) };
}

std::vector<G4TwoVector> ScaleAbout(const std::vector<G4TwoVector>& poly,
                                    const G4TwoVector& c, double s) {
    std::vector<G4TwoVector> out;
    out.reserve(poly.size());
    for (const auto& p : poly)
        out.emplace_back(c.x() + s * (p.x() - c.x()),
                         c.y() + s * (p.y() - c.y()));
    return out;
}

}  // namespace P2
