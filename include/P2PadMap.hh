#pragma once
// P2PadMap.hh  --  GENERATED, DO NOT EDIT BY HAND
//
// Regenerate with:
//     python scripts/gerber/extract_readout_pattern.py --write
//
// The P2 F.Cu readout pad field, measured feature-by-feature out of
//     design/gerbers/P2_BASKET_Apr26/Gerber/P2_BASKET-F_Cu.gbr
// (KiCad 9.0.2, 2026-01-20).  1280 pad polygons, no exceptions.
//
// The artwork is an exactly regular POLAR grid, which is what makes it
// cheap to build as real geometry:
//
//   * 42 rings on a constant radial pitch of 11.428605 mm
//     (max deviation of any ring from the fitted line: 0.0005 mm)
//   * every pad is 11.30170 mm tall radially (spread 0.00048 mm),
//     leaving a 0.12691 mm radial gap between rings
//   * within a ring the pads are uniformly spaced in phi to
//     5.92e-06 rad -- i.e. exactly, at the 1e-6 mm coordinate
//     resolution of the gerber.  Each ring is therefore ONE G4PVReplica.
//   * each pad is a true annular sector: its polygon area is 0.9997 of
//     r_c * dphi * dr, so a G4Tubs segment is not an approximation.
//
// Pads per ring grows with radius to keep the cell area constant
// ("equalized" pads): 9 on the inner ring, 51 on the outer,
// 1280 channels total = 10 connectors x 128.
//
// Angles are RADIANS, lengths are mm.  phiFirst is the CENTRE of the first
// pad of the ring; the replica envelope therefore starts at
// phiFirst - dPhiPitch/2.

namespace P2 {

struct PadRing {
    int    n;           // pads in this ring
    double rIn;         // pad inner radius  [mm]
    double rOut;        // pad outer radius  [mm]
    double dPhiPitch;   // pad-to-pad angular pitch [rad]
    double phiFirst;    // centre of the first pad  [rad]
    double dPhiPad;     // angular width of the copper pad itself [rad]
};

// Total copper in the pad field: 1697.15 cm^2 over a 1715.89 cm^2
// replica envelope, i.e. 0.9891 of it (the azimuthal gaps).  Folding in
// the radial gaps between rings gives 0.9781 of the full pad-field
// annulus -- that is the number a homogenized F.Cu sheet should carry there.
constexpr int kNPadRings = 42;
constexpr int kNPads     = 1280;

constexpr PadRing kPadRings[kNPadRings] = {
    {  9,  114.997990,  126.300854, 0.106693150, 0.096316382, 0.105679015 },
    { 11,  126.426964,  137.729473, 0.088009777, 0.083087008, 0.087080418 },
    { 12,  137.855115,  149.158001, 0.081221307, 0.076450620, 0.080362959 },
    { 13,  149.283782,  160.586569, 0.075401012, 0.070792953, 0.074605538 },
    { 14,  160.712269,  172.015111, 0.070355221, 0.065914656, 0.069612112 },
    { 15,  172.141438,  183.443779, 0.065940101, 0.061665579, 0.065244105 },
    { 16,  183.569826,  194.872244, 0.062045121, 0.057930282, 0.061389132 },
    { 17,  194.998638,  206.300854, 0.058583366, 0.054620145, 0.057964002 },
    { 18,  206.427195,  217.729415, 0.055486177, 0.051670755, 0.054899166 },
    { 19,  217.855952,  229.158112, 0.052700077, 0.049022137, 0.052142228 },
    { 20,  229.284586,  240.586624, 0.050178942, 0.046629923, 0.049648046 },
    { 21,  240.712992,  252.015178, 0.047888346, 0.044460957, 0.047381693 },
    { 22,  252.142008,  263.443787, 0.045797124, 0.042486203, 0.045311371 },
    { 23,  263.570617,  274.872348, 0.043881006, 0.040677872, 0.043415872 },
    { 24,  274.999088,  286.300941, 0.042117897, 0.039017546, 0.041671718 },
    { 25,  286.427648,  297.729537, 0.040491272, 0.037486537, 0.040061912 },
    { 26,  297.856206,  309.158055, 0.038985367, 0.036072159, 0.038572077 },
    { 27,  309.284869,  320.586651, 0.037587680, 0.034760262, 0.037188854 },
    { 28,  320.713397,  332.015211, 0.036286493, 0.033540254, 0.035901230 },
    { 29,  332.142151,  343.443852, 0.035072116, 0.032403065, 0.034699967 },
    { 30,  343.570682,  354.872399, 0.033936337, 0.031340548, 0.033576369 },
    { 31,  354.999257,  366.300901, 0.032871864, 0.030345746, 0.032522895 },
    { 32,  366.427926,  377.729497, 0.031872159, 0.029411201, 0.031533745 },
    { 33,  377.856493,  389.158131, 0.030931424, 0.028532889, 0.030602951 },
    { 34,  389.285052,  400.586627, 0.030044474, 0.027706011, 0.029725390 },
    { 35,  400.713587,  412.015176, 0.029206863, 0.026926197, 0.028896652 },
    { 36,  412.142236,  423.443810, 0.028414838, 0.026188492, 0.028112941 },
    { 37,  423.570782,  434.872393, 0.027664536, 0.025489650, 0.027370635 },
    { 38,  434.999347,  446.300954, 0.026952858, 0.024827561, 0.026666436 },
    { 39,  446.427956,  457.729540, 0.026276876, 0.024197726, 0.025997405 },
    { 40,  457.856568,  469.158072, 0.025633904, 0.023600612, 0.025361263 },
    { 41,  469.285183,  480.586668, 0.025021645, 0.023031912, 0.024755792 },
    { 42,  480.713736,  492.015243, 0.024437936, 0.022489629, 0.024178095 },
    { 43,  492.142272,  503.443799, 0.023880892, 0.021972808, 0.023627049 },
    { 44,  503.570930,  514.872350, 0.023348668, 0.021479285, 0.023100069 },
    { 45,  514.999450,  526.300939, 0.022839547, 0.021006756, 0.022596683 },
    { 46,  526.428024,  537.729513, 0.022352264, 0.020554874, 0.022114540 },
    { 47,  537.856598,  549.158107, 0.021885198, 0.020122418, 0.021652419 },
    { 48,  549.285165,  560.586651, 0.021437308, 0.019707509, 0.021209351 },
    { 49,  560.713739,  572.015196, 0.021007399, 0.019309211, 0.020783789 },
    { 50,  572.142331,  583.443788, 0.020594384, 0.018926615, 0.020375150 },
    { 51,  583.570927,  594.872373, 0.020197351, 0.018558757, 0.019982180 },
};

// ── Homogenized copper, by radial band ──────────────────────────────────────
//
// Everything that is NOT the pad field is irregular artwork and stays
// homogenized -- but per band, not as one board-wide average.  Coverage is
// the copper area fraction of the band, measured on the exact vector
// geometry (see the script; rasterized at pixel centres).
//
// Two facts these numbers encode that a single average destroys:
//   * F.Cu carries NO copper at all inside r = 115.0 mm.  The board
//     starts at r = 95 mm, so a board-wide sheet puts ~98 % copper across a
//     20 mm annulus that is bare laminate.
//   * B.Cu is not a ground plane -- it is 30805 stroked 0.125 mm
//     signal traces, and their density climbs steadily with radius as more
//     channels are gathered toward the connectors.
//
// fCu is used only when the pattern is switched off (--homogenized-readout);
// with the pattern on, the bands inside the pad field are replaced by the
// real pads and only the first and last bands are built from this table.

struct CuBand {
    double rIn;    // [mm]
    double rOut;   // [mm]
    double fCu;    // F.Cu area coverage in this band
    double bCu;    // B.Cu area coverage in this band
};

constexpr int kNCuBands = 10;
constexpr CuBand kCuBands[kNCuBands] = {
    {  95.000000, 114.997990, 0.123201, 0.123179 },
    { 114.997990, 174.982288, 0.829969, 0.048873 },
    { 174.982288, 234.966586, 0.866003, 0.073919 },
    { 234.966586, 294.950884, 0.891748, 0.105350 },
    { 294.950884, 354.935182, 0.905973, 0.133271 },
    { 354.935182, 414.919479, 0.916853, 0.161352 },
    { 414.919479, 474.903777, 0.925714, 0.189415 },
    { 474.903777, 534.888075, 0.930844, 0.216760 },
    { 534.888075, 594.872373, 0.937579, 0.246041 },
    { 594.872373, 650.000000, 0.179940, 0.234505 },
};

// Pad-field boundaries, so the geometry code does not have to re-derive them.
constexpr double kPadFieldRIn  = 114.997990;
constexpr double kPadFieldROut = 594.872373;

}  // namespace P2
