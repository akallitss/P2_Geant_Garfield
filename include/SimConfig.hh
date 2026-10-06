#pragma once
// SimConfig.hh
// Shared simulation configuration passed through the G4 user classes

#include <string>
#include "G4SystemOfUnits.hh"

enum class SimMode {
    kP2Wedge,          // P2 wedge Micromegas, full gas envelope (default)
    kVacuum,
    kFullExperiment,
    kSr90Calibration,
    kSr90NoMM,
    kLSCalib,          // source → air → 1 LS layer
    kBackScintCalib,   // source → air → 1 back plastic scint bar
};

struct SimConfig {
    std::string gas;            // Gas mixture name
    std::string particle;       // Particle type
    double      energy;         // Beam energy [Geant4 internal units = MeV]
    int         nEvents;
    std::string outFile;
    long        seed;
    int         nThreads;
    bool        verbose;
    double      alThickness_mm;
    SimMode     mode = SimMode::kP2Wedge;

    // Photon runs are >=99.9% empty; writing only events with at least one
    // gas cluster keeps the output small. The total thrown count is still
    // reported so normalization stays correct.
    bool        skipEmpty = false;

    // CF4's W-value is quoted anywhere between 35 and 52 eV, which is a
    // genuine systematic rather than a measurement to look up. Overriding the
    // component value here re-derives the mixture W, so the two ends can be
    // bracketed as run points (SIM_CAMPAIGN_PLAN P0.2). <= 0 keeps the table
    // default. Only affects mixtures containing CF4.
    double      w_cf4_eV = -1.0;

    // ── P2 wedge mode ─────────────────────────────────────────────────────
    // Values marked GUESS are not backed by design data — see docs/P2_MODEL.md
    // for the full assumptions table and what to confirm with the collaboration.
    double p2_drift_mm       = 3.964;  // drift GAS, mesh top -> drift foil. The
                                       // frame (ledge 4.0 mm, P2_Frame_V2.0.stp) sits
                                       // on the 150 um Dynamask, level with the mesh,
                                       // so gas = 4.0 - mesh slab (0.036). The nominal
                                       // "4 mm drift" (Alexandra 2026-10-06).
                                       // Campaign will scan 1..4 mm via --drift-gap.
    double p2_amp_um         = 150.0;  // amplification gap — confirmed 150 um (2026-08-05)
    double p2_mesh_wire_um   = 18.0;   // woven SS mesh 45/18: wire diameter and
    double p2_mesh_open_um   = 45.0;   //   opening, pitch 63 um (Alexandra 2026-10-06;
                                       //   supersedes "48x19"). Modelled as a 2*wire slab
                                       //   of effective-density steel
    double p2_front_gap_mm   = 3.879;  // outer window -> drift foil = frame top (8.0) -
                                       // ledge (4.0) - foil (0.121), P2_Frame_V2.0.stp
    double p2_cath_gap_mm    = 0.0;    // 0 = ONE drift-cathode foil (2026-10-06).
                                       // > 0 = the old two-foil model, foils this far apart
    double p2_back_gap_mm    = 1.0;    // PCB back -> back window = carbon back-frame depth.
                                       // Normally 1 mm, up to 3 mm in this production
                                       // (Alexandra 2026-08-05); mylar glued on its rear face.
    double p2_bulge_front_mm = 10.0;   // overpressure sag — Hencky estimate at ~3 mbar,
    double p2_bulge_back_mm  = 10.0;   //   p in 1..10 mbar -> 7..15 mm; 10 mm both
                                       //   sides (Alexandra 2026-10-06)
    double p2_window_um      = 10.0;   // outer containment mylar (Alexandra 2026-10-06)
    std::string p2_pillars   = "saclay"; // mesh-support pillar mask: "saclay" (V1,
                                       //   det1-det4), "cern" (det5), "none"
    double p2_cath_mylar_um  = 120.0;  // drift-cathode mylar (Alexandra 2026-10-06)
    double p2_cath_al_um     = 1.0;    // its aluminization, facing the drift gas
    double p2_fcu_coverage   = 0.983;  // Cu area fraction over the active area, measured
    double p2_bcu_coverage   = 0.174;  //   from the gerbers (B.Cu = signal lines, not a
                                       //   plane; radial 0.03->0.26). Layers modelled as
                                       //   full-thickness slabs of density-scaled copper.
                                       //   scripts/gerber/analyze_cu_coverage.py
    // Gun aim point, mid-active-area, apex/beam axis is x=y=0.
    //
    // This took two attempts, and the reason is worth keeping: the pad plane
    // is POLAR, so there are two independent ways to sit on a pad boundary
    // and fixing one does not fix the other.
    //
    //   * The original (307.4, 177.5) = r 354.97 mm was a round "r = 355,
    //     phi = 30" and landed 35 um from a RADIAL boundary — it is not
    //     inside any pad ring at all.
    //   * Moving it to a ring centre, (302.49, 174.64), fixed the radius and
    //     landed 61 um from the centre of the 126 um AZIMUTHAL gap between
    //     pads 14 and 15 of ring 20. Measured against the real artwork:
    //     <edepPadCu> 4740 eV with only 25.7 % of events depositing in pad
    //     copper at all, against 20383 eV and 100 % at a true pad centre.
    //
    // A pencil beam parked on a boundary makes "which pad has the most
    // charge" a coin flip and inflates pad multiplicity, invisibly -- MX17
    // lost a first result to precisely this
    // (MX17_Geant/design/RESPONSE_SIM_PLAN.md §7).
    //
    // The value below is a pad centre in BOTH coordinates, taken from the
    // gerber-derived table in include/P2PadMap.hh. That table is the
    // authority for copper (it is the artwork); the mapping files in
    // design/mapping/ put ring centres 65-81 um further out and differ in
    // pitch by 0.4 um/ring. mm_sim checks the aim point against P2PadMap.hh
    // at startup and warns separately for a radial gap, an azimuthal gap, or
    // an aim point off the pad field.
    //
    // For any pad-multiplicity, charge-sharing or positional observable use
    // --beam-spread as well: a fixed point, even a correctly chosen one,
    // describes that point and not the pad cell.
    double p2_gun_x_mm       = 305.385; // ring 20 pad centre: r 349.2215 mm,
    double p2_gun_y_mm       = 169.399; //   phi 29.0174 deg
    // Radius [mm] of a uniform disc, transverse to the beam, over which the
    // impact point is scattered per event. 0 = pencil beam. Use >= one ring
    // pitch (11.4286 mm, gerber value) to average over the pad cell.
    double p2_beam_spread_mm = 0.0;
    // Beam tilt (P0.3). theta is measured from the wedge normal (+z), phi is
    // the azimuth of the tilt in the wedge plane: phi=0 tilts toward +x, 90
    // toward +y. The beam always passes through the aim point (gun-x, gun-y)
    // at the DRIFT MID-PLANE, so changing theta rotates the track about that
    // point instead of sliding the illuminated pads across the wedge.
    // Campaign scans theta in 0,10,20,30,40 deg (SIM_CAMPAIGN_PLAN 4.1).
    double p2_gun_theta_deg  = 0.0;
    double p2_gun_phi_deg    = 0.0;
    // Distance from the aim point back to the gun, along the beam direction.
    // Must clear the front window bulge at any theta; the default is checked
    // against the built geometry at run time.
    double p2_gun_standoff_mm = 50.0;

    // ── Readout copper (P0.19) ────────────────────────────────────────────
    // Build the F.Cu pad field as REAL copper — 1280 annular-sector pads on
    // the gerber-exact polar grid, 42 G4PVReplica rings — instead of a
    // density-scaled sheet. Costs 126 logical volumes and no extra memory.
    // --homogenized-readout turns it off; both layers are zoned by radial
    // band either way.
    //
    // NOTE: p2_fcu_coverage / p2_bcu_coverage above NO LONGER REACH THE
    // GEOMETRY. The readout copper is built from the per-radial-band table
    // in include/P2PadMap.hh, measured 2026-08-07 directly from the gerbers.
    // The 0.983 in particular is wrong — it came from an aperture bug in
    // scripts/gerber/analyze_cu_coverage.py that painted the 10 connector
    // footprints as ~150 mm-radius discs (docs/P2_GEOMETRY.md §3). The two
    // scalars are retained only because meta::GeometryDigest hashes them.
    bool p2_patterned_readout = true;

    // ── Spectrum sampling (kLSCalib / kBackScintCalib) ────────────────────
    // When non-empty, PrimaryGeneratorAction samples electron energies from
    // this CSV file (two columns: energy_MeV, probability) instead of using
    // the fixed 'energy' value above.
    std::string spectrum_file;

    // ── LS cell wall parameters (from Full_Geant geometry) ───────────────
    double cfrpThickness_mm  = 2.0;    // Structural CFRP wall [mm]
    double ls_inner_cfrp_um  = 600.0;  // Inner CFRP liner [µm]
    double ls_inner_al_um    = 40.0;   // Al liner [µm]
    double ls_thick_cm       = 2.0;    // LAB layer thickness [cm]

    // ── Back plastic scintillator (kBackScintCalib) ───────────────────────
    double backscint_u_cm     = 30.0;  // Back scint face: u-width [cm]
    double backscint_v_cm     = 20.0;  // Back scint face: v-height [cm]
    double backscint_thick_cm = 2.0;   // Back scint thickness [cm]
    double backscint_tape_um  = 200.0; // Outer black mylar tape [µm]
    double backscint_al_um    = 20.0;  // Al foil on scint surface [µm]

    // ── kLSCalib / kBackScintCalib: source-to-detector air gap ───────────
    double source_to_det_mm  = 100.0;  // Air gap from gun to detector front face [mm]
};
