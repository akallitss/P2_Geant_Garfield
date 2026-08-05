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

    // ── P2 wedge mode ─────────────────────────────────────────────────────
    // Values marked GUESS are not backed by design data — see docs/P2_MODEL.md
    // for the full assumptions table and what to confirm with the collaboration.
    double p2_drift_mm       = 3.0;    // drift gap; 3 mm confirmed baseline (2026-08-05).
                                       // Campaign will scan 1..4 mm via --drift-gap.
    double p2_amp_um         = 150.0;  // amplification gap — confirmed 150 um (2026-08-05)
    double p2_mesh_wire_um   = 19.0;   // woven SS mesh wire diameter — confirmed
    double p2_mesh_open_um   = 48.0;   // mesh opening (pitch = wire + opening = 67 um) —
                                       // confirmed "48x19" (Alexandra 2026-08-05); modelled
                                       // as a 2*wire slab of effective-density steel
    double p2_front_gap_mm   = 4.0;    // outer window -> first drift foil. Window->mesh
                                       // stays 8 mm total (frame STEP: ledge 3 -> top 8);
                                       // 1 mm now sits between the two drift foils.
    double p2_cath_gap_mm    = 1.0;    // gap between the two drift-cathode mylar foils
                                       // ("maybe 1 mm", Alexandra 2026-08-05)
    double p2_back_gap_mm    = 1.0;    // PCB back -> back window = carbon back-frame depth.
                                       // Normally 1 mm, up to 3 mm in this production
                                       // (Alexandra 2026-08-05); mylar glued on its rear face.
    double p2_bulge_front_mm = 10.0;   // overpressure sag — Hencky estimate at ~3 mbar,
    double p2_bulge_back_mm  = 5.0;    //   p in 1..10 mbar -> 7..15 mm; revisit
    double p2_window_um      = 40.0;   // outer containment mylar (MX17-like) — confirm
    double p2_cath_mylar_um  = 12.0;   // each drift-cathode mylar foil — thickness GUESS
    double p2_cath_al_um     = 0.1;    // aluminization of the drift-gas-side foil — GUESS
    double p2_fcu_coverage   = 0.983;  // Cu area fraction over the active area, measured
    double p2_bcu_coverage   = 0.174;  //   from the gerbers (B.Cu = signal lines, not a
                                       //   plane; radial 0.03->0.26). Layers modelled as
                                       //   full-thickness slabs of density-scaled copper.
                                       //   scripts/gerber/analyze_cu_coverage.py
    double p2_gun_x_mm       = 307.4;  // gun aim point: r=355 mm, phi=30 deg
    double p2_gun_y_mm       = 177.5;  //   (mid-active-area; apex/beam axis is x=y=0)

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
