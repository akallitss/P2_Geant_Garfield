// mm_sim.cc
// Micromegas + full-experiment GEANT4 simulation
// Author: Dylan Neff / n_TOF X17 group

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SimConfig.hh"

#include <iostream>
#include <string>
#include <ctime>

#ifdef USE_ROOT
#include "TROOT.h"
#endif

#include "GasMixtures.hh"

void PrintUsage() {
    std::cerr << "Usage: mm_sim [options] [macro_file]\n";
    std::cerr << "Options:\n";
    std::cerr << "  -g <gas>         Gas mixture name  (default: ArIso)\n";
    std::cerr << "                   --list-gases prints the full table with compositions\n";
    std::cerr << "  --list-gases     List available gas mixtures and exit\n";
    std::cerr << "  -p <particle>    gamma, neutron, electron, positron, proton, muon, muon+,\n";
    std::cerr << "                   pion, alpha, triton  (default: electron)\n";
    std::cerr << "  -e <energy>      Particle energy [MeV]  (default: 155, MESA beam)\n";
    std::cerr << "  -n <nevents>     Number of events  (default: 10000)\n";
    std::cerr << "  -o <output>      Output file base name  (default: mm_output)\n";
    std::cerr << "  -s <seed>        Random seed  (default: time-based)\n";
    std::cerr << "  -t <nthreads>    MT threads  (default: 1)\n";
    std::cerr << "  -m <mode>        p2 | vacuum | full | sr90 | sr90nomm | lscalib | backscintcalib\n";
    std::cerr << "                   (default: p2)\n";
    std::cerr << "P2 mode options (see docs/P2_MODEL.md for defaults' provenance):\n";
    std::cerr << "  --drift-gap <mm>    Drift gas, campaign scans 1..4  (default: 3.964 = 4 mm frame)\n";
    std::cerr << "  --amp-gap <um>      Amplification gap  (default: 150)\n";
    std::cerr << "  --back-gap <mm>     Carbon back-frame depth = back gas gap,\n";
    std::cerr << "                      1 normally, up to 3 this production  (default: 1)\n";
    std::cerr << "  --bulge-front <mm>  Front window overpressure sag  (default: 10)\n";
    std::cerr << "  --bulge-back <mm>   Back window overpressure sag  (default: 10)\n";
    std::cerr << "  --pillars <set>     Mesh-support pillars: saclay (V1 mask, det1-det4),\n";
    std::cerr << "                      cern (det5), none  (default: saclay)\n";
    std::cerr << "  --gun-x <mm>        Beam aim x [gerber coords]  (default: 302.49)\n";
    std::cerr << "  --gun-y <mm>        Beam aim y [gerber coords]  (default: 174.64)\n";
    std::cerr << "                      Default is a pad-RING CENTRE (r=349.286, phi=30);\n";
    std::cerr << "                      a boundary aim point biases pad multiplicity\n";
    std::cerr << "  --beam-spread <mm>  Scatter the impact point over a disc of this radius\n";
    std::cerr << "                      (0 = pencil beam). REQUIRED for pad-multiplicity or\n";
    std::cerr << "                      positional observables; use >= 11.43 (one ring pitch)\n";
    std::cerr << "  --gun-theta <deg>   Beam tilt from the wedge normal  (default: 0)\n";
    std::cerr << "  --gun-phi <deg>     Azimuth of the tilt, 0 = toward +x  (default: 0)\n";
    std::cerr << "                      The beam pivots about (gun-x, gun-y) at the drift\n";
    std::cerr << "                      mid-plane, so theta does not move the illuminated pads\n";
    std::cerr << "  --gun-standoff <mm> Aim point to gun, along the beam  (default: 50,\n";
    std::cerr << "                      raised automatically to clear the window bulge)\n";
    std::cerr << "  --w-cf4 <eV>        Override CF4's W-value to bracket the 35-52 eV spread;\n";
    std::cerr << "                      re-derives the mixture W  (default: table value, 34)\n";
    std::cerr << "  --homogenized-readout  Build the F.Cu pad field as a density-scaled\n";
    std::cerr << "                      sheet instead of 1280 real annular-sector pads.\n";
    std::cerr << "                      Both layers stay zoned by radial band either way.\n";
    std::cerr << "  --spectrum <csv> Sample energies from Sr-90/Y-90 CSV (lscalib/backscintcalib)\n";
    std::cerr << "  --src-dist <mm>  Source-to-detector air gap [mm] (default: 100)\n";
    std::cerr << "  -a <mm>          Al shielding [mm], vacuum mode only  (default: 0)\n";
    std::cerr << "  -c <mm>          CFRP wall thickness [mm] for LS cells, full mode only  (default: 1.5)\n";
    std::cerr << "  --skip-empty     Write only events with >=1 gas cluster (photon runs;\n";
    std::cerr << "                   total thrown is still reported for normalization)\n";
    std::cerr << "  -v               Verbose output\n";
    std::cerr << "  -h               Print this help\n";
}

int main(int argc, char** argv) {
    SimConfig config;
    config.gas            = "ArIso";        // ASSUMED P2 gas — confirm (docs/P2_MODEL.md)
    config.particle       = "electron";
    config.energy         = 155.0 * MeV;    // MESA beam energy
    config.nEvents        = 10000;
    config.outFile        = "mm_output";
    config.seed           = static_cast<long>(std::time(nullptr));
    config.nThreads       = 1;
    config.verbose        = false;
    config.alThickness_mm = 0.0;
    config.mode           = SimMode::kP2Wedge;
    config.cfrpThickness_mm = 2.0;  // updated to match Full_Geant

    std::string macroFile = "";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "-h")              { PrintUsage(); return 0; }
        else if (arg == "-g" && i+1<argc)  config.gas              = argv[++i];
        else if (arg == "-p" && i+1<argc)  config.particle         = argv[++i];
        else if (arg == "-e" && i+1<argc)  config.energy           = std::stod(argv[++i]) * MeV;
        else if (arg == "-n" && i+1<argc)  config.nEvents          = std::stoi(argv[++i]);
        else if (arg == "-o" && i+1<argc)  config.outFile          = argv[++i];
        else if (arg == "-s" && i+1<argc)  config.seed             = std::stol(argv[++i]);
        else if (arg == "-t" && i+1<argc)  config.nThreads         = std::stoi(argv[++i]);
        else if (arg == "-a" && i+1<argc)  config.alThickness_mm   = std::stod(argv[++i]);
        else if (arg == "-c" && i+1<argc)  config.cfrpThickness_mm = std::stod(argv[++i]);
        else if (arg == "-v")              config.verbose           = true;
        else if (arg == "--skip-empty")    config.skipEmpty         = true;
        else if (arg == "-m" && i+1<argc) {
            std::string mode = argv[++i];
            if      (mode == "p2")     config.mode = SimMode::kP2Wedge;
            else if (mode == "vacuum") config.mode = SimMode::kVacuum;
            else if (mode == "full")   config.mode = SimMode::kFullExperiment;
            else if (mode == "sr90")    config.mode = SimMode::kSr90Calibration;
            else if (mode == "sr90nomm") config.mode = SimMode::kSr90NoMM;
            else if (mode == "lscalib")       config.mode = SimMode::kLSCalib;
            else if (mode == "backscintcalib") config.mode = SimMode::kBackScintCalib;
            else { std::cerr << "Unknown mode: " << mode << "\n"; return 1; }
        }
        else if (arg == "--spectrum" && i+1<argc) config.spectrum_file  = argv[++i];
        else if (arg == "--src-dist" && i+1<argc) config.source_to_det_mm = std::stod(argv[++i]);
        else if (arg == "--drift-gap"   && i+1<argc) config.p2_drift_mm       = std::stod(argv[++i]);
        else if (arg == "--amp-gap"     && i+1<argc) config.p2_amp_um         = std::stod(argv[++i]);
        else if (arg == "--back-gap"    && i+1<argc) config.p2_back_gap_mm    = std::stod(argv[++i]);
        else if (arg == "--bulge-front" && i+1<argc) config.p2_bulge_front_mm = std::stod(argv[++i]);
        else if (arg == "--bulge-back"  && i+1<argc) config.p2_bulge_back_mm  = std::stod(argv[++i]);
        else if (arg == "--pillars"     && i+1<argc) {
            config.p2_pillars = argv[++i];
            if (config.p2_pillars != "saclay" && config.p2_pillars != "cern" &&
                config.p2_pillars != "none") {
                std::cerr << "Unknown --pillars set: " << config.p2_pillars
                          << " (saclay | cern | none)\n";
                return 1;
            }
        }
        else if (arg == "--gun-x"       && i+1<argc) config.p2_gun_x_mm       = std::stod(argv[++i]);
        else if (arg == "--gun-y"       && i+1<argc) config.p2_gun_y_mm       = std::stod(argv[++i]);
        else if (arg == "--gun-theta"   && i+1<argc) config.p2_gun_theta_deg  = std::stod(argv[++i]);
        else if (arg == "--gun-phi"     && i+1<argc) config.p2_gun_phi_deg    = std::stod(argv[++i]);
        else if (arg == "--gun-standoff"&& i+1<argc) config.p2_gun_standoff_mm= std::stod(argv[++i]);
        else if (arg == "--beam-spread" && i+1<argc) config.p2_beam_spread_mm = std::stod(argv[++i]);
        else if (arg == "--w-cf4"       && i+1<argc) config.w_cf4_eV          = std::stod(argv[++i]);
        else if (arg == "--homogenized-readout")     config.p2_patterned_readout = false;
        else if (arg == "--list-gases") { std::cout << gas::ListMixtures(); return 0; }
        else if (arg[0] != '-') macroFile = arg;
        else { std::cerr << "Unknown option: " << arg << "\n"; PrintUsage(); return 1; }
    }

    // Validate here, not deep inside the run manager. An unknown gas or
    // particle used to surface as an uncaught std::runtime_error thrown from
    // PrimaryGeneratorAction's constructor on a worker thread, i.e. a core
    // dump after the full geometry had already been built and printed.
    if (!gas::Exists(config.gas)) {
        std::cerr << "Unknown gas: " << config.gas << "\n" << gas::ListMixtures();
        return 1;
    }
    if (config.p2_gun_theta_deg < 0.0 || config.p2_gun_theta_deg >= 90.0) {
        std::cerr << "--gun-theta must be in [0, 90) deg, got "
                  << config.p2_gun_theta_deg << "\n";
        return 1;
    }

#ifdef USE_ROOT
    // Each Geant4 worker thread builds its own TFile + TTrees in
    // RunAction::BeginOfRunAction. ROOT's global type/class registry is not
    // thread-safe by default, so concurrent TTree construction races inside
    // TROOT::GetListOfTypes and aborts. This must be called before any
    // worker starts. (Found 2026-08-05: -t 4 and -t 8 core-dumped while
    // -t 1 was fine.)
    ROOT::EnableThreadSafety();
#endif

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(config.seed);

    std::cout << "=== Micromegas Simulation ===" << "\n";
    std::cout << "  Geant4 version : " << G4Version << "\n";
    std::string modeStr = (config.mode == SimMode::kP2Wedge         ? "p2-wedge"            :
                           config.mode == SimMode::kFullExperiment  ? "full-experiment"       :
                           config.mode == SimMode::kSr90Calibration ? "sr90-calibration"    :
                           config.mode == SimMode::kSr90NoMM        ? "sr90-no-mm"          :
                           config.mode == SimMode::kLSCalib         ? "ls-calibration"      :
                           config.mode == SimMode::kBackScintCalib  ? "backscint-calibration": "vacuum");
    std::cout << "  Mode           : " << modeStr << "\n";
    const gas::Mixture& mix = gas::Find(config.gas);
    std::cout << "  Gas            : " << config.gas << "  (" << mix.label
              << ", W = " << gas::MixtureWValue(mix, config.w_cf4_eV)
              << " eV)\n";
    if (!mix.note.empty())
        std::cout << "                   note: " << mix.note << "\n";
    std::cout << "  Particle       : " << config.particle << "\n";
    std::cout << "  Energy         : " << config.energy/MeV << " MeV" << "\n";
    std::cout << "  Events         : " << config.nEvents << "\n";
    std::cout << "  Output         : " << config.outFile << "\n";
    std::cout << "  Seed           : " << config.seed << "\n";
    std::cout << "  Threads        : " << config.nThreads << "\n";
    if (config.mode == SimMode::kP2Wedge) {
        std::cout << "  Drift gap      : " << config.p2_drift_mm << " mm\n";
        std::cout << "  Amp gap        : " << config.p2_amp_um << " um\n";
        std::cout << "  Pillars        : " << config.p2_pillars << "\n";
        std::cout << "  Gun aim (x,y)  : (" << config.p2_gun_x_mm << ", "
                  << config.p2_gun_y_mm << ") mm\n";
        std::cout << "  Gun angle      : theta " << config.p2_gun_theta_deg
                  << " deg, phi " << config.p2_gun_phi_deg << " deg\n";
        std::cout << "  Beam spread    : " << config.p2_beam_spread_mm
                  << " mm" << (config.p2_beam_spread_mm <= 0.0
                               ? "  [pencil beam - not valid for pad-level observables]"
                               : "") << "\n";
    } else if (config.mode == SimMode::kVacuum)
        std::cout << "  Al shielding   : " << config.alThickness_mm << " mm\n";
    else
        std::cout << "  CFRP thickness : " << config.cfrpThickness_mm << " mm\n";
    std::cout << "=============================" << "\n";

#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(config.nThreads);
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // DetectorConstruction must outlive ActionInitialization so that
    // PrimaryGeneratorAction can read GetHe3GasCenterZ() after Construct().
    auto* detCon = new DetectorConstruction(config);
    runManager->SetUserInitialization(detCon);
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserInitialization(new ActionInitialization(config, detCon));

    runManager->Initialize();

    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (macroFile.empty()) {
        UI->ApplyCommand("/run/beamOn " + std::to_string(config.nEvents));
    } else {
        UI->ApplyCommand("/control/execute " + macroFile);
    }

    delete runManager;
    std::cout << "Simulation complete. Output: " << config.outFile << "\n";
    return 0;
}
