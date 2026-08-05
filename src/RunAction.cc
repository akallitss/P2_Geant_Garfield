// RunAction.cc
// Writes ROOT TTrees (or CSV fallback) per worker thread.

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstring>

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

struct RunAction::Impl {
#ifdef USE_ROOT
    TFile* rootFile  = nullptr;
    TTree* evtTree   = nullptr;
    TTree* clusTree  = nullptr;
    TTree* volTree   = nullptr;

    // EventTree — P2 provenance / classification (added 2026-08-05)
    Int_t    br_class;
    Char_t   br_convProcess[32], br_convVolume[32], br_convParentVol[32];
    Double_t br_convZ, br_convEnergy, br_convParentE;
    Char_t   br_firstProcess[32], br_firstVolume[32];
    Double_t br_firstX, br_firstY, br_firstZ, br_firstE;
    Int_t    br_primaryPDG;
    Double_t br_primaryCharge;
    Double_t br_edepFrontGas, br_edepCathGas, br_edepBackGas, br_edepWindowGas;
    Double_t br_edepCathMylar, br_edepCathAl, br_edepWindowMylar;
    Double_t br_edepMeshP2, br_edepPadCu, br_edepFR4P2, br_edepCuB, br_edepFrame;

    // VolumeTree — one row per (event, volume) with non-zero edep
    Int_t    vb_eventID, vb_nSteps;
    Double_t vb_edep;
    Char_t   vb_volume[64];

    // EventTree — always present
    Int_t    br_eventID;
    Double_t br_edepDrift, br_edepAmp;
    Int_t    br_nPrimDrift, br_nPrimAmp;
    Int_t    br_nClusDrift, br_nClusAmp;
    Bool_t   br_primInDrift, br_primInAmp;

    // EventTree — full/sr90 modes
    Double_t br_edepHe3Gas, br_edepResistPaste;
    Double_t br_edepMylar, br_edepCathode, br_edepMicromesh;
    Double_t br_edepPCB;
    Double_t br_edepPCBKapton, br_edepPCBCu, br_edepPCBFR4,
             br_edepPCBRohacell, br_edepPCBAlFoil;
    Double_t br_edepScintWall, br_edepScintTape, br_edepScintAlFoil;
    Double_t br_edepLS1, br_edepLS2;
    Bool_t   br_primInHe3Gas, br_primInPCB, br_primInScintWall;
    Bool_t   br_primInLS1, br_primInLS2;
    Double_t br_edepLSCFRP;
    Bool_t   br_primInLSCFRP5;

    // EventTree — kLSCalib / kBackScintCalib
    Double_t br_edepBackScint;
    Double_t br_primaryKE;
    Bool_t   br_primInBackScint;

    // ClusterTree
    Int_t    cb_eventID, cb_trackID, cb_parentID, cb_nPrimary;
    Double_t cb_x, cb_y, cb_z, cb_edep, cb_ke;
    Char_t   cb_volume[32];
    Char_t   cb_particle[32];
    // ClusterTree — provenance (added 2026-08-05)
    Char_t   cb_creatorProc[32], cb_originProc[32], cb_originVol[32],
             cb_parentVol[32];
    Double_t cb_originX, cb_originY, cb_originZ, cb_originE,
             cb_parentE, cb_time;
    Int_t    cb_ancestorID;
#else
    std::ofstream evtFile;
    std::ofstream clusFile;
    std::ofstream volFile;
#endif
};

RunAction::RunAction(const SimConfig& cfg, bool isMaster)
    : G4UserRunAction(), fConfig(cfg), fIsMaster(isMaster),
      fImpl(std::make_unique<Impl>()) {}

RunAction::~RunAction() = default;

// ============================================================
void RunAction::BeginOfRunAction(const G4Run*) {
    fTotalEvents = 0;
    fSumEdepDrift = fSumNprimDrift = 0.0;
    fSumEdepAmp  = fSumNprimAmp  = 0.0;
    fWrittenEvents = fSkippedEmpty = 0;

    if (fIsMaster) return;

    std::ostringstream ss;
    ss << fConfig.outFile;
    int tid = G4Threading::G4GetThreadId();
    if (tid >= 0) ss << "_t" << tid;

    bool isFull    = (fConfig.mode == SimMode::kFullExperiment  ||
                      fConfig.mode == SimMode::kSr90Calibration ||
                      fConfig.mode == SimMode::kSr90NoMM);
    bool isCalib   = (fConfig.mode == SimMode::kLSCalib ||
                      fConfig.mode == SimMode::kBackScintCalib);

#ifdef USE_ROOT
    std::string fname = ss.str() + ".root";
    fImpl->rootFile = TFile::Open(fname.c_str(), "RECREATE");
    if (!fImpl->rootFile || fImpl->rootFile->IsZombie()) {
        G4cerr << "ERROR: Cannot open ROOT file " << fname << G4endl;
        return;
    }

    std::string modeStr = (fConfig.mode == SimMode::kP2Wedge         ? "p2"           :
                           fConfig.mode == SimMode::kFullExperiment  ? "full"         :
                           fConfig.mode == SimMode::kSr90Calibration ? "sr90"         :
                           fConfig.mode == SimMode::kSr90NoMM        ? "sr90nomm"     :
                           fConfig.mode == SimMode::kLSCalib         ? "lscalib"      :
                           fConfig.mode == SimMode::kBackScintCalib  ? "backscintcalib": "vacuum");
    std::string tag = "mode=" + modeStr
                    + "_gas=" + fConfig.gas
                    + "_particle=" + fConfig.particle
                    + "_E=" + std::to_string(fConfig.energy/MeV) + "MeV";

    fImpl->evtTree = new TTree("EventTree", "Per-event summary");
    fImpl->evtTree->SetTitle(tag.c_str());

    fImpl->evtTree->Branch("eventID",    &fImpl->br_eventID);
    fImpl->evtTree->Branch("edepDrift",  &fImpl->br_edepDrift);
    fImpl->evtTree->Branch("edepAmp",    &fImpl->br_edepAmp);
    fImpl->evtTree->Branch("nPrimDrift", &fImpl->br_nPrimDrift);
    fImpl->evtTree->Branch("nPrimAmp",   &fImpl->br_nPrimAmp);
    fImpl->evtTree->Branch("nClusDrift", &fImpl->br_nClusDrift);
    fImpl->evtTree->Branch("nClusAmp",   &fImpl->br_nClusAmp);
    fImpl->evtTree->Branch("primInDrift",&fImpl->br_primInDrift);
    fImpl->evtTree->Branch("primInAmp",  &fImpl->br_primInAmp);

    // ── Provenance / classification — always written ─────────────────────
    // "interactionClass" is the branch the headline plot cuts on:
    //   EventTree->Draw("edepDrift", "interactionClass==2")  etc.
    // See InteractionClass in include/EventData.hh for the codes.
    fImpl->evtTree->Branch("interactionClass", &fImpl->br_class);
    fImpl->evtTree->Branch("convProcess", fImpl->br_convProcess, "convProcess[32]/C");
    fImpl->evtTree->Branch("convVolume",  fImpl->br_convVolume,  "convVolume[32]/C");
    fImpl->evtTree->Branch("convParentBirthVolume", fImpl->br_convParentVol,
                           "convParentBirthVolume[32]/C");
    fImpl->evtTree->Branch("convZ",              &fImpl->br_convZ);
    fImpl->evtTree->Branch("convEnergy",         &fImpl->br_convEnergy);
    fImpl->evtTree->Branch("convParentEnergy",   &fImpl->br_convParentE);
    fImpl->evtTree->Branch("firstIntProcess", fImpl->br_firstProcess,
                           "firstIntProcess[32]/C");
    fImpl->evtTree->Branch("firstIntVolume",  fImpl->br_firstVolume,
                           "firstIntVolume[32]/C");
    fImpl->evtTree->Branch("firstIntX",      &fImpl->br_firstX);
    fImpl->evtTree->Branch("firstIntY",      &fImpl->br_firstY);
    fImpl->evtTree->Branch("firstIntZ",      &fImpl->br_firstZ);
    fImpl->evtTree->Branch("firstIntEnergy", &fImpl->br_firstE);
    fImpl->evtTree->Branch("primaryPDG",     &fImpl->br_primaryPDG);
    fImpl->evtTree->Branch("primaryCharge",  &fImpl->br_primaryCharge);
    fImpl->evtTree->Branch("primaryKE_MeV",  &fImpl->br_primaryKE);

    // Per-layer sums for the P2 stack — the conversion budget.
    // The *Gas ones are regions where a deposit produces NO signal; keeping
    // them is what separates "converted somewhere useless" from "no
    // conversion at all".
    fImpl->evtTree->Branch("edepFrontGas",    &fImpl->br_edepFrontGas);
    fImpl->evtTree->Branch("edepCathGas",     &fImpl->br_edepCathGas);
    fImpl->evtTree->Branch("edepBackGas",     &fImpl->br_edepBackGas);
    fImpl->evtTree->Branch("edepWindowGas",   &fImpl->br_edepWindowGas);
    fImpl->evtTree->Branch("edepCathMylar",   &fImpl->br_edepCathMylar);
    fImpl->evtTree->Branch("edepCathAl",      &fImpl->br_edepCathAl);
    fImpl->evtTree->Branch("edepWindowMylar", &fImpl->br_edepWindowMylar);
    fImpl->evtTree->Branch("edepMesh",        &fImpl->br_edepMeshP2);
    fImpl->evtTree->Branch("edepPadCu",       &fImpl->br_edepPadCu);
    fImpl->evtTree->Branch("edepFR4",         &fImpl->br_edepFR4P2);
    fImpl->evtTree->Branch("edepCuB",         &fImpl->br_edepCuB);
    fImpl->evtTree->Branch("edepFrame",       &fImpl->br_edepFrame);

    if (isFull) {
        fImpl->evtTree->Branch("edepHe3Gas",      &fImpl->br_edepHe3Gas);
        fImpl->evtTree->Branch("edepResistPaste", &fImpl->br_edepResistPaste);
        fImpl->evtTree->Branch("edepMylar",       &fImpl->br_edepMylar);
        fImpl->evtTree->Branch("edepCathode",     &fImpl->br_edepCathode);
        fImpl->evtTree->Branch("edepMicromesh",   &fImpl->br_edepMicromesh);
        fImpl->evtTree->Branch("edepPCB",         &fImpl->br_edepPCB);
        fImpl->evtTree->Branch("edepPCBKapton",   &fImpl->br_edepPCBKapton);
        fImpl->evtTree->Branch("edepPCBCu",       &fImpl->br_edepPCBCu);
        fImpl->evtTree->Branch("edepPCBFR4",      &fImpl->br_edepPCBFR4);
        fImpl->evtTree->Branch("edepPCBRohacell", &fImpl->br_edepPCBRohacell);
        fImpl->evtTree->Branch("edepPCBAlFoil",   &fImpl->br_edepPCBAlFoil);
        fImpl->evtTree->Branch("edepScintWall",   &fImpl->br_edepScintWall);
        fImpl->evtTree->Branch("edepScintTape",   &fImpl->br_edepScintTape);
        fImpl->evtTree->Branch("edepScintAlFoil", &fImpl->br_edepScintAlFoil);
        fImpl->evtTree->Branch("edepLS1",          &fImpl->br_edepLS1);
        fImpl->evtTree->Branch("edepLS2",          &fImpl->br_edepLS2);
        fImpl->evtTree->Branch("primInHe3Gas",     &fImpl->br_primInHe3Gas);
        fImpl->evtTree->Branch("primInPCB",        &fImpl->br_primInPCB);
        fImpl->evtTree->Branch("primInScintWall",  &fImpl->br_primInScintWall);
        fImpl->evtTree->Branch("primInLS1",        &fImpl->br_primInLS1);
        fImpl->evtTree->Branch("primInLS2",        &fImpl->br_primInLS2);
        fImpl->evtTree->Branch("edepLSCFRP",       &fImpl->br_edepLSCFRP);
        fImpl->evtTree->Branch("primInLSCFRP5",    &fImpl->br_primInLSCFRP5);
    }

    if (isCalib) {
        // primaryKE_MeV is now branched unconditionally above.
        fImpl->evtTree->Branch("edepLS1",          &fImpl->br_edepLS1);
        fImpl->evtTree->Branch("edepLSCFRP",       &fImpl->br_edepLSCFRP);
        fImpl->evtTree->Branch("primInLS1",        &fImpl->br_primInLS1);
        fImpl->evtTree->Branch("edepBackScint",    &fImpl->br_edepBackScint);
        fImpl->evtTree->Branch("primInBackScint",  &fImpl->br_primInBackScint);
    }

    fImpl->clusTree = new TTree("ClusterTree", "Per-cluster ionization detail");
    fImpl->clusTree->SetTitle(tag.c_str());
    fImpl->clusTree->Branch("eventID",  &fImpl->cb_eventID);
    fImpl->clusTree->Branch("trackID",  &fImpl->cb_trackID);
    fImpl->clusTree->Branch("parentID", &fImpl->cb_parentID);
    fImpl->clusTree->Branch("x",        &fImpl->cb_x);
    fImpl->clusTree->Branch("y",        &fImpl->cb_y);
    fImpl->clusTree->Branch("z",        &fImpl->cb_z);
    fImpl->clusTree->Branch("edep",     &fImpl->cb_edep);
    fImpl->clusTree->Branch("nPrimary", &fImpl->cb_nPrimary);
    fImpl->clusTree->Branch("ke",       &fImpl->cb_ke);
    fImpl->clusTree->Branch("volume",   fImpl->cb_volume,  "volume[32]/C");
    fImpl->clusTree->Branch("particle", fImpl->cb_particle,"particle[32]/C");
    fImpl->clusTree->Branch("time",     &fImpl->cb_time);
    // Provenance: creatorProcess is this track's own process; origin* is the
    // interaction that started the chain (see include/TrackOrigin.hh).
    fImpl->clusTree->Branch("creatorProcess", fImpl->cb_creatorProc,
                            "creatorProcess[32]/C");
    fImpl->clusTree->Branch("originProcess",  fImpl->cb_originProc,
                            "originProcess[32]/C");
    fImpl->clusTree->Branch("originVolume",   fImpl->cb_originVol,
                            "originVolume[32]/C");
    fImpl->clusTree->Branch("parentBirthVolume", fImpl->cb_parentVol,
                            "parentBirthVolume[32]/C");
    fImpl->clusTree->Branch("originX",      &fImpl->cb_originX);
    fImpl->clusTree->Branch("originY",      &fImpl->cb_originY);
    fImpl->clusTree->Branch("originZ",      &fImpl->cb_originZ);
    fImpl->clusTree->Branch("originEnergy", &fImpl->cb_originE);
    fImpl->clusTree->Branch("parentEnergy", &fImpl->cb_parentE);
    fImpl->clusTree->Branch("ancestorID",   &fImpl->cb_ancestorID);

    // VolumeTree: every volume that received energy, named, one row each.
    fImpl->volTree = new TTree("VolumeTree", "Per-event per-volume edep");
    fImpl->volTree->SetTitle(tag.c_str());
    fImpl->volTree->Branch("eventID", &fImpl->vb_eventID);
    fImpl->volTree->Branch("volume",  fImpl->vb_volume, "volume[64]/C");
    fImpl->volTree->Branch("edep",    &fImpl->vb_edep);
    fImpl->volTree->Branch("nSteps",  &fImpl->vb_nSteps);

    G4cout << "RunAction: Opened " << fname << G4endl;

#else
    std::string evtFname  = ss.str() + "_events.csv";
    std::string clusFname = ss.str() + "_clusters.csv";

    fImpl->evtFile.open(evtFname);
    fImpl->evtFile << "eventID,edepDrift_eV,edepAmp_eV,nPrimDrift,nPrimAmp,"
                      "nClusDrift,nClusAmp,primInDrift,primInAmp"
                      ",interactionClass,convProcess,convVolume"
                      ",convParentBirthVolume,convZ_mm,convEnergy_MeV"
                      ",convParentEnergy_MeV"
                      ",firstIntProcess,firstIntVolume"
                      ",firstIntX_mm,firstIntY_mm,firstIntZ_mm,firstIntE_MeV"
                      ",primaryPDG,primaryCharge,primaryKE_MeV"
                      ",edepFrontGas_eV,edepCathGas_eV,edepBackGas_eV"
                      ",edepWindowGas_eV,edepCathMylar_eV,edepCathAl_eV"
                      ",edepWindowMylar_eV,edepMesh_eV,edepPadCu_eV"
                      ",edepFR4_eV,edepCuB_eV,edepFrame_eV";
    if (isFull) {
        fImpl->evtFile << ",edepHe3Gas_eV,edepResistPaste_eV"
                          ",edepMylar_eV,edepCathode_eV,edepMicromesh_eV"
                          ",edepPCB_eV,edepPCBKapton_eV,edepPCBCu_eV"
                          ",edepPCBFR4_eV,edepPCBRohacell_eV,edepPCBAlFoil_eV"
                          ",edepScintWall_eV,edepScintTape_eV,edepScintAlFoil_eV"
                          ",edepLS1_eV,edepLS2_eV"
                          ",primInHe3Gas,primInPCB,primInScintWall"
                          ",primInLS1,primInLS2"
                          ",edepLSCFRP_eV,primInLSCFRP5";
    }
    if (isCalib) {
        fImpl->evtFile << ",edepLS1_eV,edepLSCFRP_eV,primInLS1"
                          ",edepBackScint_eV,primInBackScint";
    }
    fImpl->evtFile << "\n";

    fImpl->clusFile.open(clusFname);
    fImpl->clusFile << "eventID,trackID,parentID,x_mm,y_mm,z_mm,"
                       "edep_eV,nPrimary,ke_MeV,volume,particle,time_ns,"
                       "creatorProcess,originProcess,originVolume,"
                       "parentBirthVolume,originX_mm,originY_mm,originZ_mm,"
                       "originEnergy_MeV,parentEnergy_MeV,ancestorID\n";

    fImpl->volFile.open(ss.str() + "_volumes.csv");
    fImpl->volFile << "eventID,volume,edep_eV,nSteps\n";

    G4cout << "RunAction: Writing " << evtFname << G4endl;
#endif
}

// ============================================================
void RunAction::EndOfRunAction(const G4Run* run) {
    if (!fIsMaster) {
#ifdef USE_ROOT
        if (fImpl->rootFile) {
            fImpl->rootFile->Write();
            fImpl->rootFile->Close();
            delete fImpl->rootFile;
            fImpl->rootFile = nullptr;
        }
#else
        fImpl->evtFile.close();
        fImpl->clusFile.close();
        fImpl->volFile.close();
#endif
        if (fConfig.skipEmpty)
            G4cout << "  [worker] thrown=" << fTotalEvents
                   << " written=" << fWrittenEvents
                   << " skippedEmpty=" << fSkippedEmpty << G4endl;
    }

    if (fIsMaster || !G4Threading::IsMultithreadedApplication()) {
        G4int nev = run->GetNumberOfEvent();
        std::string modeStr = (fConfig.mode == SimMode::kP2Wedge         ? "p2-wedge"            :
                               fConfig.mode == SimMode::kFullExperiment  ? "full-experiment"       :
                               fConfig.mode == SimMode::kSr90Calibration ? "sr90-calibration"    :
                               fConfig.mode == SimMode::kSr90NoMM        ? "sr90-no-mm"          :
                               fConfig.mode == SimMode::kLSCalib         ? "ls-calibration"      :
                               fConfig.mode == SimMode::kBackScintCalib  ? "backscint-calibration": "vacuum");
        G4cout << "\n========= Run Summary =========" << G4endl;
        G4cout << "  Mode   : " << modeStr         << G4endl;
        G4cout << "  Gas    : " << fConfig.gas      << G4endl;
        G4cout << "  Events : " << nev              << G4endl;
        if (fTotalEvents > 0) {
            G4cout << std::fixed << std::setprecision(2);
            G4cout << "  <Edep_drift>  : " << fSumEdepDrift  / fTotalEvents << " eV" << G4endl;
            G4cout << "  <Edep_amp>    : " << fSumEdepAmp    / fTotalEvents << " eV" << G4endl;
        }
        // Normalization bookkeeping: with --skip-empty the trees hold fewer
        // rows than events thrown, and every per-photon probability must be
        // divided by the THROWN count. nev is the merged total and is
        // correct on the master; the per-worker written/skipped counts are
        // printed by each worker below (they do not merge automatically).
        if (fConfig.skipEmpty)
            G4cout << "  Thrown (norm) : " << nev
                   << "   [--skip-empty: divide rates by this]" << G4endl;
        G4cout << "===============================" << G4endl;
    }
}

// ============================================================
void RunAction::RecordEvent(const EventData& data) {
    fTotalEvents++;
    fSumEdepDrift  += data.edepDrift;
    fSumNprimDrift += data.nPrimaryDrift;
    fSumEdepAmp    += data.edepAmp;
    fSumNprimAmp   += data.nPrimaryAmp;

    // --skip-empty: photon runs are >=99.9 % empty. Accumulators above have
    // already counted this event, so normalization stays correct; we simply
    // do not write a row. fSkippedEmpty is reported at end of run.
    if (fConfig.skipEmpty &&
        data.driftClusters.empty() && data.ampClusters.empty()) {
        fSkippedEmpty++;
        return;
    }
    fWrittenEvents++;

    bool isFull    = (fConfig.mode == SimMode::kFullExperiment  ||
                      fConfig.mode == SimMode::kSr90Calibration ||
                      fConfig.mode == SimMode::kSr90NoMM);
    bool isCalib   = (fConfig.mode == SimMode::kLSCalib ||
                      fConfig.mode == SimMode::kBackScintCalib);

#ifdef USE_ROOT
    if (!fImpl->evtTree) return;

    fImpl->br_eventID    = data.eventID;
    fImpl->br_edepDrift  = data.edepDrift;
    fImpl->br_edepAmp    = data.edepAmp;
    fImpl->br_nPrimDrift = data.nPrimaryDrift;
    fImpl->br_nPrimAmp   = data.nPrimaryAmp;
    fImpl->br_nClusDrift = static_cast<int>(data.driftClusters.size());
    fImpl->br_nClusAmp   = static_cast<int>(data.ampClusters.size());
    fImpl->br_primInDrift = data.primaryInDrift;
    fImpl->br_primInAmp   = data.primaryInAmp;

    auto setStr = [](Char_t* dst, size_t n, const std::string& s) {
        std::strncpy(dst, s.c_str(), n - 1);
        dst[n - 1] = '\0';
    };

    fImpl->br_class = static_cast<Int_t>(data.interactionClass);
    setStr(fImpl->br_convProcess, 32, data.convProcess);
    setStr(fImpl->br_convVolume,  32, data.convVolume);
    setStr(fImpl->br_convParentVol, 32, data.convParentBirthVolume);
    fImpl->br_convZ        = data.convZ;
    fImpl->br_convEnergy   = data.convEnergy;
    fImpl->br_convParentE  = data.convParentEnergy;
    setStr(fImpl->br_firstProcess, 32, data.firstIntProcess);
    setStr(fImpl->br_firstVolume,  32, data.firstIntVolume);
    fImpl->br_firstX = data.firstIntX;
    fImpl->br_firstY = data.firstIntY;
    fImpl->br_firstZ = data.firstIntZ;
    fImpl->br_firstE = data.firstIntEnergy;
    fImpl->br_primaryPDG    = data.primaryPDG;
    fImpl->br_primaryCharge = data.primaryCharge;
    fImpl->br_primaryKE     = data.primaryKE_MeV;

    fImpl->br_edepFrontGas    = data.edepFrontGas;
    fImpl->br_edepCathGas     = data.edepCathGas;
    fImpl->br_edepBackGas     = data.edepBackGas;
    fImpl->br_edepWindowGas   = data.edepWindowGas;
    fImpl->br_edepCathMylar   = data.edepCathMylar;
    fImpl->br_edepCathAl      = data.edepCathAl;
    fImpl->br_edepWindowMylar = data.edepWindowMylar;
    fImpl->br_edepMeshP2      = data.edepMeshP2;
    fImpl->br_edepPadCu       = data.edepPadCu;
    fImpl->br_edepFR4P2       = data.edepFR4P2;
    fImpl->br_edepCuB         = data.edepCuB;
    fImpl->br_edepFrame       = data.edepFrame;

    if (isFull) {
        fImpl->br_edepHe3Gas      = data.edepHe3Gas;
        fImpl->br_edepResistPaste = data.edepResistPaste;
        fImpl->br_edepMylar       = data.edepMylar;
        fImpl->br_edepCathode     = data.edepCathode;
        fImpl->br_edepMicromesh   = data.edepMicromesh;
        fImpl->br_edepPCB         = data.edepPCB;
        fImpl->br_edepPCBKapton   = data.edepPCBKapton;
        fImpl->br_edepPCBCu       = data.edepPCBCu;
        fImpl->br_edepPCBFR4      = data.edepPCBFR4;
        fImpl->br_edepPCBRohacell = data.edepPCBRohacell;
        fImpl->br_edepPCBAlFoil   = data.edepPCBAlFoil;
        fImpl->br_edepScintWall   = data.edepScintWall;
        fImpl->br_edepScintTape   = data.edepScintTape;
        fImpl->br_edepScintAlFoil = data.edepScintAlFoil;
        fImpl->br_edepLS1         = data.edepLS1;
        fImpl->br_edepLS2         = data.edepLS2;
        fImpl->br_primInHe3Gas    = data.primInHe3Gas;
        fImpl->br_primInPCB       = data.primInPCB;
        fImpl->br_primInScintWall = data.primInScintWall;
        fImpl->br_primInLS1       = data.primInLS1;
        fImpl->br_primInLS2       = data.primInLS2;
        fImpl->br_edepLSCFRP      = data.edepLSCFRP;
        fImpl->br_primInLSCFRP5   = data.primInLSCFRP5;
    }

    if (isCalib) {
        fImpl->br_edepLS1         = data.edepLS1;
        fImpl->br_edepLSCFRP      = data.edepLSCFRP;
        fImpl->br_primInLS1       = data.primInLS1;
        fImpl->br_edepBackScint   = data.edepBackScint;
        fImpl->br_primInBackScint = data.primInBackScint;
    }

    fImpl->evtTree->Fill();

    auto fillClusters = [&](const std::vector<IonizationCluster>& clusters) {
        for (const auto& c : clusters) {
            fImpl->cb_eventID  = data.eventID;
            fImpl->cb_trackID  = c.trackID;
            fImpl->cb_parentID = c.parentID;
            fImpl->cb_x        = c.x;
            fImpl->cb_y        = c.y;
            fImpl->cb_z        = c.z;
            fImpl->cb_edep     = c.edep;
            fImpl->cb_nPrimary = c.nPrimary;
            fImpl->cb_ke       = c.kineticEnergy;
            fImpl->cb_time     = c.globalTime;
            setStr(fImpl->cb_volume,   32, c.volumeName);
            setStr(fImpl->cb_particle, 32, c.particleName);
            setStr(fImpl->cb_creatorProc, 32, c.creatorProcess);
            setStr(fImpl->cb_originProc,  32, c.originProcess);
            setStr(fImpl->cb_originVol,   32, c.originVolume);
            setStr(fImpl->cb_parentVol,   32, c.parentBirthVolume);
            fImpl->cb_originX    = c.originX;
            fImpl->cb_originY    = c.originY;
            fImpl->cb_originZ    = c.originZ;
            fImpl->cb_originE    = c.originEnergy;
            fImpl->cb_parentE    = c.parentEnergy;
            fImpl->cb_ancestorID = c.ancestorID;
            fImpl->clusTree->Fill();
        }
    };
    fillClusters(data.driftClusters);
    fillClusters(data.ampClusters);

    // Every volume that received energy this event.
    for (const auto& kv : data.edepByVolume) {
        fImpl->vb_eventID = data.eventID;
        setStr(fImpl->vb_volume, 64, kv.first);
        fImpl->vb_edep = kv.second;
        auto ns = data.nStepsByVolume.find(kv.first);
        fImpl->vb_nSteps = (ns != data.nStepsByVolume.end()) ? ns->second : 0;
        fImpl->volTree->Fill();
    }

#else
    fImpl->evtFile << data.eventID << ","
                   << data.edepDrift     << "," << data.edepAmp     << ","
                   << data.nPrimaryDrift << "," << data.nPrimaryAmp << ","
                   << data.driftClusters.size() << "," << data.ampClusters.size() << ","
                   << data.primaryInDrift << "," << data.primaryInAmp
                   << "," << static_cast<int>(data.interactionClass)
                   << "," << data.convProcess
                   << "," << data.convVolume
                   << "," << data.convParentBirthVolume
                   << "," << data.convZ
                   << "," << data.convEnergy
                   << "," << data.convParentEnergy
                   << "," << data.firstIntProcess
                   << "," << data.firstIntVolume
                   << "," << data.firstIntX
                   << "," << data.firstIntY
                   << "," << data.firstIntZ
                   << "," << data.firstIntEnergy
                   << "," << data.primaryPDG
                   << "," << data.primaryCharge
                   << "," << data.primaryKE_MeV
                   << "," << data.edepFrontGas
                   << "," << data.edepCathGas
                   << "," << data.edepBackGas
                   << "," << data.edepWindowGas
                   << "," << data.edepCathMylar
                   << "," << data.edepCathAl
                   << "," << data.edepWindowMylar
                   << "," << data.edepMeshP2
                   << "," << data.edepPadCu
                   << "," << data.edepFR4P2
                   << "," << data.edepCuB
                   << "," << data.edepFrame;
    if (isFull) {
        fImpl->evtFile << "," << data.edepHe3Gas
                       << "," << data.edepResistPaste
                       << "," << data.edepMylar
                       << "," << data.edepCathode
                       << "," << data.edepMicromesh
                       << "," << data.edepPCB
                       << "," << data.edepPCBKapton
                       << "," << data.edepPCBCu
                       << "," << data.edepPCBFR4
                       << "," << data.edepPCBRohacell
                       << "," << data.edepPCBAlFoil
                       << "," << data.edepScintWall
                       << "," << data.edepScintTape
                       << "," << data.edepScintAlFoil
                       << "," << data.edepLS1
                       << "," << data.edepLS2
                       << "," << data.primInHe3Gas
                       << "," << data.primInPCB
                       << "," << data.primInScintWall
                       << "," << data.primInLS1
                       << "," << data.primInLS2
                       << "," << data.edepLSCFRP
                       << "," << data.primInLSCFRP5;
    }
    if (isCalib) {
        fImpl->evtFile << "," << data.edepLS1
                       << "," << data.edepLSCFRP
                       << "," << data.primInLS1
                       << "," << data.edepBackScint
                       << "," << data.primInBackScint;
    }
    fImpl->evtFile << "\n";

    auto writeClusters = [&](const std::vector<IonizationCluster>& clusters) {
        for (const auto& c : clusters) {
            fImpl->clusFile << data.eventID << "," << c.trackID << "," << c.parentID << ","
                            << c.x << "," << c.y << "," << c.z << ","
                            << c.edep << "," << c.nPrimary << ","
                            << c.kineticEnergy << ","
                            << c.volumeName << "," << c.particleName << ","
                            << c.globalTime << ","
                            << c.creatorProcess << "," << c.originProcess << ","
                            << c.originVolume << "," << c.parentBirthVolume << ","
                            << c.originX << "," << c.originY << "," << c.originZ << ","
                            << c.originEnergy << "," << c.parentEnergy << ","
                            << c.ancestorID << "\n";
        }
    };
    writeClusters(data.driftClusters);
    writeClusters(data.ampClusters);

    for (const auto& kv : data.edepByVolume) {
        auto ns = data.nStepsByVolume.find(kv.first);
        fImpl->volFile << data.eventID << "," << kv.first << ","
                       << kv.second << ","
                       << (ns != data.nStepsByVolume.end() ? ns->second : 0)
                       << "\n";
    }
#endif
}
