#pragma once
// EventData.hh — per-event data container

#include "G4ThreeVector.hh"
#include "TrackOrigin.hh"

#include <map>
#include <string>
#include <vector>

// ── Interaction classes ──────────────────────────────────────────────────
// Coarse per-event label, computed in EventAction::EndOfEventAction from the
// provenance of whichever ancestor deposited the most energy in the active
// gas. The raw origin fields are written out too, so analysis can always
// re-cut more finely than this; the enum exists so the headline plot
// ("event-by-event edep separated by interaction type") is one TTree::Draw.
// See docs/research/PHOTON_DISCRIMINATION_NOTES.md §6.2.
enum class InteractionClass : int {
    kNone            = 0,  // nothing deposited in the active gas
    kPrimaryTrack    = 1,  // a charged primary crossed the gas — the MIP signal
    kGasPhotoelectric= 2,  // incident photon photoabsorbed in DriftGas/AmpGas
    kGasCompton      = 3,  // Compton scatter in the active gas
    kGasOther        = 4,  // pair production / other conversion in the gas
    kWallPhotoelectric=5,  // photoabsorption in a solid, electron reached gas
    kWallCompton     = 6,  // Compton in a solid, electron reached gas
    kWallOther       = 7,
    kFluorescence    = 8,  // soft fluorescence photon born in a solid,
                           //   re-absorbed in the gas (Cu-K 8.05, Fe-K 6.4 keV)
    kOther           = 9,
};

const char* InteractionClassName(InteractionClass c);

// True for the two volumes that can actually produce a signal.
inline bool IsActiveGas(const std::string& v) {
    return v == "DriftGas" || v == "AmpGas";
}

struct IonizationCluster {
    double x, y, z;
    double edep;
    int    nPrimary;
    int    trackID;
    int    parentID;
    std::string volumeName;
    std::string particleName;
    double kineticEnergy;

    // ── Provenance (added 2026-08-05) ─────────────────────────────────────
    std::string creatorProcess;      // process that made THIS track
    std::string originProcess;       // process that started the chain
    std::string originVolume;        // where that happened
    double      originX = 0, originY = 0, originZ = 0;
    double      originEnergy = 0.0;  // [MeV]
    int         ancestorID   = 0;
    double      parentEnergy = 0.0;  // converting neutral's creation KE [MeV]
    std::string parentBirthVolume;
    double      globalTime = 0.0;    // [ns]
};

struct EventData {
    int eventID = -1;

    // ── Micromegas gas scoring ─────────────────────────────────────────────
    double edepDrift = 0.0;
    double edepAmp   = 0.0;

    std::vector<IonizationCluster> driftClusters;
    std::vector<IonizationCluster> ampClusters;

    int  nPrimaryDrift = 0;
    int  nPrimaryAmp   = 0;
    bool primaryInDrift = false;
    bool primaryInAmp   = false;

    // ── Score EVERYTHING (added 2026-08-05) ───────────────────────────────
    // Every volume that receives energy, by name — nothing is silently
    // dropped. Written to the VolumeTree, one row per (event, volume).
    // The named P2 sums below are conveniences carved out of this map.
    std::map<std::string, double> edepByVolume;   // [eV]
    std::map<std::string, int>    nStepsByVolume;

    // Named P2-stack sums [eV] — the layers the conversion budget needs.
    double edepFrontGas   = 0.0;   // upstream of the drift cathode: no signal
    double edepCathGas    = 0.0;   // between the two cathode foils: no signal
    double edepBackGas    = 0.0;   // behind the PCB: no signal
    double edepWindowGas  = 0.0;   // inside the bulged windows: no signal
    double edepCathMylar  = 0.0;
    double edepCathAl     = 0.0;
    double edepWindowMylar= 0.0;
    double edepMeshP2     = 0.0;
    double edepPadCu      = 0.0;   // PCB_Cu_F — faces the amplification gap
    double edepPadGapGas  = 0.0;   // gas in the etched recesses of both Cu
                                   // layers (inter-pad grooves, fan-out): no signal
    double edepFR4P2      = 0.0;
    double edepCuB        = 0.0;
    double edepFrame      = 0.0;

    // ── Provenance bookkeeping ────────────────────────────────────────────
    std::map<int, TrackOrigin> trackOrigins;   // trackID -> origin, per event

    // First interaction of the primary (photon runs: where the X-ray went).
    bool        firstIntRecorded = false;
    std::string firstIntProcess;
    std::string firstIntVolume;
    double      firstIntX = 0, firstIntY = 0, firstIntZ = 0;
    double      firstIntEnergy = 0.0;   // primary KE before the interaction [MeV]

    // Per-event classification.
    InteractionClass interactionClass = InteractionClass::kNone;
    std::string      convProcess;       // origin process of the dominant chain
    std::string      convVolume;        // origin volume of it
    double           convZ = 0.0;       // [mm]
    double           convEnergy = 0.0;  // [MeV]
    double           convParentEnergy = 0.0;
    std::string      convParentBirthVolume;

    int    primaryPDG    = 0;
    double primaryCharge = 0.0;

    // ── Full-experiment per-layer edep [eV] ───────────────────────────────
    double edepHe3Gas      = 0.0;
    double edepResistPaste = 0.0;

    // ── MM entrance / dead layers ─────────────────────────────────────────
    double edepMylar     = 0.0;
    double edepCathode   = 0.0;
    double edepMicromesh = 0.0;

    // ── PCB stack ─────────────────────────────────────────────────────────
    double edepPCB         = 0.0;
    double edepPCBKapton   = 0.0;
    double edepPCBCu       = 0.0;
    double edepPCBFR4      = 0.0;
    double edepPCBRohacell = 0.0;
    double edepPCBAlFoil   = 0.0;

    // ── Scintillator wall ─────────────────────────────────────────────────
    double edepScintWall   = 0.0;
    double edepScintTape   = 0.0;
    double edepScintAlFoil = 0.0;

    // ── Liquid scintillator stack (2 layers × 2 cm) ───────────────────────
    double edepLS1    = 0.0;
    double edepLS2    = 0.0;
    double edepLSCFRP = 0.0;  // all structural + inner CFRP/Al walls combined

    // ── Back plastic scintillator (kBackScintCalib) ────────────────────────
    double edepBackScint = 0.0;  // PVT active volume [eV]

    // ── Primary beam kinetic energy (kLSCalib / kBackScintCalib) ──────────
    // Stored so the analysis can see what spectrum was sampled.
    double primaryKE_MeV = 0.0;

    // ── Transmission flags ─────────────────────────────────────────────────
    bool primInHe3Gas    = false;
    bool primInPCB       = false;
    bool primInScintWall = false;
    bool primInLS1       = false;
    bool primInLS2       = false;
    bool primInLSCFRP5   = false;  // primary exited LS stack (reached back CFRP wall)
    bool primInBackScint = false;

    void Reset() {
        eventID = -1;
        edepDrift = edepAmp = 0.0;
        nPrimaryDrift = nPrimaryAmp = 0;
        primaryInDrift = primaryInAmp = false;
        driftClusters.clear();
        ampClusters.clear();

        edepByVolume.clear();
        nStepsByVolume.clear();
        edepFrontGas = edepCathGas = edepBackGas = edepWindowGas = 0.0;
        edepCathMylar = edepCathAl = edepWindowMylar = 0.0;
        edepMeshP2 = edepPadCu = edepFR4P2 = edepCuB = edepFrame = 0.0;
        edepPadGapGas = 0.0;

        trackOrigins.clear();
        firstIntRecorded = false;
        firstIntProcess.clear();
        firstIntVolume.clear();
        firstIntX = firstIntY = firstIntZ = 0.0;
        firstIntEnergy = 0.0;

        interactionClass = InteractionClass::kNone;
        convProcess.clear();
        convVolume.clear();
        convParentBirthVolume.clear();
        convZ = convEnergy = convParentEnergy = 0.0;
        primaryPDG = 0;
        primaryCharge = 0.0;

        edepHe3Gas = edepResistPaste = 0.0;
        edepMylar = edepCathode = edepMicromesh = 0.0;
        edepPCB = edepPCBKapton = edepPCBCu = edepPCBFR4
                = edepPCBRohacell = edepPCBAlFoil = 0.0;
        edepScintWall = edepScintTape = edepScintAlFoil = 0.0;
        edepLS1 = edepLS2 = edepLSCFRP = 0.0;
        edepBackScint = 0.0;
        primaryKE_MeV = 0.0;
        primInHe3Gas = primInPCB = primInScintWall = false;
        primInLS1 = primInLS2 = primInLSCFRP5 = primInBackScint = false;
    }
};
