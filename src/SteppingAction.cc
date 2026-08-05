// SteppingAction.cc
// Gas volumes: compute primary ion pairs via W-value, record per-cluster detail.
// Full-experiment extra volumes: record per-layer total edep and transmission flags.
// kLSCalib mode: score LS layer + back scint bar only.

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "EventData.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"

#include <cmath>

const std::map<std::string, double> SteppingAction::kWValues = {
    {"ArCF4",    34.0},
    {"HeEth",    27.0},
    {"ArCO2",    27.0},
    {"ArCF4Iso", 33.0},
    {"NeIso",    27.0},
    {"NeCF4",    30.0},
    {"ArCF4CO2", 34.0},
    {"ArIso",    26.0},
    {"PureCF4",  34.0},
    {"PureAr",   26.4},
    {"PureHe",   41.3},
    {"PureNe",   36.4},
    {"PureEthane",26.0},
    {"PureIso",  26.0},
    {"PureCO2",  33.0},
};

SteppingAction::SteppingAction(const SimConfig& cfg, EventAction* eventAction)
    : G4UserSteppingAction(), fConfig(cfg), fEventAction(eventAction) {}

double SteppingAction::GetWValue(const std::string& gas) const {
    auto it = kWValues.find(gas);
    return (it != kWValues.end()) ? it->second : 26.4;
}

// Named sums for the P2 stack, carved out of the generic edepByVolume map.
// Only the p2-mode volume names appear here; other modes fall through
// harmlessly and are still caught by the generic map.
void SteppingAction::ScoreP2Layer(EventData& data, const std::string& v,
                                  double edep_eV) const {
    auto starts = [&v](const char* p) {
        const std::string s(p);
        return v.size() >= s.size() && v.compare(0, s.size(), s) == 0;
    };

    // Gas that is outside the drift/amplification field: deposits here can
    // never produce a signal. Recording them is the whole point — it
    // separates "converted somewhere useless" from "did not convert".
    if      (v == "FrontGas")          data.edepFrontGas += edep_eV;
    else if (v == "DriftCathode_Gas")  data.edepCathGas  += edep_eV;
    else if (v == "BackGas")           data.edepBackGas  += edep_eV;
    else if (starts("FrontWindow_Gas") || starts("BackWindow_Gas"))
                                       data.edepWindowGas += edep_eV;
    // Solids.
    else if (starts("DriftCathode_Mylar")) data.edepCathMylar += edep_eV;
    else if (v == "DriftCathode_Al")       data.edepCathAl    += edep_eV;
    else if (starts("FrontWindow_Mylar") || starts("BackWindow_Mylar"))
                                           data.edepWindowMylar += edep_eV;
    else if (v == "Micromesh")             data.edepMeshP2    += edep_eV;
    else if (v == "PCB_Cu_F")              data.edepPadCu     += edep_eV;
    else if (v == "PCB_FR4")               data.edepFR4P2     += edep_eV;
    else if (v == "PCB_Cu_B")              data.edepCuB       += edep_eV;
    else if (starts("GasFrame"))           data.edepFrame     += edep_eV;
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    const G4Track* track = step->GetTrack();
    EventData& data = fEventAction->GetEventData();

    // ── First interaction of the primary ─────────────────────────────────
    // Must run BEFORE the edep<=0 return: a photon's photoabsorption or
    // Compton scatter deposits nothing on the photon's own step (the energy
    // leaves with the secondary), so an edep-gated hook would never see it.
    // This is what turns a photon run into a per-layer conversion budget.
    if (!data.firstIntRecorded && track->GetTrackID() == 1) {
        const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
        if (proc) {
            const std::string pn = proc->GetProcessName();
            if (pn != "Transportation" && pn != "CoupledTransportation" &&
                pn != "msc" && pn != "Rayl" && pn != "electronNuclear") {
                const G4VPhysicalVolume* ipv =
                    step->GetPostStepPoint()->GetPhysicalVolume();
                const G4ThreeVector ip = step->GetPostStepPoint()->GetPosition();
                data.firstIntRecorded = true;
                data.firstIntProcess  = pn;
                data.firstIntVolume   = ipv ? ipv->GetName() : "OutOfWorld";
                data.firstIntX = ip.x()/mm;
                data.firstIntY = ip.y()/mm;
                data.firstIntZ = ip.z()/mm;
                data.firstIntEnergy =
                    step->GetPreStepPoint()->GetKineticEnergy()/MeV;
            }
        }
    }

    G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.0) return;

    const G4VPhysicalVolume* pv = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!pv) return;

    const std::string volName = pv->GetName();
    int trackID  = track->GetTrackID();
    bool isPrimary = (trackID == 1);

    // ── Score EVERY volume ───────────────────────────────────────────────
    // Nothing is dropped on the floor: the generic map catches any volume,
    // named or not, so "converted somewhere that makes no signal" is
    // distinguishable from "did not convert".
    data.edepByVolume[volName]   += edep/eV;
    data.nStepsByVolume[volName] += 1;
    ScoreP2Layer(data, volName, edep/eV);

    // ── Gas volumes: ionisation cluster scoring (all non-vacuum modes) ───
    bool inDrift = (volName == "DriftGas");
    bool inAmp   = (volName == "AmpGas");

    if (inDrift || inAmp) {
        int parentID = track->GetParentID();
        const std::string pName = track->GetDefinition()->GetParticleName();
        double ke = step->GetPreStepPoint()->GetKineticEnergy();
        G4ThreeVector pos = 0.5*(step->GetPreStepPoint()->GetPosition() +
                                  step->GetPostStepPoint()->GetPosition());

        double W = GetWValue(fConfig.gas) * eV;
        int nPrimary = static_cast<int>(std::floor(edep / W));
        double frac  = edep / W - nPrimary;
        if (G4UniformRand() < frac) nPrimary++;

        IonizationCluster cluster;
        cluster.x            = pos.x()/mm;
        cluster.y            = pos.y()/mm;
        cluster.z            = pos.z()/mm;
        cluster.edep         = edep/eV;
        cluster.nPrimary     = nPrimary;
        cluster.trackID      = trackID;
        cluster.parentID     = parentID;
        cluster.volumeName   = volName;
        cluster.particleName = pName;
        cluster.kineticEnergy = ke/MeV;
        cluster.globalTime   = track->GetGlobalTime()/ns;

        const G4VProcess* cproc = track->GetCreatorProcess();
        cluster.creatorProcess = cproc ? cproc->GetProcessName() : "primary";

        // Provenance: which interaction, in which layer, started this chain.
        auto oit = data.trackOrigins.find(trackID);
        if (oit != data.trackOrigins.end()) {
            const TrackOrigin& o = oit->second;
            cluster.originProcess    = o.process;
            cluster.originVolume     = o.volume;
            cluster.originX          = o.x;
            cluster.originY          = o.y;
            cluster.originZ          = o.z;
            cluster.originEnergy     = o.energy;
            cluster.ancestorID       = o.ancestorID;
            cluster.parentEnergy     = o.parentEnergy;
            cluster.parentBirthVolume= o.parentBirthVolume;
        }

        if (inDrift) {
            data.edepDrift     += edep/eV;
            data.nPrimaryDrift += nPrimary;
            data.driftClusters.push_back(cluster);
            if (isPrimary) data.primaryInDrift = true;
        } else {
            data.edepAmp     += edep/eV;
            data.nPrimaryAmp += nPrimary;
            data.ampClusters.push_back(cluster);
            if (isPrimary) data.primaryInAmp = true;
        }
        return;
    }

    if (fConfig.mode == SimMode::kVacuum) return;

    // ── kLSCalib: score LS layer only ───────────────────────────────────
    if (fConfig.mode == SimMode::kLSCalib) {
        if (volName == "LiqScint_1") {
            data.edepLS1 += edep/eV;
            if (isPrimary) data.primInLS1 = true;
            return;
        }
        if ((volName.size() >= 8  && volName.substr(0,8)  == "LS_CFRP_")    ||
            (volName.size() >= 11 && volName.substr(0,11) == "LS_InnerCFR") ||
            (volName.size() >= 6  && volName.substr(0,6)  == "LS_Al_")) {
            data.edepLSCFRP += edep/eV;
            return;
        }
        return;
    }

    // ── kBackScintCalib: score back scint bar only ───────────────────────
    if (fConfig.mode == SimMode::kBackScintCalib) {
        if (volName == "BackScint") {
            data.edepBackScint += edep/eV;
            if (isPrimary) data.primInBackScint = true;
            return;
        }
        return;
    }

    // ── Full / Sr90 modes: per-layer edep ───────────────────────────────
    if (volName == "ResistivePaste") {
        data.edepResistPaste += edep/eV;
        return;
    }
    if (volName == "He3Gas") {
        data.edepHe3Gas += edep/eV;
        if (isPrimary) data.primInHe3Gas = true;
        return;
    }
    if (volName == "GasWindow_Mylar") {
        data.edepMylar += edep/eV;
        return;
    }
    if (volName == "GasWindow_Al" ||
        volName == "DriftCathode_Kapton" ||
        volName == "DriftCathode_Cu") {
        data.edepCathode += edep/eV;
        return;
    }
    if (volName == "Micromesh") {
        data.edepMicromesh += edep/eV;
        return;
    }
    if (volName.size() >= 4 && volName.substr(0,4) == "PCB_") {
        data.edepPCB += edep/eV;
        if (isPrimary) data.primInPCB = true;
        if      (volName == "PCB_Kapton")   data.edepPCBKapton   += edep/eV;
        else if (volName.size()>=6 && volName.substr(0,6)=="PCB_Cu")  data.edepPCBCu       += edep/eV;
        else if (volName.size()>=7 && volName.substr(0,7)=="PCB_FR4") data.edepPCBFR4      += edep/eV;
        else if (volName == "PCB_Rohacell") data.edepPCBRohacell += edep/eV;
        else if (volName == "PCB_AlFoil")   data.edepPCBAlFoil   += edep/eV;
        return;
    }
    if (volName == "ScintWall_BlackTape1" || volName == "ScintWall_BlackTape2") {
        data.edepScintTape += edep/eV;
        return;
    }
    if (volName == "ScintWall_AlFoil") {
        data.edepScintAlFoil += edep/eV;
        return;
    }
    if (volName == "PlasticScint") {
        data.edepScintWall += edep/eV;
        if (isPrimary) data.primInScintWall = true;
        return;
    }
    if (volName == "LiqScint_1") {
        data.edepLS1 += edep/eV;
        if (isPrimary) data.primInLS1 = true;
        return;
    }
    if (volName == "LiqScint_2") {
        data.edepLS2 += edep/eV;
        if (isPrimary) data.primInLS2 = true;
        return;
    }
    // LS structural CFRP walls + inner CFRP liners + Al liners → all into edepLSCFRP
    if (volName.size() >= 8 && volName.substr(0,8) == "LS_CFRP_") {
        data.edepLSCFRP += edep/eV;
        // LS_CFRP_3 is the back wall (exited LS2) — reusing primInLSCFRP5 flag
        if (isPrimary && volName == "LS_CFRP_3") data.primInLSCFRP5 = true;
        return;
    }
    if (volName.size() >= 11 && volName.substr(0,11) == "LS_InnerCFR") {
        data.edepLSCFRP += edep/eV;
        return;
    }
    if (volName.size() >= 6 && volName.substr(0,6) == "LS_Al_") {
        data.edepLSCFRP += edep/eV;
        return;
    }
}
