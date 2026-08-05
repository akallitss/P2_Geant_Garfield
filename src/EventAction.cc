// EventAction.cc

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <map>

const char* InteractionClassName(InteractionClass c) {
    switch (c) {
        case InteractionClass::kNone:             return "none";
        case InteractionClass::kPrimaryTrack:     return "primary-track";
        case InteractionClass::kGasPhotoelectric: return "gas-photoelectric";
        case InteractionClass::kGasCompton:       return "gas-compton";
        case InteractionClass::kGasOther:         return "gas-other";
        case InteractionClass::kWallPhotoelectric:return "wall-photoelectric";
        case InteractionClass::kWallCompton:      return "wall-compton";
        case InteractionClass::kWallOther:        return "wall-other";
        case InteractionClass::kFluorescence:     return "fluorescence-reabs";
        default:                                  return "other";
    }
}

EventAction::EventAction(const SimConfig& cfg, RunAction* runAction)
    : G4UserEventAction(), fConfig(cfg), fRunAction(runAction) {}

// Label the event by the provenance of whichever ancestor chain deposited
// the most energy in the active gas. Raw origin fields are written out
// alongside, so analysis can always re-cut more finely than this label.
void EventAction::ClassifyEvent() {
    std::map<int, double> edepByAncestor;
    const IonizationCluster* dominant = nullptr;

    auto scan = [&](const std::vector<IonizationCluster>& cl) {
        for (const auto& c : cl) edepByAncestor[c.ancestorID] += c.edep;
    };
    scan(fData.driftClusters);
    scan(fData.ampClusters);

    if (edepByAncestor.empty()) {
        fData.interactionClass = InteractionClass::kNone;
        return;
    }

    int bestID = edepByAncestor.begin()->first;
    double bestE = -1.0;
    for (const auto& kv : edepByAncestor)
        if (kv.second > bestE) { bestE = kv.second; bestID = kv.first; }

    auto pick = [&](const std::vector<IonizationCluster>& cl) {
        for (const auto& c : cl)
            if (c.ancestorID == bestID && !dominant) dominant = &c;
    };
    pick(fData.driftClusters);
    pick(fData.ampClusters);
    if (!dominant) { fData.interactionClass = InteractionClass::kOther; return; }

    fData.convProcess           = dominant->originProcess;
    fData.convVolume            = dominant->originVolume;
    fData.convZ                 = dominant->originZ;
    fData.convEnergy            = dominant->originEnergy;
    fData.convParentEnergy      = dominant->parentEnergy;
    fData.convParentBirthVolume = dominant->parentBirthVolume;

    const std::string& proc = dominant->originProcess;
    const bool inGas = IsActiveGas(dominant->originVolume);

    if (proc == "primary") {
        fData.interactionClass = InteractionClass::kPrimaryTrack;
    } else if (inGas) {
        // A soft photon born in a solid and re-absorbed in the gas is the
        // fluorescence chain (Cu-K 8.05 keV, Fe-K 6.4 keV), not a direct
        // conversion of the incident X-ray. Separating it matters: it is
        // the one gas-sensitive wall channel and the easiest fake to reject.
        const bool softFromSolid =
            proc == "phot" &&
            !dominant->parentBirthVolume.empty() &&
            !IsActiveGas(dominant->parentBirthVolume) &&
            dominant->parentEnergy > 0.0 &&
            dominant->parentEnergy < 0.020;      // < 20 keV [MeV units]
        if      (softFromSolid)    fData.interactionClass = InteractionClass::kFluorescence;
        else if (proc == "phot")   fData.interactionClass = InteractionClass::kGasPhotoelectric;
        else if (proc == "compt")  fData.interactionClass = InteractionClass::kGasCompton;
        else                       fData.interactionClass = InteractionClass::kGasOther;
    } else {
        if      (proc == "phot")   fData.interactionClass = InteractionClass::kWallPhotoelectric;
        else if (proc == "compt")  fData.interactionClass = InteractionClass::kWallCompton;
        else                       fData.interactionClass = InteractionClass::kWallOther;
    }
}

void EventAction::BeginOfEventAction(const G4Event* event) {
    fData.Reset();
    fData.eventID = event->GetEventID();

    // Record the sampled primary KE (useful when spectrum sampling is on)
    if (event->GetNumberOfPrimaryVertex() > 0) {
        const G4PrimaryVertex* vtx = event->GetPrimaryVertex(0);
        if (vtx && vtx->GetNumberOfParticle() > 0) {
            const G4PrimaryParticle* p = vtx->GetPrimary(0);
            if (p) {
                fData.primaryKE_MeV = p->GetKineticEnergy() / MeV;
                fData.primaryPDG    = p->GetPDGcode();
                if (p->GetParticleDefinition())
                    fData.primaryCharge = p->GetParticleDefinition()->GetPDGCharge();
            }
        }
    }
}

void EventAction::EndOfEventAction(const G4Event* event) {
    ClassifyEvent();

    // Pass event summary to RunAction for accumulation / writing
    fRunAction->RecordEvent(fData);

    if (fConfig.verbose && (fData.eventID % 1000 == 0)) {
        G4cout << "[Event " << fData.eventID << "]"
               << "  Edep_drift=" << fData.edepDrift / eV << " eV"
               << "  N_prim_drift=" << fData.nPrimaryDrift
               << "  Edep_amp=" << fData.edepAmp / eV << " eV"
               << "  N_prim_amp=" << fData.nPrimaryAmp
               << "  class=" << InteractionClassName(fData.interactionClass)
               << (fData.convVolume.empty() ? "" : " @ " + fData.convVolume)
               << G4endl;
    }
}
