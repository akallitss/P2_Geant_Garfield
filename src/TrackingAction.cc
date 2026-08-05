// TrackingAction.cc
//
// Geant4 tracks depth-first, so a parent is always tracked before its
// secondaries. That lets us keep the provenance chain in a plain
// map<trackID, TrackOrigin> (owned by EventData, cleared each event) and
// look the parent up here — no G4VUserTrackInformation and no const_cast on
// the secondary list.

#include "TrackingAction.hh"
#include "EventAction.hh"
#include "EventData.hh"

#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

TrackingAction::TrackingAction(EventAction* eventAction)
    : G4UserTrackingAction(), fEventAction(eventAction) {}

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
    EventData& data = fEventAction->GetEventData();

    const int    tid    = track->GetTrackID();
    const int    pid    = track->GetParentID();
    const double charge = track->GetDefinition()->GetPDGCharge();
    const int    pdg    = track->GetDefinition()->GetPDGEncoding();
    const auto   pos    = track->GetPosition();
    const double ke     = track->GetKineticEnergy();

    const G4VPhysicalVolume* pv = track->GetVolume();
    const std::string vol = pv ? pv->GetName() : "OutOfWorld";

    TrackOrigin org;
    org.ownCharge = charge;

    if (pid == 0) {
        // Primary: it is its own origin.
        org.process    = "primary";
        org.volume     = vol;
        org.x = pos.x()/mm; org.y = pos.y()/mm; org.z = pos.z()/mm;
        org.energy     = ke/MeV;
        org.pdg        = pdg;
        org.ancestorID = tid;
        org.parentEnergy      = ke/MeV;
        org.parentBirthVolume = vol;
    } else {
        auto it = data.trackOrigins.find(pid);
        const bool haveParent   = (it != data.trackOrigins.end());
        const bool parentNeutral = haveParent ? (it->second.ownCharge == 0.0)
                                              : true;

        if (haveParent && !parentNeutral) {
            // Parent is charged: this is a delta ray, a brems photon, or any
            // other descendant of an already-converted particle. Keep the
            // attribution of whatever made the parent.
            org = it->second;
            org.ownCharge = charge;          // but record OUR charge
        } else {
            // Parent is neutral (or unknown): this track is the point at
            // which a neutral became charged — i.e. the conversion we bin by.
            const G4VProcess* cp = track->GetCreatorProcess();
            org.process    = cp ? cp->GetProcessName() : "unknown";
            org.volume     = vol;
            org.x = pos.x()/mm; org.y = pos.y()/mm; org.z = pos.z()/mm;
            org.energy     = ke/MeV;
            org.pdg        = pdg;
            org.ancestorID = tid;
            // The neutral that converted: its energy AT ITS OWN CREATION and
            // where it was born. For the Cu-K fluorescence chain this reads
            // 0.00805 MeV / PCB_Cu_F, which is how the analysis separates it
            // from a direct conversion of the incident X-ray.
            if (haveParent) {
                org.parentEnergy      = it->second.energy;
                org.parentBirthVolume = it->second.volume;
            }
        }
    }

    data.trackOrigins[tid] = org;
}
