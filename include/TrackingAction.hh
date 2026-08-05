#pragma once
// TrackingAction.hh
// Builds the per-track provenance chain (TrackOrigin) so that every scored
// ionization can be traced back to the interaction that started it.
// See TrackOrigin.hh for the inheritance rule and
// docs/research/PHOTON_DISCRIMINATION_NOTES.md §6.2 for why it is needed.

#include "G4UserTrackingAction.hh"

class EventAction;
class G4Track;

class TrackingAction : public G4UserTrackingAction {
public:
    explicit TrackingAction(EventAction* eventAction);
    ~TrackingAction() override = default;

    void PreUserTrackingAction(const G4Track* track) override;

private:
    EventAction* fEventAction;
};
