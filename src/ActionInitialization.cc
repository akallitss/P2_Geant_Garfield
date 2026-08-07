// ActionInitialization.cc

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

ActionInitialization::ActionInitialization(const SimConfig& cfg,
                                            const DetectorConstruction* detCon)
    : G4VUserActionInitialization(), fConfig(cfg), fDetCon(detCon) {}

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction(fConfig, true, fDetCon));
}

void ActionInitialization::Build() const {
    SetUserAction(new PrimaryGeneratorAction(fConfig, fDetCon));

    auto* runAction   = new RunAction(fConfig, false, fDetCon);
    auto* eventAction = new EventAction(fConfig, runAction);

    SetUserAction(runAction);
    SetUserAction(eventAction);
    // TrackingAction must be registered: it builds the provenance chain that
    // SteppingAction reads when it writes clusters.
    SetUserAction(new TrackingAction(eventAction));
    SetUserAction(new SteppingAction(fConfig, eventAction));
}
