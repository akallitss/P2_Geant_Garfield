#pragma once
// SteppingAction.hh
// Intercepts every step in the gas volumes to:
//   1. Record energy deposition
//   2. Estimate number of primary ion pairs via W-value
//   3. Record cluster position

#include "G4UserSteppingAction.hh"
#include "SimConfig.hh"
#include <map>
#include <string>

class EventAction;
class G4Step;
struct EventData;

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(const SimConfig& cfg, EventAction* eventAction);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

private:
    // W-value (mean energy per ion pair) [eV] for the configured mixture,
    // resolved from gas::MixtureWValue at construction. See GasMixtures.cc
    // for the weighting and for the Penning-transfer caveat that makes every
    // neon-mixture value an upper bound.
    double GetWValue(const std::string& gas) const;

    // Named P2-stack sums, carved out of the generic per-volume map.
    void ScoreP2Layer(EventData& data, const std::string& volName,
                      double edep_eV) const;

    const SimConfig& fConfig;
    EventAction*     fEventAction;

    double fWValue = 26.4;   // eV, resolved once in the constructor
};
