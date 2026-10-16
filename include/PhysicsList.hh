#pragma once
// PhysicsList.hh
// Physics list for Micromegas simulation.
// Uses FTFP_BERT as hadronic base + G4EmStandardPhysics_option4 for EM
// (option4 = best EM, required for accurate low-energy electron/photon ionization)
// + G4RadioactiveDecay + G4StepLimiterPhysics

#include "G4VModularPhysicsList.hh"
#include "SimConfig.hh"

class PhysicsList : public G4VModularPhysicsList {
public:
    explicit PhysicsList(const SimConfig& cfg);
    ~PhysicsList() override = default;
    void ConstructProcess() override;
    void SetCuts() override;

private:
    std::string fPAI;   // SimConfig::pai
    bool        fPixe;  // SimConfig::pixe
};
