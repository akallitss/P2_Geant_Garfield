#pragma once
// HeedGasModel.hh — primary ionization in the drift gas from Garfield++/Heed.
//
// Why: Geant4 turns a gas step's energy into electrons as edep / W. That has
// the right mean but no notion of a collision, so a 20 eV collision can give
// zero electrons, and events with few collisions come out with too few
// electrons. On a bare 4 mm Ar/iso slab this doubles the fraction of 100 MeV
// e- with < 10 primaries (Geant4 4.1 %, Heed 1.8 %; docs/results/heed/pai).
// Heed simulates each ionizing collision (>= 1 electron each) and transports
// its delta electrons, as validated against the Heed standalone runs.
//
// How: a G4VFastSimulationModel on the region P2HeedGas (DriftGas only). A
// charged particle above a per-species threshold that is in the drift gas is
// handed to TrackHeed along a straight line to the exit of the volume. Every
// Heed electron becomes one IonizationCluster (nPrimary = 1) at its own
// position and time. The particle then continues in Geant4 from the exit
// point, same direction, kinetic energy reduced by the Heed energy loss.
// Below threshold (soft electrons: photoelectrons, deltas from the solids)
// Geant4 keeps tracking and SteppingAction converts edep / W as before.
//
// Limits: straight line (no multiple scattering inside the gas; 3 deg for a
// 1 MeV e- over 4 mm, hence the 1 MeV e+- threshold); Heed's own W and
// Fano; Heed delta electrons that would leave the gas stop at its boundary.
// Not thread-safe (Heed): mm_sim refuses --heed with more than one thread.

#include <memory>
#include <string>

#include "G4VFastSimulationModel.hh"
#include "SimConfig.hh"

class G4VSolid;
namespace Garfield {
class MediumMagboltz;
class ComponentConstant;
class Sensor;
class TrackHeed;
}

class HeedGasModel : public G4VFastSimulationModel {
public:
    HeedGasModel(const G4String& name, G4Region* envelope, const SimConfig& cfg,
                 const G4VSolid* gasSolid);
    ~HeedGasModel() override;

    G4bool IsApplicable(const G4ParticleDefinition& p) override;
    G4bool ModelTrigger(const G4FastTrack& track) override;
    void   DoIt(const G4FastTrack& track, G4FastStep& step) override;

    // Heed's name for a Geant4 particle and its minimum kinetic energy [MeV];
    // empty name = not handled.
    static std::string HeedName(const std::string& g4name);
    double MinEnergyMeV(const std::string& g4name) const;

private:
    const SimConfig&  fConfig;
    const G4VSolid*   fSolid;
    std::unique_ptr<Garfield::MediumMagboltz>    fGas;
    std::unique_ptr<Garfield::ComponentConstant> fCmp;
    std::unique_ptr<Garfield::Sensor>            fSensor;
    std::unique_ptr<Garfield::TrackHeed>         fHeed;
    static constexpr double kHalfXY = 100.0;   // mm, Heed box half-width around the entry
    double      fZlo = 0, fZhi = 0;            // local z range of the gas volume
    std::string fHeedParticle;   // what TrackHeed is set up for now
    double      fHeedEnergy = -1;
    int         fFailures = 0;                 // NewTrack failures reported
};
