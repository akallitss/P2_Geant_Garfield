#pragma once
// PrimaryGeneratorAction.hh

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "SimConfig.hh"
#include <memory>
#include <vector>
#include <string>

class G4Event;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(const SimConfig& cfg,
                           const DetectorConstruction* detCon = nullptr);
    ~PrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event* event) override;

private:
    void LoadSpectrum(const std::string& filepath);
    double SampleSpectrum() const;

    std::unique_ptr<G4ParticleGun> fGun;
    const SimConfig&               fConfig;
    const DetectorConstruction*    fDetCon;

    // P2 beam geometry, resolved once in the constructor. Kept as members
    // because --beam-spread re-randomises the impact point every event, so
    // the gun position cannot be set once and forgotten.
    G4ThreeVector fP2Aim{0, 0, 0};    // aim point at the drift mid-plane
    G4ThreeVector fP2Dir{0, 0, 1};    // beam direction
    G4double      fP2Standoff = 0.0;  // aim point -> gun, along the beam
    bool          fP2Mode     = false;

    // Sr-90/Y-90 spectrum: CDF table for inverse-transform sampling
    std::vector<double> fSpecEnergies;  // energy values [MeV]
    std::vector<double> fSpecCDF;       // cumulative probabilities [0,1]
    bool fUseSpectrum = false;
};
