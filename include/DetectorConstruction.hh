#pragma once
// DetectorConstruction.hh

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "SimConfig.hh"
#include <map>
#include <string>

class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(const SimConfig& cfg);
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    G4LogicalVolume* GetDriftGasLV()    const { return fDriftGasLV; }
    G4LogicalVolume* GetAmpGasLV()      const { return fAmpGasLV; }
    G4double         GetHe3GasCenterZ() const { return fHe3GasCenterZ; }
    // Upstream-most z of the P2 stack (front window bulge apex), kP2Wedge mode.
    G4double         GetP2FrontZ()      const { return fP2FrontZ; }
    // Drift-gap mid-plane z — the aim point the angled gun pivots about.
    G4double         GetP2DriftCenterZ() const { return fP2DriftCenterZ; }

private:
    void DefineMaterials();
    G4Material* GetGasMixture(const std::string& name);
    G4VPhysicalVolume* ConstructP2();

    const SimConfig& fConfig;
    G4LogicalVolume* fDriftGasLV   = nullptr;
    G4LogicalVolume* fAmpGasLV     = nullptr;
    G4LogicalVolume* fHe3GasLV     = nullptr;
    G4LogicalVolume* fBackScintLV  = nullptr;  // back plastic scintillator bar (kLSCalib)

    G4double fHe3GasCenterZ  = 0.0;
    G4double fP2FrontZ       = 0.0;
    G4double fP2DriftCenterZ = 0.0;

    std::map<std::string, G4Material*> fGasMaterials;
};
