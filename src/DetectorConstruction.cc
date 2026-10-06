// DetectorConstruction.cc
//
// Modes:
//  kVacuum          : Micromegas in vacuum, optional Al shielding upstream.
//  kFullExperiment  : He-3 target → air → MM → PCB → air → scint wall → air → LS stack.
//  kSr90Calibration : Sr-90 source in air → MM → PCB → air → scint wall → air → LS stack.
//  kSr90NoMM        : Sr-90 source in air → scint wall → air → LS stack (no MM/PCB).
//  kLSCalib         : Sr-90 source capsule → air → 1 LS layer → air → back scint bar.
//
// Geometry updated to match Full_Geant (4-arm X17 sim):
//  - He-3 target: r=1.5 cm, L=8 cm (was r=2.5 cm, L=15 cm)
//  - LS: 2 layers × 2 cm LAB, each preceded by inner CFRP liner + Al liner
//  - Scint wall: BlackMylar tape (200 µm) instead of PVC tape (165 µm)

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4PVReplica.hh"

#include "SensitiveDetector.hh"
#include "P2Wedge.hh"
#include "P2PadMap.hh"
#include "P2Pillars.hh"
#include "GasMixtures.hh"

#include <stdexcept>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>

// ============================================================
DetectorConstruction::DetectorConstruction(const SimConfig& cfg)
    : G4VUserDetectorConstruction(), fConfig(cfg) {}

// ============================================================
void DetectorConstruction::DefineMaterials() {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* elH  = nist->FindOrBuildElement("H");
    G4Element* elC  = nist->FindOrBuildElement("C");
    G4Element* elN  = nist->FindOrBuildElement("N");
    G4Element* elO  = nist->FindOrBuildElement("O");
    G4Element* elSi = nist->FindOrBuildElement("Si");
    G4Element* elAr = nist->FindOrBuildElement("Ar");
    G4Element* elNe = nist->FindOrBuildElement("Ne");
    G4Element* elHe = nist->FindOrBuildElement("He");
    G4Element* elF  = nist->FindOrBuildElement("F");

    // ── Pure gases ───────────────────────────────────────────
    //
    // Densities are computed from the molar mass by the ideal-gas law at the
    // *declared* conditions (20 C, 1 atm) rather than hard-coded.
    //
    // Fixed 2026-08-05: every density here used to be the 0 C STP value
    // (Ar 1.782, Ne 0.8999, CO2 1.977 mg/cm3 ...) while the G4Material was
    // declared at 293.15 K -- so every gas was 7.3 % too dense, which feeds
    // straight into dE/dx, primary-ionization counts and photon interaction
    // probability. See docs/research/GAS_FIX_NOTES.md.
    //
    // Ideal gas is good to <1 % for the noble gases and CO2/CF4 at 1 atm.
    // The worst real-gas deviation is isobutane (~2 % denser than ideal at
    // 20 C); at <=20 % quencher that is a <0.4 % error on any mixture here,
    // well below the P/T control of a real gas system.
    const G4double kGasT = 293.15*kelvin, kGasP = 1*atmosphere;
    auto idealRho = [](G4double molarMass_g_mol) {
        const G4double R = 8.314462618;      // J/(mol K)
        const G4double T = 293.15, P = 101325.0;
        return P * (molarMass_g_mol * 1e-3) / (R * T) * 1e-3;   // g/cm3
    };

    G4Material* isobutane = new G4Material("Isobutane", idealRho(58.1222)*g/cm3, 2,
                                            kStateGas, kGasT, kGasP);
    isobutane->AddElement(elC, 4);
    isobutane->AddElement(elH, 10);

    G4Material* ethane = new G4Material("Ethane", idealRho(30.069)*g/cm3, 2,
                                         kStateGas, kGasT, kGasP);
    ethane->AddElement(elC, 2);
    ethane->AddElement(elH, 6);

    G4Material* methane = new G4Material("Methane", idealRho(16.0425)*g/cm3, 2,
                                          kStateGas, kGasT, kGasP);
    methane->AddElement(elC, 1);
    methane->AddElement(elH, 4);

    G4Material* nitrogen = new G4Material("N2_gas", idealRho(28.0134)*g/cm3, 1,
                                           kStateGas, kGasT, kGasP);
    nitrogen->AddElement(elN, 2);

    G4Material* CO2 = new G4Material("CO2_gas", idealRho(44.0095)*g/cm3, 2,
                                      kStateGas, kGasT, kGasP);
    CO2->AddElement(elC, 1);
    CO2->AddElement(elO, 2);

    G4Material* CF4 = new G4Material("CF4_gas", idealRho(88.0043)*g/cm3, 2,
                                      kStateGas, kGasT, kGasP);
    CF4->AddElement(elC, 1);
    CF4->AddElement(elF, 4);

    G4Material* purAr = new G4Material("PureArgon", idealRho(39.948)*g/cm3, 1,
                                        kStateGas, kGasT, kGasP);
    purAr->AddElement(elAr, 1);

    G4Material* pureHe = new G4Material("PureHe", idealRho(4.002602)*g/cm3, 1,
                                         kStateGas, kGasT, kGasP);
    pureHe->AddElement(elHe, 1);

    G4Material* pureNe = new G4Material("PureNe", idealRho(20.1797)*g/cm3, 1,
                                         kStateGas, kGasT, kGasP);
    pureNe->AddElement(elNe, 1);

    // ── Gas mixtures ─────────────────────────────────────────
    //
    // Composition now comes from gas::All() (include/GasMixtures.hh), which
    // is also what SteppingAction reads its W-value from and what the run
    // metadata records. Nothing about a mixture is written down twice.
    //
    // Mixtures are specified by VOLUME fraction, which is how gas systems are
    // actually mixed and how every mixture in the literature is quoted.
    //
    // Fixed 2026-08-05: these fractions used to be handed straight to
    // G4Material::AddMaterial, which takes **mass** fractions. For Ar/Iso the
    // error is modest (95/5 by volume = 92.9/7.1 by mass) but for the neon
    // mixtures it is severe -- Ne/Iso 95/5 by volume is **86.8/13.2 by
    // mass**, so the simulated gas had roughly 2.6x the intended isobutane
    // content by mass. The mixture density was also computed separately at
    // each call site, giving two places to get it wrong.
    //
    // Now both follow from the volume fractions: for ideal gases at a common
    // T and P the partial densities add, so
    //     rho_mix     = sum_i f_vol_i * rho_i
    //     f_mass_i    = f_vol_i * rho_i / rho_mix
    std::map<std::string, G4Material*> pure = {
        {"Ar",     purAr},
        {"Ne",     pureNe},
        {"He",     pureHe},
        {"CO2",    CO2},
        {"CF4",    CF4},
        {"CH4",    methane},
        {"C2H6",   ethane},
        {"iC4H10", isobutane},
        {"N2",     nitrogen},
    };

    for (const auto& mix : gas::All()) {
        // A single-component "mixture" is the pure gas itself; G4Material
        // refuses a one-component AddMaterial fraction of 1 in some builds,
        // and there is nothing to mix anyway.
        if (mix.parts.size() == 1) {
            fGasMaterials[mix.name] = pure.at(mix.parts[0].component);
            continue;
        }
        const G4double rho = gas::MixtureDensity(mix) * g/cm3;
        auto* m = new G4Material(mix.name, rho,
                                 static_cast<G4int>(mix.parts.size()),
                                 kStateGas, kGasT, kGasP);
        for (std::size_t i = 0; i < mix.parts.size(); ++i)
            m->AddMaterial(pure.at(mix.parts[i].component),
                           gas::MassFraction(mix, i));
        fGasMaterials[mix.name] = m;
    }

    // ── He-3 at 300 bar ──────────────────────────────────────
    {
        auto* isoHe3 = new G4Isotope("He3_iso", 2, 3, 3.0160293*g/mole);
        auto* elHe3  = new G4Element("Helium3_elem", "3He", 1);
        elHe3->AddIsotope(isoHe3, 1.0);
        auto* m = new G4Material("He3Gas_300bar", 37.6e-3*g/cm3, 1,
                                  kStateGas, 293.15*kelvin, 300*atmosphere);
        m->AddElement(elHe3, 1);
        fGasMaterials["He3Gas_300bar"] = m;
    }

    // ── Structural / detector materials ─────────────────────
    {
        auto* m = new G4Material("CFRP", 1.55*g/cm3, 3);
        m->AddElement(elC, 0.8968); m->AddElement(elH, 0.0207); m->AddElement(elO, 0.0826);
        fGasMaterials["CFRP"] = m;
    }
    {
        auto* m = new G4Material("ResistivePaste", 1.4*g/cm3, 3);
        m->AddElement(elC, 0.65); m->AddElement(elH, 0.08); m->AddElement(elO, 0.27);
        fGasMaterials["ResistivePaste"] = m;
    }
    {
        auto* m = new G4Material("FR4", 1.85*g/cm3, 4);
        m->AddElement(elSi, 0.2805); m->AddElement(elO,  0.4195);
        m->AddElement(elC,  0.2750); m->AddElement(elH,  0.0250);
        fGasMaterials["FR4"] = m;
    }
    {
        auto* m = new G4Material("Rohacell51", 0.052*g/cm3, 4);
        m->AddElement(elC, 0.5783); m->AddElement(elH, 0.0602);
        m->AddElement(elN, 0.1687); m->AddElement(elO, 0.1928);
        fGasMaterials["Rohacell51"] = m;
    }
    {
        auto* m = new G4Material("LAB_LiqScint", 0.86*g/cm3, 2);
        m->AddElement(elC, 0.8780); m->AddElement(elH, 0.1220);
        fGasMaterials["LAB_LiqScint"] = m;
    }
    // BlackMylar: used for scint wall wrapping and back scint tape (G4_MYLAR = PET).
    fGasMaterials["BlackMylar"] = nist->FindOrBuildMaterial("G4_MYLAR");
}

// ============================================================
G4Material* DetectorConstruction::GetGasMixture(const std::string& name) {
    auto it = fGasMaterials.find(name);
    if (it == fGasMaterials.end())
        throw std::runtime_error("Unknown gas/material: " + name + "\n" +
                                 gas::ListMixtures());
    return it->second;
}

// ============================================================
G4VPhysicalVolume* DetectorConstruction::Construct() {
    DefineMaterials();

    if (fConfig.mode == SimMode::kP2Wedge)
        return ConstructP2();

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* matAir     = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* matMylar   = nist->FindOrBuildMaterial("G4_MYLAR");
    G4Material* matAl      = nist->FindOrBuildMaterial("G4_Al");
    G4Material* matKapton  = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Material* matCu      = nist->FindOrBuildMaterial("G4_Cu");
    G4Material* matSteel   = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* matGas     = GetGasMixture(fConfig.gas);
    G4Material* matHe3     = fGasMaterials.at("He3Gas_300bar");
    G4Material* matCFRP    = fGasMaterials.at("CFRP");
    G4Material* matResPaste= fGasMaterials.at("ResistivePaste");
    G4Material* matFR4     = fGasMaterials.at("FR4");
    G4Material* matRohacell= fGasMaterials.at("Rohacell51");
    G4Material* matLAB     = fGasMaterials.at("LAB_LiqScint");
    G4Material* matBlkMylar= fGasMaterials.at("BlackMylar");
    G4Material* matPlScint = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    G4double detXY = 40.0 * cm;

    // ── MM stack layer thicknesses ───────────────────────────
    G4double tMylar    = 40.0  * um;
    G4double tAlWin    = 0.1   * um;
    G4double tKapCath  = 50.0  * um;
    G4double tCuCath   = 9.0   * um;
    G4double tDrift    = 3.0   * cm;
    G4double tMesh     = 30.0  * um;
    G4double tAmp      = 150.0 * um;
    G4double tResPaste = 100.0 * um;
    G4double mmTotalZ  = tMylar + tAlWin + tKapCath + tCuCath
                       + tDrift + tMesh + tAmp + tResPaste;

    // ── PCB stack ─────────────────────────────────────────────
    G4double tPCB_Kap  = 50.0  * um;
    G4double tPCB_Cu   = 26.0  * um;
    G4double tPCB_FR4  = 100.0 * um;
    G4double tPCB_Roh  = 5.0   * mm;
    G4double tPCB_Al   = 50.0  * um;
    G4double pcbTotalZ = tPCB_Kap + 4*(tPCB_Cu + tPCB_FR4) + tPCB_Roh + tPCB_Al;

    // ── Scint wall (BlackMylar tape, from Full_Geant) ─────────
    G4double tBlkTape  = 200.0 * um;
    G4double tPlScint  = 3.0   * mm;
    G4double tScAl     = 50.0  * um;
    G4double scintWallZ = 2*tBlkTape + tPlScint + tScAl;

    // ── LS stack (from Full_Geant: 2 layers × 2cm with inner liners) ─────
    G4double tLSCfrp      = fConfig.cfrpThickness_mm  * mm;   // structural CFRP wall
    G4double tLSInnerCfrp = fConfig.ls_inner_cfrp_um  * um;   // inner CFRP liner
    G4double tLSInnerAl   = fConfig.ls_inner_al_um    * um;   // Al liner
    G4double tLS          = fConfig.ls_thick_cm        * cm;   // LAB layer
    // 3 CFRP walls + 2 inner CFRP liners + 2 Al liners + 2 LAB layers
    G4double lsStackZ = 3*tLSCfrp + 2*(tLSInnerCfrp + tLSInnerAl + tLS);

    // ── Vis attributes ───────────────────────────────────────
    auto visMylar     = new G4VisAttributes(G4Color(0.7, 0.9, 0.7, 0.5));
    auto visAl        = new G4VisAttributes(G4Color(0.7, 0.7, 0.7, 0.8));
    auto visKapton    = new G4VisAttributes(G4Color(0.9, 0.7, 0.0, 0.7));
    auto visCu        = new G4VisAttributes(G4Color(0.8, 0.4, 0.1, 0.8));
    auto visDrift     = new G4VisAttributes(G4Color(0.2, 0.5, 1.0, 0.3));
    auto visMesh      = new G4VisAttributes(G4Color(0.5, 0.5, 0.5, 0.9));
    auto visAmp       = new G4VisAttributes(G4Color(1.0, 0.3, 0.3, 0.3));
    auto visResPaste  = new G4VisAttributes(G4Color(0.2, 0.2, 0.2, 0.8));
    auto visHe3       = new G4VisAttributes(G4Color(0.6, 0.9, 1.0, 0.4));
    auto visCFRP      = new G4VisAttributes(G4Color(0.15, 0.15, 0.15, 0.9));
    auto visFR4       = new G4VisAttributes(G4Color(0.2, 0.6, 0.2, 0.8));
    auto visRohacell  = new G4VisAttributes(G4Color(0.9, 0.9, 0.6, 0.5));
    auto visScint     = new G4VisAttributes(G4Color(0.9, 0.9, 0.2, 0.7));
    auto visLAB       = new G4VisAttributes(G4Color(0.3, 0.8, 0.9, 0.4));
    auto visBlkMylar  = new G4VisAttributes(G4Color(0.1, 0.1, 0.1, 0.9));

    G4LogicalVolume*   worldLV = nullptr;
    G4VPhysicalVolume* worldPV = nullptr;

    // ── PlaceSlab helper ─────────────────────────────────────
    // Placing slab with its front face at *zFront, centred at zFront+t/2.
    // Increments zFront by t. hx,hy may differ from detXY/2 for larger detectors.
    auto MakeSlab = [&](const std::string& name, G4double t,
                         G4double hx, G4double hy,
                         G4Material* mat, G4VisAttributes* vis) -> G4LogicalVolume* {
        auto* box = new G4Box(name, hx, hy, t/2);
        auto* lv  = new G4LogicalVolume(box, mat, name);
        if (vis) lv->SetVisAttributes(vis);
        return lv;
    };

    // ═══════════════════════════════════════════════════════════
    // VACUUM MODE
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kVacuum) {

        G4double alThickness = fConfig.alThickness_mm * mm;
        G4double alGap       = 2.0 * cm;

        G4double gunZ = 10.0 * cm;
        G4double minUpstream = 2.0*gunZ - mmTotalZ - 2.5*cm;
        G4double upstreamMargin = std::max(alGap + alThickness + 0.5*cm, minUpstream);
        G4double worldZ = mmTotalZ + upstreamMargin + 2.5*cm;

        auto* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, worldZ/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        if (alThickness > 0) {
            G4double alZCenter = -mmTotalZ/2.0 - alGap - alThickness/2.0;
            auto* alBox = new G4Box("AlShield", detXY/2, detXY/2, alThickness/2);
            auto* alLV  = new G4LogicalVolume(alBox, matAl, "AlShield");
            auto* visAlS = new G4VisAttributes(G4Color(0.75, 0.75, 0.75, 0.9));
            visAlS->SetForceSolid(true);
            alLV->SetVisAttributes(visAlS);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,alZCenter),
                              alLV, "AlShield", worldLV, false, 0, true);
        }

        G4cout << "\n=== Vacuum mode ===" << G4endl;
        G4cout << "  Gas: " << fConfig.gas
               << "  (rho=" << matGas->GetDensity()/(mg/cm3) << " mg/cm3)" << G4endl;

        // Place MM layers
        G4double zFront = -mmTotalZ / 2.0;
        auto PlaceSlab = [&](const std::string& name, G4double t,
                              G4Material* mat, G4VisAttributes* vis,
                              G4LogicalVolume*& outLV) {
            auto* lv = MakeSlab(name, t, detXY/2, detXY/2, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zFront+t/2), lv, name, worldLV, false, 0, true);
            zFront += t;
            outLV = lv;
        };
        G4LogicalVolume* dummy = nullptr;
        PlaceSlab("GasWindow_Mylar",     tMylar,    matMylar,    visMylar,    dummy);
        PlaceSlab("GasWindow_Al",        tAlWin,    matAl,       visAl,       dummy);
        PlaceSlab("DriftCathode_Kapton", tKapCath,  matKapton,   visKapton,   dummy);
        PlaceSlab("DriftCathode_Cu",     tCuCath,   matCu,       visCu,       dummy);
        PlaceSlab("DriftGas",            tDrift,    matGas,      visDrift,    fDriftGasLV);
        PlaceSlab("Micromesh",           tMesh,     matSteel,    visMesh,     dummy);
        PlaceSlab("AmpGas",              tAmp,      matGas,      visAmp,      fAmpGasLV);
        PlaceSlab("ResistivePaste",      tResPaste, matResPaste, visResPaste, dummy);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // FULL EXPERIMENT MODE
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kFullExperiment) {

        // He-3 capsule (from Full_Geant: r=1.5 cm, L=8 cm)
        G4double he3R          = 1.5  * cm;
        G4double he3HalfL      = 4.0  * cm;
        G4double alWallT       = 0.5  * mm;
        G4double cfrpWallT     = 0.9  * mm;
        G4double alR           = he3R + alWallT;
        G4double cfrpR         = alR  + cfrpWallT;
        G4double alHalfL       = he3HalfL + alWallT;
        G4double cfrpHalfL     = alHalfL  + cfrpWallT;
        G4double capsuleZExtent = 2.0 * cfrpR;

        G4double airGap1   = 200.0 * mm;
        G4double airGap2   = 20.0  * mm;
        G4double airGap3   = 20.0  * mm;

        G4double totalFullZ = capsuleZExtent + airGap1 + mmTotalZ + pcbTotalZ
                            + airGap2 + scintWallZ + airGap3 + lsStackZ;

        auto* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, (totalFullZ+2*cm)/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4double capsuleZCenter = -totalFullZ/2.0 + cfrpR;
        fHe3GasCenterZ = capsuleZCenter;

        auto* capRot = new G4RotationMatrix();
        capRot->rotateX(-90.*deg);

        auto* cfrpSolid = new G4Tubs("He3Capsule_CFRP", 0, cfrpR, cfrpHalfL, 0, 360.*deg);
        auto* cfrpLV    = new G4LogicalVolume(cfrpSolid, matCFRP, "He3Capsule_CFRP");
        cfrpLV->SetVisAttributes(visCFRP);
        new G4PVPlacement(capRot, G4ThreeVector(0,0,capsuleZCenter), cfrpLV, "He3Capsule_CFRP", worldLV, false, 0, true);

        auto* alSolid = new G4Tubs("He3Capsule_Al", 0, alR, alHalfL, 0, 360.*deg);
        auto* alLV    = new G4LogicalVolume(alSolid, matAl, "He3Capsule_Al");
        alLV->SetVisAttributes(visAl);
        new G4PVPlacement(nullptr, G4ThreeVector(), alLV, "He3Capsule_Al", cfrpLV, false, 0, true);

        auto* he3Solid = new G4Tubs("He3Gas", 0, he3R, he3HalfL, 0, 360.*deg);
        fHe3GasLV = new G4LogicalVolume(he3Solid, matHe3, "He3Gas");
        fHe3GasLV->SetVisAttributes(visHe3);
        new G4PVPlacement(nullptr, G4ThreeVector(), fHe3GasLV, "He3Gas", alLV, false, 0, true);

        G4cout << "\n=== Full-experiment geometry ===" << G4endl;
        G4cout << "  He-3: r=" << he3R/cm << " cm, L=" << 2*he3HalfL/cm << " cm, 300 bar" << G4endl;
        G4cout << "  LS stack (2×" << fConfig.ls_thick_cm << " cm LAB, "
               << fConfig.cfrpThickness_mm << " mm CFRP walls)" << G4endl;
        G4cout << "  Total Z: " << totalFullZ/mm << " mm" << G4endl;

        G4double zF = -totalFullZ/2.0 + capsuleZExtent;
        auto Place = [&](const std::string& name, G4double t,
                          G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
            out = MakeSlab(name, t, detXY/2, detXY/2, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
            zF += t;
        };
        G4LogicalVolume* dL = nullptr;

        Place("AirGap1",             airGap1,   matAir,      nullptr,      dL);
        Place("GasWindow_Mylar",     tMylar,    matMylar,    visMylar,     dL);
        Place("GasWindow_Al",        tAlWin,    matAl,       visAl,        dL);
        Place("DriftCathode_Kapton", tKapCath,  matKapton,   visKapton,    dL);
        Place("DriftCathode_Cu",     tCuCath,   matCu,       visCu,        dL);
        Place("DriftGas",            tDrift,    matGas,      visDrift,     fDriftGasLV);
        Place("Micromesh",           tMesh,     matSteel,    visMesh,      dL);
        Place("AmpGas",              tAmp,      matGas,      visAmp,       fAmpGasLV);
        Place("ResistivePaste",      tResPaste, matResPaste, visResPaste,  dL);

        Place("PCB_Kapton",   tPCB_Kap, matKapton,  visKapton,  dL);
        Place("PCB_Cu_1",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_1",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_2",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_2",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_3",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_3",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_4",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_4",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Rohacell", tPCB_Roh, matRohacell,visRohacell,dL);
        Place("PCB_AlFoil",   tPCB_Al,  matAl,      visAl,      dL);

        Place("AirGap2",               airGap2,   matAir,      nullptr,     dL);
        Place("ScintWall_BlackTape1",  tBlkTape,  matBlkMylar, visBlkMylar, dL);
        Place("PlasticScint",          tPlScint,  matPlScint,  visScint,    dL);
        Place("ScintWall_BlackTape2",  tBlkTape,  matBlkMylar, visBlkMylar, dL);
        Place("ScintWall_AlFoil",      tScAl,     matAl,       visAl,       dL);

        Place("AirGap3",          airGap3,    matAir, nullptr, dL);
        Place("LS_CFRP_1",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_1",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_1",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_1",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_2",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_2",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_2",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_2",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_3",        tLSCfrp,    matCFRP, visCFRP, dL);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // SR-90 CALIBRATION MODE
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kSr90Calibration) {

        G4double airToMM  = 226.5 * mm;
        G4double airGap2  = 20.0  * mm;
        G4double airGap3  = 20.0  * mm;

        G4double totalZ = airToMM + mmTotalZ + pcbTotalZ
                        + airGap2 + scintWallZ + airGap3 + lsStackZ;

        auto* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, (totalZ+2*cm)/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        fHe3GasCenterZ = -totalZ / 2.0;

        G4cout << "\n=== Sr-90 calibration geometry ===" << G4endl;
        G4cout << "  Air source-to-MM: " << airToMM/mm << " mm" << G4endl;
        G4cout << "  LS (2×" << fConfig.ls_thick_cm << " cm, CFRP " << fConfig.cfrpThickness_mm << " mm)" << G4endl;
        G4cout << "  Total Z: " << totalZ/mm << " mm" << G4endl;

        G4double zF = -totalZ / 2.0;
        auto Place = [&](const std::string& name, G4double t,
                          G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
            out = MakeSlab(name, t, detXY/2, detXY/2, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
            zF += t;
        };
        G4LogicalVolume* dL = nullptr;

        Place("AirGap1",             airToMM,   matAir,      nullptr,      dL);
        Place("GasWindow_Mylar",     tMylar,    matMylar,    visMylar,     dL);
        Place("GasWindow_Al",        tAlWin,    matAl,       visAl,        dL);
        Place("DriftCathode_Kapton", tKapCath,  matKapton,   visKapton,    dL);
        Place("DriftCathode_Cu",     tCuCath,   matCu,       visCu,        dL);
        Place("DriftGas",            tDrift,    matGas,      visDrift,     fDriftGasLV);
        Place("Micromesh",           tMesh,     matSteel,    visMesh,      dL);
        Place("AmpGas",              tAmp,      matGas,      visAmp,       fAmpGasLV);
        Place("ResistivePaste",      tResPaste, matResPaste, visResPaste,  dL);

        Place("PCB_Kapton",   tPCB_Kap, matKapton,  visKapton,  dL);
        Place("PCB_Cu_1",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_1",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_2",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_2",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_3",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_3",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Cu_4",     tPCB_Cu,  matCu,      visCu,      dL);
        Place("PCB_FR4_4",    tPCB_FR4, matFR4,     visFR4,     dL);
        Place("PCB_Rohacell", tPCB_Roh, matRohacell,visRohacell,dL);
        Place("PCB_AlFoil",   tPCB_Al,  matAl,      visAl,      dL);

        Place("AirGap2",               airGap2,  matAir,      nullptr,     dL);
        Place("ScintWall_BlackTape1",  tBlkTape, matBlkMylar, visBlkMylar, dL);
        Place("PlasticScint",          tPlScint, matPlScint,  visScint,    dL);
        Place("ScintWall_BlackTape2",  tBlkTape, matBlkMylar, visBlkMylar, dL);
        Place("ScintWall_AlFoil",      tScAl,    matAl,       visAl,       dL);

        Place("AirGap3",          airGap3,    matAir, nullptr, dL);
        Place("LS_CFRP_1",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_1",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_1",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_1",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_2",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_2",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_2",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_2",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_3",        tLSCfrp,    matCFRP, visCFRP, dL);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // SR-90 NO-MM MODE
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kSr90NoMM) {

        G4double airToScint = 226.5*mm + mmTotalZ + pcbTotalZ + 20.0*mm;
        G4double airGap3    = 20.0 * mm;

        G4double totalZ = airToScint + scintWallZ + airGap3 + lsStackZ;

        auto* worldSolid = new G4Box("World", detXY/2+2*cm, detXY/2+2*cm, (totalZ+2*cm)/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        fHe3GasCenterZ = -totalZ / 2.0;

        G4cout << "\n=== Sr-90 no-MM geometry ===" << G4endl;
        G4cout << "  Air source-to-scint: " << airToScint/mm << " mm" << G4endl;
        G4cout << "  Total Z: " << totalZ/mm << " mm" << G4endl;

        G4double zF = -totalZ / 2.0;
        auto Place = [&](const std::string& name, G4double t,
                          G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
            out = MakeSlab(name, t, detXY/2, detXY/2, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
            zF += t;
        };
        G4LogicalVolume* dL = nullptr;

        Place("AirGap1",              airToScint, matAir,      nullptr,     dL);
        Place("ScintWall_BlackTape1", tBlkTape,   matBlkMylar, visBlkMylar, dL);
        Place("PlasticScint",         tPlScint,   matPlScint,  visScint,    dL);
        Place("ScintWall_BlackTape2", tBlkTape,   matBlkMylar, visBlkMylar, dL);
        Place("ScintWall_AlFoil",     tScAl,      matAl,       visAl,       dL);

        Place("AirGap3",          airGap3,    matAir, nullptr, dL);
        Place("LS_CFRP_1",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_1",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_1",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_1",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_2",        tLSCfrp,    matCFRP, visCFRP, dL);
        Place("LS_InnerCFRP_2",   tLSInnerCfrp, matCFRP, visCFRP, dL);
        Place("LS_Al_2",          tLSInnerAl, matAl,   visAl,   dL);
        Place("LiqScint_2",       tLS,        matLAB,  visLAB,  dL);
        Place("LS_CFRP_3",        tLSCfrp,    matCFRP, visCFRP, dL);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // LS CALIBRATION MODE
    // Bare source gun → air gap → 1 LS layer only.
    // No source capsule, no back scint.
    // ═══════════════════════════════════════════════════════════
    if (fConfig.mode == SimMode::kLSCalib) {

        G4double airToLS  = fConfig.source_to_det_mm * mm;
        // LS: CFRP_front | InnerCFRP | Al | LAB | InnerCFRP | Al | CFRP_back
        G4double lsZ    = 2*tLSCfrp + 2*(tLSInnerCfrp + tLSInnerAl) + tLS;
        G4double totalZ = airToLS + lsZ;

        G4double lsHX = 22.5*cm;  // 45×45 cm LS face
        G4double lsHY = 22.5*cm;

        auto* worldSolid = new G4Box("World", lsHX+2*cm, lsHY+2*cm, (totalZ+2*cm)/2);
        worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
        worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
        worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

        fHe3GasCenterZ = -totalZ / 2.0;  // gun at front of world

        G4cout << "\n=== LS Calibration geometry ===" << G4endl;
        G4cout << "  Source-to-LS air gap: " << airToLS/mm << " mm" << G4endl;
        G4cout << "  LS: " << tLS/cm << " cm LAB,  CFRP walls " << fConfig.cfrpThickness_mm << " mm" << G4endl;
        G4cout << "  Total Z: " << totalZ/mm << " mm" << G4endl;

        G4double zF = -totalZ / 2.0;

        auto Place = [&](const std::string& name, G4double t,
                          G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
            out = MakeSlab(name, t, lsHX, lsHY, mat, vis);
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
            zF += t;
        };
        G4LogicalVolume* dL = nullptr;

        Place("AirGap1",       airToLS,      matAir,  nullptr,  dL);
        Place("LS_CFRP_1",     tLSCfrp,      matCFRP, visCFRP,  dL);
        Place("LS_InnerCFRP_1",tLSInnerCfrp, matCFRP, visCFRP,  dL);
        Place("LS_Al_1",       tLSInnerAl,   matAl,   visAl,    dL);
        Place("LiqScint_1",    tLS,          matLAB,  visLAB,   dL);
        Place("LS_InnerCFRP_2",tLSInnerCfrp, matCFRP, visCFRP,  dL);
        Place("LS_Al_2",       tLSInnerAl,   matAl,   visAl,    dL);
        Place("LS_CFRP_2",     tLSCfrp,      matCFRP, visCFRP,  dL);

        return worldPV;
    }

    // ═══════════════════════════════════════════════════════════
    // BACK SCINT CALIBRATION MODE
    // Bare source gun → air gap → 1 back scint bar only.
    // No source capsule, no LS layer.
    // ═══════════════════════════════════════════════════════════
    // else kBackScintCalib

    G4double airToBSc = fConfig.source_to_det_mm * mm;
    G4double tBscTape = fConfig.backscint_tape_um  * um;
    G4double tBscAl   = fConfig.backscint_al_um    * um;
    G4double tBscPVT  = fConfig.backscint_thick_cm * cm;
    G4double bscZ     = 2*tBscTape + 2*tBscAl + tBscPVT;
    G4double totalZ   = airToBSc + bscZ;

    G4double bsHX = fConfig.backscint_u_cm/2 * cm;  // 15 cm half-width (30 cm total)
    G4double bsHY = fConfig.backscint_v_cm/2 * cm;  // 10 cm half-height (20 cm total)

    auto* worldSolid = new G4Box("World", bsHX+2*cm, bsHY+2*cm, (totalZ+2*cm)/2);
    worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
    worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    fHe3GasCenterZ = -totalZ / 2.0;

    G4cout << "\n=== Back Scint Calibration geometry ===" << G4endl;
    G4cout << "  Source-to-scint air gap: " << airToBSc/mm << " mm" << G4endl;
    G4cout << "  Back scint: " << tBscPVT/cm << " cm PVT, "
           << fConfig.backscint_u_cm << "×" << fConfig.backscint_v_cm << " cm face" << G4endl;
    G4cout << "  Total Z: " << totalZ/mm << " mm" << G4endl;

    G4double zF = -totalZ / 2.0;

    auto visBscPVT = new G4VisAttributes(G4Color(0.9, 0.5, 0.1, 0.8));

    auto PlaceB = [&](const std::string& name, G4double t,
                       G4Material* mat, G4VisAttributes* vis, G4LogicalVolume*& out) {
        out = MakeSlab(name, t, bsHX, bsHY, mat, vis);
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+t/2), out, name, worldLV, false, 0, true);
        zF += t;
    };
    G4LogicalVolume* dL = nullptr;

    {
        // Air gap: use world-sized slab so gun position is inside air
        auto* airBox = new G4Box("AirGap1", bsHX+2*cm, bsHY+2*cm, airToBSc/2);
        auto* airLV  = new G4LogicalVolume(airBox, matAir, "AirGap1");
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,zF+airToBSc/2), airLV, "AirGap1", worldLV, false, 0, true);
        zF += airToBSc;
    }

    PlaceB("BackScintWrap_Tape1", tBscTape, matBlkMylar, visBlkMylar, dL);
    PlaceB("BackScintWrap_Al1",   tBscAl,   matAl,       visAl,       dL);
    PlaceB("BackScint",           tBscPVT,  matPlScint,  visBscPVT,   fBackScintLV);
    PlaceB("BackScintWrap_Al2",   tBscAl,   matAl,       visAl,       dL);
    PlaceB("BackScintWrap_Tape2", tBscTape, matBlkMylar, visBlkMylar, dL);

    return worldPV;
}

// ============================================================
// P2 WEDGE MODE
//
// Full P2 wedge Micromegas with its gas envelope, in gerber coordinates:
// the wedge apex (= MESA beam axis) is at world (x,y) = (0,0), and all
// layers are placed at their true transverse position. z = 0 is the FRONT
// WINDOW plane (top of the gas frame); the beam travels along +z:
//
//   front window: 10 um mylar, bulged (terraced dome, apex upstream) z < 0
//   front gas gap                        z = 0 .. 3.879 mm
//   drift cathode: one 120 um mylar foil + 1 um Al facing the drift gas,
//     lying on the frame ledge (Alexandra 2026-10-06)
//   drift gas   [sensitive]              3.964 mm (4 mm frame - mesh)
//   micromesh   (woven SS 45 um opening / 18 um wire -> 36 um slab of
//     effective-density steel; areal mass matches the woven mesh)
//   amplification gas [sensitive]        150 um, with the Dynamask pillars
//   readout PCB: F.Cu 18um / FR4 200um / B.Cu 18um   (Stack_Up_P2.txt);
//     real F.Cu pads, B.Cu density-scaled by radial band
//   back gas gap (= carbon back-frame depth, 1 mm), pad-zone outline
//   back window: 10 um mylar, bulged (dome apex downstream)
//   + frame rings: front = plastic 8 mm frame (P2_Frame_V2.0.stp, two
//     rings either side of the foil ledge) standing on the 150 um bulk
//     Dynamask; back = 1 mm carbon frame, board outline -> pad zone.
//
// The window bulge from the few-mbar overpressure is modelled as a
// terraced dome: N stacked gas prisms whose wedge profile shrinks about
// the window centroid following a spherical-cap profile, each step capped
// by a flat mylar ring, so every vertical path crosses exactly one mylar
// thickness. See docs/P2_MODEL.md for the assumptions table.
// ============================================================
G4VPhysicalVolume* DetectorConstruction::ConstructP2() {
    using namespace P2;

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* matAir   = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* matMylar = nist->FindOrBuildMaterial("G4_MYLAR");
    G4Material* matAl    = nist->FindOrBuildMaterial("G4_Al");
    G4Material* matCu    = nist->FindOrBuildMaterial("G4_Cu");
    G4Material* matSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* matPlast = nist->FindOrBuildMaterial("G4_POLYCARBONATE");
    G4Material* matFR4   = fGasMaterials.at("FR4");
    G4Material* matGas   = GetGasMixture(fConfig.gas);

    // Carbon-fiber back frame (confirmed carbon; 1.6 g/cm3 typical CF plate)
    G4Material* matCarbon = G4Material::GetMaterial("CarbonFrame", false);
    if (!matCarbon) {
        matCarbon = new G4Material("CarbonFrame", 1.60*g/cm3, 1);
        matCarbon->AddMaterial(nist->FindOrBuildMaterial("G4_C"), 1.0);
    }

    // Pillars: Dynamask photoimageable dry film, rho 1.2-1.4 g/cm3. An
    // acrylic, so PMMA's C5H8O2 composition at the middle of that density
    // range stands in for the (unpublished) exact formulation.
    G4Material* matDynamask = G4Material::GetMaterial("Dynamask", false);
    if (!matDynamask) {
        matDynamask = new G4Material("Dynamask", 1.30*g/cm3, 1);
        matDynamask->AddMaterial(nist->FindOrBuildMaterial("G4_PLEXIGLASS"), 1.0);
    }

    // ── Thicknesses (config in mm/um; P2:: constants in mm) ──────────────
    const G4double tCuF     = kTCuF * mm;
    const G4double tFR4     = kTFR4 * mm;
    const G4double tCuB     = kTCuB * mm;
    const G4double tAmp     = fConfig.p2_amp_um       * um;
    const G4double tDrift   = fConfig.p2_drift_mm     * mm;
    const G4double tCathMy  = fConfig.p2_cath_mylar_um* um;
    const G4double tCathAl  = fConfig.p2_cath_al_um   * um;
    const G4double cathGap  = fConfig.p2_cath_gap_mm  * mm;
    const G4double frontGap = fConfig.p2_front_gap_mm * mm;

    // Woven mesh -> slab of effective-density steel: slab thickness = 2 wire
    // diameters (weave height); steel volume per unit area of a plain weave
    // with wire d, pitch p (two orthogonal wire sets) = pi*d^2/(2p), so the
    // fill fraction of the slab is pi*d/(4p). Optical transparency
    // (opening/pitch)^2 ~ 51% is NOT modelled — only the material budget.
    const G4double dWire    = fConfig.p2_mesh_wire_um * um;
    const G4double meshPitch= dWire + fConfig.p2_mesh_open_um * um;
    const G4double tMesh    = 2.0 * dWire;
    const G4double meshFill = M_PI * dWire / (4.0 * meshPitch);
    G4Material* matMesh = G4Material::GetMaterial("MeshSteelEff", false);
    if (!matMesh) {
        matMesh = new G4Material("MeshSteelEff",
                                 matSteel->GetDensity() * meshFill, 1);
        matMesh->AddMaterial(matSteel, 1.0);
    }
    const G4double backGap  = fConfig.p2_back_gap_mm  * mm;
    const G4double tWin     = fConfig.p2_window_um    * um;
    const G4double hFront   = fConfig.p2_bulge_front_mm * mm;
    const G4double hBack    = fConfig.p2_bulge_back_mm  * mm;

    // Density-scaled copper: a full-thickness sheet whose density carries a
    // measured area coverage. Cached by name so bands sharing a coverage
    // share a material. Never exactly zero — G4Material rejects that.
    auto EffCu = [&](double frac) {
        char nm[32];
        std::snprintf(nm, sizeof nm, "CuEff_%06d",
                      static_cast<int>(std::lround(frac * 1e5)));
        G4Material* m = G4Material::GetMaterial(nm, false);
        if (!m) {
            m = new G4Material(nm, matCu->GetDensity() * std::max(frac, 1e-6), 1);
            m->AddMaterial(matCu, 1.0);
        }
        return m;
    };

    // ── Profiles (WedgeOutline returns mm-valued coords = G4 length) ─────
    const auto polyBoard   = WedgeOutline(kBoardRIn, kBoardROut,
                                          kBoardEdgeOffset, kBoardTopCut);
    const auto polyFrame   = WedgeOutline(kFrameRIn, kFrameROut,
                                          kFrameEdgeOffset, kFrameTopCut);
    const auto polyOpening = WedgeOutline(kOpenRIn, kOpenROut,
                                          kOpenEdgeOffset, kOpenTopCut);
    const auto polyDriftOpening = WedgeOutline(kDriftOpenRIn, kDriftOpenROut,
                                               kDriftOpenEdgeOffset, kBoardTopCut);
    const auto polyZone    = WedgeOutline(kZoneRIn, kZoneROut,
                                          kZoneEdgeOffset, kBoardTopCut);
    const G4TwoVector openC = Centroid(polyOpening);
    const G4TwoVector zoneC = Centroid(polyZone);

    // ── Vis ───────────────────────────────────────────────────────────────
    auto visMylar = new G4VisAttributes(G4Color(0.70, 0.90, 0.70, 0.50));
    auto visWin   = new G4VisAttributes(G4Color(0.55, 0.85, 0.95, 0.45));
    auto visAl    = new G4VisAttributes(G4Color(0.70, 0.70, 0.70, 0.80));
    auto visCu    = new G4VisAttributes(G4Color(0.80, 0.40, 0.10, 0.80));
    auto visDrift = new G4VisAttributes(G4Color(0.20, 0.50, 1.00, 0.30));
    auto visMesh  = new G4VisAttributes(G4Color(0.50, 0.50, 0.50, 0.90));
    auto visAmp   = new G4VisAttributes(G4Color(1.00, 0.30, 0.30, 0.30));
    auto visFR4   = new G4VisAttributes(G4Color(0.20, 0.60, 0.20, 0.80));
    auto visFrame = new G4VisAttributes(G4Color(0.90, 0.90, 0.85, 0.90));
    auto visGasV  = new G4VisAttributes(G4Color(0.55, 0.85, 0.95, 0.15));

    // ── World: contains the full wedge in gerber coordinates ─────────────
    // cathGap <= 0 selects the single-foil cathode (the baseline since
    // 2026-10-06); > 0 the old two-foil model with gas between the foils.
    const bool     twoFoils = cathGap > 0.0;
    const G4double tCath    = (twoFoils ? tCathMy + cathGap : 0.0) + tCathMy + tCathAl;
    // The frame steps out at the drift-foil ledge (P2_Frame_V2.0.stp), so it
    // is built as two rings: window side down to the foil's board face, and
    // board side from there to the Dynamask. The frame sits on the 150 um
    // Dynamask left on the PCB by the bulk process, and the mesh is cut
    // inside it (Alexandra 2026-10-06), so the board-side ring spans drift
    // gas + mesh only: 3.964 + 0.036 = 4.000 mm, the STEP ledge height.
    const G4double frameHWin   = frontGap + tCath;
    const G4double frameHBoard = tDrift + tMesh;
    const G4double frameH      = frameHWin + frameHBoard + tAmp;   // to the pad plane
    const G4double zPCBEnd  = frameH + tCuF + tFR4 + tCuB;
    const G4double zBackWin = zPCBEnd + backGap;
    const G4double zMax     = zBackWin + hBack + tWin;
    fP2FrontZ = -(hFront + tWin);

    const G4double worldHX = (kBoardROut + 40.0) * mm;
    const G4double worldHY = (kBoardTopCut + 70.0) * mm;
    const G4double worldHZ = std::max(hFront + tWin + 80.0*mm, zMax + 20.0*mm);

    auto* worldSolid = new G4Box("World", worldHX, worldHY, worldHZ);
    auto* worldLV = new G4LogicalVolume(worldSolid, matAir, "World");
    auto* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV,
                                      "World", nullptr, false, 0, true);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ── Helpers ───────────────────────────────────────────────────────────
    auto Prism = [](const std::string& name,
                    const std::vector<G4TwoVector>& poly, G4double t) {
        return new G4ExtrudedSolid(name, poly, t/2,
                                   G4TwoVector(), 1.0, G4TwoVector(), 1.0);
    };

    G4double zF = 0.0;   // running front face, starts at the front window plane
    auto PlaceLayer = [&](const std::string& name,
                          const std::vector<G4TwoVector>& poly, G4double t,
                          G4Material* mat, G4VisAttributes* vis,
                          G4LogicalVolume** outLV = nullptr) {
        auto* lv = new G4LogicalVolume(Prism(name, poly, t), mat, name);
        if (vis) lv->SetVisAttributes(vis);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zF + t/2), lv, name,
                          worldLV, false, 0, true);
        zF += t;
        if (outLV) *outLV = lv;
    };

    // Terraced bulged window. zBase: window plane; sign: -1 = dome rises
    // upstream (front window), +1 = downstream (back window); H: sag height.
    auto BuildWindow = [&](const std::string& tag, G4double zBase,
                           int sign, G4double H,
                           const std::vector<G4TwoVector>& outline,
                           const G4TwoVector& c) {
        const int N = 6;
        std::vector<double> sig(N + 1), h(N + 1);
        for (int k = 0; k <= N; ++k) {
            const double a = k * M_PI / (2.0 * N);
            sig[k] = std::max(std::cos(a), 0.04);   // spherical-cap profile
            h[k]   = H * std::sin(a);
        }
        for (int k = 0; k < N; ++k) {
            const auto poly = ScaleAbout(outline, c, sig[k]);
            // gas step
            const G4double dz = h[k+1] - h[k];
            const std::string gnm = tag + "_Gas" + std::to_string(k);
            auto* glv = new G4LogicalVolume(Prism(gnm, poly, dz), matGas, gnm);
            glv->SetVisAttributes(visGasV);
            new G4PVPlacement(nullptr,
                G4ThreeVector(0, 0, zBase + sign * (h[k] + dz/2)),
                glv, gnm, worldLV, false, 0, true);
            // mylar terrace: annular ring on top of this step, full cap on the last
            const std::string mnm = tag + "_Mylar" + std::to_string(k);
            G4VSolid* msolid = nullptr;
            if (k < N - 1) {
                auto inner = ScaleAbout(outline, c, sig[k+1]);
                msolid = new G4SubtractionSolid(mnm,
                            Prism(mnm + "_o", poly, tWin),
                            Prism(mnm + "_i", inner, 4*tWin));
            } else {
                msolid = Prism(mnm, poly, tWin);
            }
            auto* mlv = new G4LogicalVolume(msolid, matMylar, mnm);
            mlv->SetVisAttributes(visWin);
            new G4PVPlacement(nullptr,
                G4ThreeVector(0, 0, zBase + sign * (h[k+1] + tWin/2)),
                mlv, mnm, worldLV, false, 0, true);
        }
    };

    // ── Stack, front window plane (z=0) downstream ────────────────────────
    // Drift cathode: one aluminized mylar foil lying on the frame ledge, Al
    // facing the drift gas. The two-foil variant keeps its old volume names.
    PlaceLayer("FrontGas",           polyOpening, frontGap, matGas,   visGasV);
    if (twoFoils) {
        PlaceLayer("DriftCathode_Mylar1",polyOpening, tCathMy,  matMylar, visMylar);
        PlaceLayer("DriftCathode_Gas",   polyOpening, cathGap,  matGas,   visGasV);
        PlaceLayer("DriftCathode_Mylar2",polyOpening, tCathMy,  matMylar, visMylar);
    } else {
        PlaceLayer("DriftCathode_Mylar", polyOpening, tCathMy,  matMylar, visMylar);
    }
    PlaceLayer("DriftCathode_Al",    polyOpening, tCathAl,  matAl,    visAl);
    G4LogicalVolume* driftLV = nullptr;
    // zF is the running front face, so the drift mid-plane is fixed before the
    // call advances it. The angled gun (--gun-theta) aims here rather than at
    // the window, so a tilted beam illuminates the same pads as a normal one
    // instead of walking across the plane with theta.
    fP2DriftCenterZ = zF + tDrift/2;
    // Drift frame for Stage B. World z is NOT the drift depth -- the drift
    // gas sits at positive world z, several mm downstream of the window, and
    // electrons drift toward +z (the mesh is downstream of the gas). A
    // digitizer that treated z as "distance drifted" would be wrong by an
    // offset AND a sign, and would fail quietly: MX17 documents exactly this
    // trap in response/digitizer/clusters.py, where half the gap came out at
    // negative depth. So record the frame here and let Stage B convert:
    //     depth = fP2MeshZ - z_world     (mm, 0 at the mesh, +tDrift at the cathode)
    fP2DriftEntryZ = zF;                // cathode side of the drift gas
    PlaceLayer("DriftGas",          polyDriftOpening, tDrift, matGas, visDrift, &driftLV);
    fDriftGasLV = driftLV;
    fP2MeshZ = zF;                      // drift gas ends / mesh begins
    PlaceLayer("Micromesh",         polyDriftOpening, tMesh,  matMesh, visMesh);
    G4LogicalVolume* ampLV = nullptr;
    fP2AmpEntryZ = zF;                  // mesh ends / amplification gap begins
    PlaceLayer("AmpGas",            polyDriftOpening, tAmp,   matGas,  visAmp, &ampLV);
    fP2PadPlaneZ = zF;                  // amp gap ends / pad plane
    fAmpGasLV = ampLV;
    // ── Mesh-support pillars, daughters of the amp gas ────────────────────
    // Full amp-gap height, so the gas they displace leaves the sensitive
    // volume. AmpGas is an extrusion of the absolute gerber outline placed at
    // x = y = 0, so its local x/y are gerber x/y. include/P2Pillars.hh is
    // generated with every pillar checked inside the outline and none
    // overlapping (scripts/gerber/extract_pillars.py), which is why the
    // ~12k-41k placements skip the per-volume overlap check.
    int nPillars = 0;
    if (fConfig.p2_pillars != "none") {
        const bool cern = fConfig.p2_pillars == "cern";
        auto PlacePillars = [&](const char* name, double r,
                                const PillarXY* xy, int n) {
            auto* lv = new G4LogicalVolume(
                new G4Tubs(name, 0, r*mm, tAmp/2, 0, 360.*deg), matDynamask, name);
            lv->SetVisAttributes(visFR4);
            for (int i = 0; i < n; ++i)
                new G4PVPlacement(nullptr, G4ThreeVector(xy[i].x*mm, xy[i].y*mm, 0),
                                  lv, name, ampLV, false, i, false);
            nPillars += n;
        };
        if (cern) {
            PlacePillars("Pillar",    kPillarRCern,    kPillarsCern,    kNPillarsCern);
            PlacePillars("PillarBig", kPillarRBigCern, kPillarsBigCern, kNPillarsBigCern);
        } else {
            PlacePillars("Pillar",    kPillarRSaclay,    kPillarsSaclay,    kNPillarsSaclay);
            PlacePillars("PillarBig", kPillarRBigSaclay, kPillarsBigSaclay, kNPillarsBigSaclay);
        }
    }
    // ── Readout copper ────────────────────────────────────────────────────
    // Each 18 um copper layer is built as a GAS-FILLED envelope with the
    // copper placed inside it, rather than as one board-spanning sheet of
    // reduced-density copper. Two reasons, and the second is the bigger one:
    //
    //  1. Where the etch removed copper there is no copper — there is an
    //     18 um recess, and what fills it is chamber gas: on F.Cu the
    //     amplification gap continues down into the inter-pad grooves, on
    //     B.Cu the back gas gap continues up into them. The old model made
    //     those grooves out of thin copper.
    //  2. The copper is nowhere near uniformly distributed, and a single
    //     board-wide (or active-area-wide) coverage smears it radially. On
    //     F.Cu the pad field is ~98 % copper but the fan-out beyond it is
    //     18 % and the inner margin is 12 %; a 0.983 sheet over the whole
    //     board carries 25 % more copper than the artwork has, most of it
    //     dumped where the artwork has almost none.
    //
    // Coverage per radial band, and the pad grid itself, are measured
    // feature-by-feature from the production gerbers — include/P2PadMap.hh,
    // regenerated by scripts/gerber/extract_readout_pattern.py.
    constexpr double kBandArcMargin = 0.020;   // mm; > the 9.4 um arc sagitta
    auto BuildCuLayer = [&](const std::string& tag, G4double t, bool isFront) {
        using namespace P2;
        const bool pattern = isFront && fConfig.p2_patterned_readout;
        const std::string gapName = tag + "_Gap";

        G4LogicalVolume* envLV = nullptr;
        PlaceLayer(gapName, polyBoard, t, matGas, visGasV, &envLV);

        // Homogenized bands. With the pattern on, the bands covering the pad
        // field are skipped — real pads go there instead — so only the inner
        // margin and the fan-out are built from the table.
        for (int i = 0; i < kNCuBands; ++i) {
            const CuBand& b = kCuBands[i];
            const double f = isFront ? b.fCu : b.bCu;
            if (f <= 0.0) continue;
            if (pattern && b.rOut > kPadFieldRIn + 1e-6
                        && b.rIn  < kPadFieldROut - 1e-6) continue;

            // A band profile is a POLYGON: its 96-segment arc inscribes the
            // true circle and so dips ~9 um inside it (r*(1-cos(dphi/2nArc))).
            // Bands abutting each other are unaffected — they share vertices
            // — but a band abutting the pad field meets the outermost ring's
            // exact G4Tubs circle and overlaps it. Back that one boundary off
            // by more than the sagitta; both sides are the envelope's gas, so
            // the 20 um sliver is gas either way.
            const double rIn = (pattern && b.rIn >= kPadFieldROut - 1e-6)
                             ? b.rIn + kBandArcMargin : b.rIn;
            // WedgeOutline requires topCut < rOut, but most bands end well
            // inside the y = 540.01 chord and are never clipped by it. Clamp
            // to the y at which the outer arc meets the 60 deg offset edge:
            // that is exactly the "not clipped" case, and it keeps the shared
            // boundary between consecutive bands identical.
            const double aOutEdge = kSectorDeg * M_PI / 180.0
                                  + std::asin(kBoardEdgeOffset / b.rOut);
            const auto poly = WedgeOutline(
                rIn, b.rOut, kBoardEdgeOffset,
                std::min(kBoardTopCut, b.rOut * std::sin(aOutEdge)));
            const std::string nm = tag + "_Band" + std::to_string(i);
            auto* lv = new G4LogicalVolume(Prism(nm, poly, t), EffCu(f), nm);
            lv->SetVisAttributes(visCu);
            new G4PVPlacement(nullptr, G4ThreeVector(), lv, tag, envLV,
                              false, i, false);
        }
        if (!pattern) return;

        // The real pad field: 42 rings, each an annular envelope that ONE
        // G4PVReplica divides into its n pads. That is 3 logical volumes per
        // ring — 126 for all 1280 pads — and the navigator locates a pad by
        // index arithmetic instead of searching a daughter list, so the cost
        // is flat in the number of pads.
        //
        // G4 kPhi replica convention (G4ReplicaNavigation::ComputeTransformation):
        // cell i is rotated by -(offset + (i+0.5)*width), so it occupies
        // mother-phi [offset + i*w, offset + (i+1)*w] while its own solid is
        // defined SYMMETRIC about phi = 0. Hence the envelope starts half a
        // pitch below the first pad centre, that same value is the replica
        // offset, and the cell and pad solids are centred on zero.
        for (int k = 0; k < kNPadRings; ++k) {
            const PadRing& R = kPadRings[k];
            const G4double w    = R.dPhiPitch * rad;
            const G4double sPhi = R.phiFirst * rad - w / 2;
            const G4double rIn  = R.rIn * mm;
            const G4double rOut = R.rOut * mm;
            const std::string ks = std::to_string(k);

            auto* ringLV = new G4LogicalVolume(
                new G4Tubs(tag + "_Ring" + ks, rIn, rOut, t/2, sPhi, R.n * w),
                matGas, tag + "_Ring" + ks);
            ringLV->SetVisAttributes(G4VisAttributes::GetInvisible());
            new G4PVPlacement(nullptr, G4ThreeVector(), ringLV, gapName,
                              envLV, false, k, false);

            auto* cellLV = new G4LogicalVolume(
                new G4Tubs(tag + "_Cell" + ks, rIn, rOut, t/2, -w/2, w),
                matGas, tag + "_Cell" + ks);
            cellLV->SetVisAttributes(G4VisAttributes::GetInvisible());
            new G4PVReplica(gapName, cellLV, ringLV, kPhi, R.n, w, sPhi);

            // The pad spans the ring's full radial height; only the azimuthal
            // gap is left as gas. The 0.127 mm radial gaps between rings sit
            // outside the ring envelopes, in the layer envelope's gas.
            auto* padLV = new G4LogicalVolume(
                new G4Tubs(tag + "_Pad" + ks, rIn, rOut, t/2,
                           -R.dPhiPad * rad / 2, R.dPhiPad * rad),
                matCu, tag + "_Pad" + ks);
            padLV->SetVisAttributes(visCu);
            new G4PVPlacement(nullptr, G4ThreeVector(), padLV, tag,
                              cellLV, false, 0, false);
        }
    };

    BuildCuLayer("PCB_Cu_F", tCuF, true);
    PlaceLayer("PCB_FR4",           polyBoard,   tFR4,     matFR4,   visFR4);
    BuildCuLayer("PCB_Cu_B", tCuB, false);
    PlaceLayer("BackGas",           polyZone,    backGap,  matGas,   visGasV);

    // ── Gas frame rings (front: window plane -> Dynamask; back: mirror) ──
    // Front frame: plastic, confirmed 2026-08-05 (polycarbonate assumed for
    // the exact type), two rings either side of the drift-foil ledge,
    // standing on the bulk Dynamask. Back ring: carbon frame glued to the
    // PCB, back mylar glued to its rear face (confirmed, photos 2026-08-05);
    // it covers the bare board from the board edge in to the pad zone
    // (Alexandra 2026-10-06).
    //
    // Not modelled: the ~5 mm Dynamask border between the pad zone and the
    // drift opening, inside the amp layer -- it stays gas there. It lies
    // outside every pad, and the Saclay pillars at r_in reach 0.12 mm past
    // the pad-zone arc, so a border solid would need per-pillar clipping.
    {
        auto* dynSolid = new G4SubtractionSolid("BulkDynamask",
            Prism("BulkDynamask_o", polyBoard, tAmp),
            Prism("BulkDynamask_i", polyDriftOpening, tAmp + 2.0*mm));
        auto* dlv = new G4LogicalVolume(dynSolid, matDynamask, "BulkDynamask");
        dlv->SetVisAttributes(visFR4);
        new G4PVPlacement(nullptr,
            G4ThreeVector(0, 0, frameHWin + frameHBoard + tAmp/2), dlv,
            "BulkDynamask", worldLV, false, 0, true);

        auto* ringSolid = new G4SubtractionSolid("GasFrame",
            Prism("GasFrame_o", polyFrame, frameHWin),
            Prism("GasFrame_i", polyOpening, frameHWin + 2.0*mm));
        auto* lv = new G4LogicalVolume(ringSolid, matPlast, "GasFrame");
        lv->SetVisAttributes(visFrame);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, frameHWin/2), lv,
                          "GasFrame", worldLV, false, 0, true);

        auto* ledgeSolid = new G4SubtractionSolid("GasFrameLedge",
            Prism("GasFrameLedge_o", polyFrame, frameHBoard),
            Prism("GasFrameLedge_i", polyDriftOpening, frameHBoard + 2.0*mm));
        auto* llv = new G4LogicalVolume(ledgeSolid, matPlast, "GasFrameLedge");
        llv->SetVisAttributes(visFrame);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, frameHWin + frameHBoard/2),
                          llv, "GasFrameLedge", worldLV, false, 0, true);

        auto* backRing = new G4SubtractionSolid("GasFrameBack",
            Prism("GasFrameBack_o", polyBoard, backGap),
            Prism("GasFrameBack_i", polyZone, backGap + 2.0*mm));
        auto* blv = new G4LogicalVolume(backRing, matCarbon, "GasFrameBack");
        blv->SetVisAttributes(visFrame);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPCBEnd + backGap/2),
                          blv, "GasFrameBack", worldLV, false, 0, true);
    }

    // ── Bulged windows ────────────────────────────────────────────────────
    BuildWindow("FrontWindow", 0.0,      -1, hFront, polyOpening, openC);
    BuildWindow("BackWindow",  zBackWin, +1, hBack,  polyZone,    zoneC);

    // ── Fine-cut region: the gas plus the thin layers bounding it ─────────
    // Photon conversions in the mesh and the pad copper only matter if the
    // electron escapes into the gas, which is decided within a few microns
    // of the surface. The global 10 um e- range cut is ~80 keV in copper, so
    // without a finer cut here the escape physics is coarse exactly where it
    // is load-bearing. PhysicsList::SetCuts() attaches 1 um cuts to this
    // region. See docs/research/PHOTON_DISCRIMINATION_NOTES.md §6.1.
    {
        // "PCB_Cu_F_Gap" is the F.Cu layer envelope; a region root propagates
        // to every daughter, so naming it here covers the homogenized bands,
        // the 42 ring envelopes, the replica cells and all 1280 pads without
        // enumerating them.
        auto* fine = new G4Region("P2FineCut");
        for (const char* n : {"DriftGas", "AmpGas", "Micromesh",
                              "PCB_Cu_F_Gap", "DriftCathode_Al",
                              "DriftCathode_Mylar", "DriftCathode_Mylar2"}) {
            if (auto* lv = G4LogicalVolumeStore::GetInstance()->GetVolume(n, false))
                fine->AddRootLogicalVolume(lv);
        }
    }

    G4cout << "\n=== P2 wedge geometry ===" << G4endl;
    G4cout << "  Gas: " << fConfig.gas
           << "  (rho=" << matGas->GetDensity()/(mg/cm3) << " mg/cm3)" << G4endl;
    G4cout << "  Drift gap      : " << tDrift/mm  << " mm" << G4endl;
    G4cout << "  Amp gap        : " << tAmp/um    << " um" << G4endl;
    G4cout << "  Front/back gap : " << frontGap/mm << " / " << backGap/mm << " mm" << G4endl;
    if (twoFoils)
        G4cout << "  Drift cathode  : 2 x " << tCathMy/um << " um mylar, "
               << cathGap/mm << " mm apart, Al on drift side" << G4endl;
    else
        G4cout << "  Drift cathode  : " << tCathMy/um << " um mylar + "
               << tCathAl/um << " um Al, on the frame ledge" << G4endl;
    G4cout << "  Frame          : " << frameHWin/mm << " + " << frameHBoard/mm
           << " mm (window side + board side of the ledge; STEP: 4.0 + 4.0),"
           << " on " << tAmp/um << " um Dynamask" << G4endl;
    G4cout << "  Back frame     : carbon, " << backGap/mm
           << " mm, board outline -> pad zone" << G4endl;
    G4cout << "  Mesh           : woven SS "
           << fConfig.p2_mesh_open_um << "/" << fConfig.p2_mesh_wire_um
           << " um -> " << tMesh/um << " um slab, fill "
           << meshFill << G4endl;
    G4cout << "  Window bulge   : " << hFront/mm  << " / " << hBack/mm << " mm sag" << G4endl;
    G4cout << "  Pillars        : " << fConfig.p2_pillars << " mask, "
           << nPillars << " placed in the amp gap" << G4endl;
    G4cout << "  PCB            : " << tCuF/um << "um F.Cu / " << tFR4/um
           << "um FR4 / " << tCuB/um << "um B.Cu" << G4endl;
    G4cout << "  Readout copper : "
           << (fConfig.p2_patterned_readout
                 ? "REAL pads — " + std::to_string(P2::kNPads) + " annular "
                   "sectors in " + std::to_string(P2::kNPadRings) + " replicated rings"
                 : std::string("homogenized"))
           << ", " << P2::kNCuBands << " radial coverage bands" << G4endl;
    G4cout << "  z span         : " << fP2FrontZ/mm << " .. " << zMax/mm
           << " mm (front window plane = 0)" << G4endl;

    return worldPV;
}

// ============================================================
void DetectorConstruction::ConstructSDandField() {
    if (fDriftGasLV) {
        auto* sd = new SensitiveDetector("DriftGasSD", "DriftGasHits", "DriftGas", fConfig);
        G4SDManager::GetSDMpointer()->AddNewDetector(sd);
        SetSensitiveDetector(fDriftGasLV, sd);
        fDriftGasLV->SetUserLimits(new G4UserLimits(100*um));
    }
    if (fAmpGasLV) {
        auto* sd = new SensitiveDetector("AmpGasSD", "AmpGasHits", "AmpGas", fConfig);
        G4SDManager::GetSDMpointer()->AddNewDetector(sd);
        SetSensitiveDetector(fAmpGasLV, sd);
        fAmpGasLV->SetUserLimits(new G4UserLimits(100*um));
    }
    if (fHe3GasLV) {
        fHe3GasLV->SetUserLimits(new G4UserLimits(1.0*mm));
    }
}
