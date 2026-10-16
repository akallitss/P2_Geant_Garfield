// PhysicsList.cc
// Builds the physics for the Micromegas simulation.
// Key considerations:
//   - Gammas  : photoelectric, Compton, pair production (all via EM option4)
//   - Neutrons: elastic + inelastic via FTFP_BERT (QGSP or FTFP + Bertini cascade)
//              + thermal neutrons via NeutronHP
//   - Electrons: accurate low-energy EM (Livermore model in option4 goes down to eV)
//   - Ionization: G4ionIonisation tracks energy loss in gas steps

#include "PhysicsList.hh"

// Modular physics components
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4PAIModel.hh"
#include "G4ParticleTable.hh"
#include "G4FastSimulationPhysics.hh"

PhysicsList::PhysicsList(const SimConfig& cfg)
    : G4VModularPhysicsList(), fPAI(cfg.pai), fPixe(cfg.pixe), fHeed(cfg.heed) {
    SetVerboseLevel(0);

    // EM physics: option4 uses Livermore models below 100 keV for e-/gamma
    // This gives accurate photoelectric + Auger + ionization in low-Z gas
    RegisterPhysics(new G4EmStandardPhysics_option4(0));

    // EM extras: synchrotron, GDR etc (harmless to include)
    RegisterPhysics(new G4EmExtraPhysics(0));

    // Hadronic elastic with HP (high-precision neutron data, <20 MeV)
    RegisterPhysics(new G4HadronElasticPhysicsHP(0));

    // Hadronic inelastic with HP neutrons
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP(0));

    // Decay
    RegisterPhysics(new G4DecayPhysics(0));

    // Radioactive decay (useful for activation studies)
    RegisterPhysics(new G4RadioactiveDecayPhysics(0));

    // Step limiter (respects G4UserLimits set in detector volumes)
    RegisterPhysics(new G4StepLimiterPhysics());

    // Kill neutrons below 1 eV after tracking to avoid infinite loops
    RegisterPhysics(new G4NeutronTrackingCut(0));

    // --heed: the fast-simulation process that hands charged tracks in the
    // drift gas to HeedGasModel (DetectorConstruction::ConstructSDandField).
    if (fHeed) {
        auto* fsp = new G4FastSimulationPhysics();
        for (const char* p : {"e-", "e+", "mu-", "mu+", "pi+", "pi-", "proton"})
            fsp->ActivateFastSimulation(p);
        RegisterPhysics(fsp);
    }

    // ── Atomic deexcitation ───────────────────────────────────────────────
    // Set AFTER registering the EM constructor, whose own ctor writes to the
    // same G4EmParameters singleton.
    //
    // deexcitationIgnoreCut is the load-bearing one. Production cuts are
    // *range* cuts converted per material, and the 0.1 mm gamma cut below is
    // ~5-10 keV inside copper -- so without this flag the Cu-K fluorescence
    // line at 8.05 keV (yield 0.44) sits at or under threshold and is
    // silently not produced. That line is the dominant gas-sensitive photon
    // channel in the P2 stack: a soft photon leaving the pad plane and being
    // re-absorbed in the gas is ~10x more likely in argon than in neon, and
    // it is the fake class the charge cut rejects best. Losing it would bias
    // the whole Ar-vs-Ne comparison.
    // See docs/research/PHOTON_DISCRIMINATION_NOTES.md §4.3, §6.1.
    auto* emp = G4EmParameters::Instance();
    emp->SetFluo(true);
    emp->SetAuger(true);
    emp->SetPixe(fPixe);
    emp->SetDeexcitationIgnoreCut(true);

    // The documented route to PAI. In Geant4 11.3 it does nothing with
    // option4 (G4EmModelActivator skips processes whose isIonisation flag is
    // not yet set); kept as --pai params to test it on our 11.2.
    if (fPAI == "params") {
        for (const char* p : {"e-", "e+", "mu-", "mu+", "proton", "pi+", "pi-"})
            emp->AddPAIModel(p, "P2GasPAI", "pai");
    }
}

// ── PAI in the gas only (--pai config) ──────────────────────────────────
// Option4's ionisation uses the Urban fluctuation model, which puts the most
// probable energy loss in a 4 mm Ar layer far above Heed's (Dylan's MX17
// check, 2026-10-15). G4PAIModel is built from the gas photo-absorption
// cross section, as Heed is. It replaces the ionisation model of charged
// particles in the region P2GasPAI (DriftGas + AmpGas) only; the solids keep
// option4. Done here, after the processes exist, as in Geant4 TestEm8.
void PhysicsList::ConstructProcess() {
    G4VModularPhysicsList::ConstructProcess();
    if (fPAI != "config") return;
    auto* conf = G4LossTableManager::Instance()->EmConfigurator();
    const std::pair<const char*, const char*> procs[] = {
        {"e-", "eIoni"}, {"e+", "eIoni"}, {"mu-", "muIoni"}, {"mu+", "muIoni"},
        {"proton", "hIoni"}, {"pi+", "hIoni"}, {"pi-", "hIoni"}};
    // With --heed, tracks below the Heed thresholds still step through
    // P2HeedGas in Geant4: give them PAI there too.
    std::vector<const char*> regions = {"P2GasPAI"};
    if (fHeed) regions.push_back("P2HeedGas");
    for (const char* reg : regions) {
        for (const auto& [name, proc] : procs) {
            auto* p = G4ParticleTable::GetParticleTable()->FindParticle(name);
            auto* pai = new G4PAIModel(p, "PAIModel");
            conf->SetExtraEmModel(name, proc, pai, reg, 0.0, 100.0 * TeV, pai);
        }
    }
}

void PhysicsList::SetCuts() {
    // Global range cuts. These are *range* cuts: the same number means very
    // different energies in gas and in metal (10 um is ~1 keV in argon at
    // 1 atm but ~80 keV in copper), which is why the thin high-Z layers get
    // their own region below.
    SetCutValue(1.0 * mm,  "proton");
    SetCutValue(10.0 * um, "e-");       // short cut for electrons -- capture delta rays
    SetCutValue(10.0 * um, "e+");
    SetCutValue(0.1  * mm, "gamma");    // 0.1 mm photon production threshold

    // ── Fine cuts over the gas and the thin layers bounding it ────────────
    // The mesh and the copper pad plane are the two biggest photon-conversion
    // sources in the stack, and what matters is whether the conversion
    // electron *escapes into the gas*. That is decided within a few microns
    // of the surface, so the secondary-production threshold there has to be
    // well below the global one. The region is built in
    // DetectorConstruction::ConstructP2(); absent (other modes) we skip.
    if (auto* reg = G4RegionStore::GetInstance()->GetRegion("P2FineCut", false)) {
        auto* cuts = new G4ProductionCuts();
        cuts->SetProductionCut(1.0 * um, G4Gamma::Gamma());
        cuts->SetProductionCut(1.0 * um, G4Electron::Electron());
        cuts->SetProductionCut(1.0 * um, G4Positron::Positron());
        reg->SetProductionCuts(cuts);
    }

    // The PAI and Heed gas regions need a cuts object of their own: G4EmModelManager
    // assigns couples to regions by G4ProductionCuts POINTER, so a region
    // left on the default cuts would capture every couple in the world and
    // put PAI in all the solids too. Same 1 um as P2FineCut.
    for (const char* name : {"P2GasPAI", "P2HeedGas"}) {
        auto* reg = G4RegionStore::GetInstance()->GetRegion(name, false);
        if (!reg) continue;
        auto* cuts = new G4ProductionCuts();
        cuts->SetProductionCut(1.0 * um, G4Gamma::Gamma());
        cuts->SetProductionCut(1.0 * um, G4Electron::Electron());
        cuts->SetProductionCut(1.0 * um, G4Positron::Positron());
        reg->SetProductionCuts(cuts);
    }
}
