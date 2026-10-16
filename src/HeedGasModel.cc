// HeedGasModel.cc — see include/HeedGasModel.hh for what and why.

#include "HeedGasModel.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "G4AffineTransform.hh"
#include "G4EventManager.hh"
#include "G4FastStep.hh"
#include "G4FastTrack.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VSolid.hh"

#include "Garfield/ComponentConstant.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"

#include "EventAction.hh"
#include "EventData.hh"
#include "GasMixtures.hh"

namespace {
// Garfield's names for our gas components (MediumGas accepts these).
std::string GarfieldGasName(const std::string& c) {
    if (c == "iC4H10") return "iC4H10";
    if (c == "C2H6")   return "C2H6";
    return c;                       // Ar, Ne, CO2, CF4, He, ...
}
}  // namespace

HeedGasModel::HeedGasModel(const G4String& name, G4Region* envelope,
                           const SimConfig& cfg, const G4VSolid* gasSolid)
    : G4VFastSimulationModel(name, envelope), fConfig(cfg), fSolid(gasSolid) {
    // Heed reads its photo-absorption data from $HEED_DATABASE; default to
    // the database of the Garfield++ we were built against (batch jobs).
#ifdef P2_HEED_DATABASE
    setenv("HEED_DATABASE", P2_HEED_DATABASE, 0);
#endif
    // ── Gas: the mm_sim mixture, 20 C, 1 atm (the Geant4 density is the
    //    ideal gas at the same point; see GasMixtures.cc) ───────────────────
    const gas::Mixture& mix = gas::Find(cfg.gas);
    if (mix.parts.empty() || mix.parts.size() > 6)
        throw std::runtime_error("HeedGasModel: gas " + cfg.gas + " has " +
                                 std::to_string(mix.parts.size()) + " components");
    std::string n[6];
    double f[6] = {0, 0, 0, 0, 0, 0};
    for (std::size_t i = 0; i < mix.parts.size(); ++i) {
        n[i] = GarfieldGasName(mix.parts[i].component);
        f[i] = 100.0 * mix.parts[i].volFrac;
    }
    fGas = std::make_unique<Garfield::MediumMagboltz>();
    fGas->SetComposition(n[0], f[0], n[1], f[1], n[2], f[2], n[3], f[3],
                         n[4], f[4], n[5], f[5]);
    fGas->SetTemperature(293.15);
    fGas->SetPressure(760.0);

    // ── Geometry: Heed sees a box of +-kHalfXY around the track's entry
    //    point in x and y (DoIt shifts every track to x = y = 0) and the gas
    //    volume's full local z range. Handing it the whole P2 drift gas (a
    //    55 x 52 cm wedge, 4 mm thick, off-centre in its frame) made Heed
    //    return no clusters at all. The real outline is enforced in DoIt:
    //    electrons outside the solid are dropped. The field is irrelevant
    //    for the ionization. ────────────────────────────────────────────────
    G4ThreeVector lo, hi;
    fSolid->BoundingLimits(lo, hi);
    // 1 um margin in z: a track handed over on the volume's surface is, after
    // the global->local transform, often a rounding error outside it, and
    // Heed silently makes no clusters for a track that starts outside its
    // gas. Clusters outside the real solid are dropped in DoIt anyway.
    fZlo = lo.z() - 1.0 * um;
    fZhi = hi.z() + 1.0 * um;
    const double h = kHalfXY / cm;
    fCmp = std::make_unique<Garfield::ComponentConstant>();
    fCmp->SetElectricField(0., 0., 100.);
    fCmp->SetArea(-h, -h, fZlo / cm, h, h, fZhi / cm);
    fCmp->SetMedium(fGas.get());
    fSensor = std::make_unique<Garfield::Sensor>();
    fSensor->AddComponent(fCmp.get());
    fSensor->SetArea(-h, -h, fZlo / cm, h, h, fZhi / cm);

    fHeed = std::make_unique<Garfield::TrackHeed>(fSensor.get());
    fHeed->EnableDeltaElectronTransport();

    // Heed draws from Garfield's own engine: seed it from the run seed so a
    // job is reproducible.
#ifdef GARFIELD_OLD_RANDOM
    Garfield::randomEngine.Seed(static_cast<unsigned int>(cfg.seed));
#else
    Garfield::Random::SetEngine(Garfield::RandomEngineRoot(static_cast<unsigned int>(cfg.seed)));
#endif
}

HeedGasModel::~HeedGasModel() = default;

std::string HeedGasModel::HeedName(const std::string& g) {
    if (g == "e-" || g == "e+" || g == "mu-" || g == "mu+" ||
        g == "pi-" || g == "pi+")
        return g;
    if (g == "proton") return "p";
    return "";
}

// Lower limits as in the Garfield++ Geant4 interface example, except e+-:
// 1 MeV instead of 60 keV, so that ignoring multiple scattering over 4 mm of
// gas stays a few-degree effect.
double HeedGasModel::MinEnergyMeV(const std::string& g) const {
    if (g == "e-" || g == "e+")   return fConfig.heed_emin_e_MeV;
    if (g == "mu-" || g == "mu+") return 10.0;
    if (g == "pi-" || g == "pi+") return 20.0;
    if (g == "proton")            return 90.0;
    return 1e30;
}

G4bool HeedGasModel::IsApplicable(const G4ParticleDefinition& p) {
    return !HeedName(p.GetParticleName()).empty();
}

G4bool HeedGasModel::ModelTrigger(const G4FastTrack& ft) {
    const G4Track* t = ft.GetPrimaryTrack();
    if (t->GetKineticEnergy() / MeV < MinEnergyMeV(t->GetDefinition()->GetParticleName()))
        return false;
    // A track already sitting on the exit surface (just put there by DoIt)
    // must not be picked up again.
    const double L = fSolid->DistanceToOut(ft.GetPrimaryTrackLocalPosition(),
                                           ft.GetPrimaryTrackLocalDirection());
    return L > 1.0 * nm;
}

void HeedGasModel::DoIt(const G4FastTrack& ft, G4FastStep& fs) {
    const G4Track* track = ft.GetPrimaryTrack();
    const std::string g4name = track->GetDefinition()->GetParticleName();
    const double ekin = track->GetKineticEnergy();
    const double t0 = track->GetGlobalTime();
    const G4ThreeVector p0 = ft.GetPrimaryTrackLocalPosition();
    const G4ThreeVector d0 = ft.GetPrimaryTrackLocalDirection();
    const double L = fSolid->DistanceToOut(p0, d0);

    // TrackHeed rebuilds its cross-section tables when the particle or the
    // energy changes. The tables depend on beta-gamma only through slowly
    // varying terms, so the energy given to Heed is binned in 1 % steps.
    const std::string hname = HeedName(g4name);
    if (hname != fHeedParticle) {
        fHeed->SetParticle(hname);
        fHeedParticle = hname;
        fHeedEnergy = -1;
    }
    const double lstep = std::log(1.01);
    const double ebin = std::exp(std::round(std::log(ekin / eV) / lstep) * lstep);
    if (ebin != fHeedEnergy) {
        fHeed->SetKineticEnergy(ebin);
        fHeedEnergy = ebin;
    }
    // Heed frame = local frame shifted so the track starts at x = y = 0.
    const G4ThreeVector shift(p0.x(), p0.y(), 0.);
    const bool ok = fHeed->NewTrack(0., 0., p0.z() / cm, t0 / ns, d0.x(), d0.y(), d0.z());

    // ── Every Heed electron -> one cluster record ────────────────────────
    auto* ea = static_cast<EventAction*>(
        G4EventManager::GetEventManager()->GetUserEventAction());
    EventData& data = ea->GetEventData();
    const G4AffineTransform* toGlobal = ft.GetInverseAffineTransformation();
    const int trackID = track->GetTrackID();
    const G4VProcess* cproc = track->GetCreatorProcess();
    const G4LogicalVolume* blv = track->GetLogicalVolumeAtVertex();
    const TrackOrigin* origin = nullptr;
    if (auto it = data.trackOrigins.find(trackID); it != data.trackOrigins.end())
        origin = &it->second;

    double eloss = 0.0;   // [eV]
    int nElectrons = 0;
    for (const auto& c : fHeed->GetClusters()) {
        if (fSolid->Inside(G4ThreeVector(c.x * cm, c.y * cm, c.z * cm) + shift) == kOutside)
            continue;
        eloss += c.energy;
        const double eEach = c.electrons.empty() ? 0.0 : c.energy / c.electrons.size();
        for (const auto& e : c.electrons) {
            const G4ThreeVector local = G4ThreeVector(e.x * cm, e.y * cm, e.z * cm) + shift;
            if (fSolid->Inside(local) == kOutside) continue;
            const G4ThreeVector pos = toGlobal->TransformPoint(local);
            IonizationCluster cl;
            cl.x = pos.x() / mm;
            cl.y = pos.y() / mm;
            cl.z = pos.z() / mm;
            cl.edep = eEach;
            cl.nPrimary = 1;
            cl.trackID = trackID;
            cl.parentID = track->GetParentID();
            cl.volumeName = "DriftGas";
            cl.particleName = g4name;
            cl.kineticEnergy = ekin / MeV;
            cl.globalTime = e.t;
            cl.creatorProcess = cproc ? cproc->GetProcessName() : "primary";
            cl.birthVolume = blv ? blv->GetName() : "";
            if (origin) {
                cl.originProcess = origin->process;
                cl.originVolume = origin->volume;
                cl.originX = origin->x;
                cl.originY = origin->y;
                cl.originZ = origin->z;
                cl.originEnergy = origin->energy;
                cl.ancestorID = origin->ancestorID;
                cl.parentEnergy = origin->parentEnergy;
                cl.parentBirthVolume = origin->parentBirthVolume;
            }
            data.driftClusters.push_back(cl);
            ++nElectrons;
        }
    }
    if (!ok && ++fFailures <= 10)
        G4cerr << "HeedGasModel: TrackHeed::NewTrack failed at local z = " << p0.z() / mm
               << " mm (gas " << fZlo / mm << " .. " << fZhi / mm << ")" << G4endl;
    if (std::getenv("P2_HEED_DEBUG")) {
        G4cout << "HeedGasModel: " << g4name << " " << ekin / MeV << " MeV, local "
               << p0 / mm << " dir " << d0 << " L " << L / mm << " mm, Heed clusters "
               << fHeed->GetClusters().size() << ", electrons kept " << nElectrons
               << (ok ? "" : "  NEWTRACK FAILED") << G4endl;
    }
    data.nPrimaryDrift += nElectrons;
    if (trackID == 1) data.primaryInDrift = true;
    // edepDrift is added by SteppingAction from the deposit proposed below.

    // ── The particle continues from the exit point ───────────────────────
    const double dE = std::min(eloss * eV, ekin);
    fs.ProposeTotalEnergyDeposited(dE);
    if (ekin - dE <= 0.0) {
        fs.KillPrimaryTrack();
        return;
    }
    const double beta = track->GetVelocity() / c_light;
    fs.ProposePrimaryTrackFinalPosition(p0 + L * d0);   // local frame
    fs.ProposePrimaryTrackFinalMomentumDirection(d0);
    fs.ProposePrimaryTrackFinalKineticEnergy(ekin - dE);
    fs.ProposePrimaryTrackFinalTime(t0 + L / (beta * c_light));
    fs.ProposePrimaryTrackPathLength(L);
}
