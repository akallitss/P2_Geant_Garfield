// PrimaryGeneratorAction.cc
// Vacuum mode          : pencil beam at z = -10 cm, along +z.
// Full / sr90 modes    : beam from detector stack front, along +z.
// kLSCalib / kBackScintCalib : bare electron gun at detector front, optionally
//                              sampling the Sr-90/Y-90 beta spectrum from a CSV file.

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <map>

PrimaryGeneratorAction::PrimaryGeneratorAction(const SimConfig& cfg,
                                               const DetectorConstruction* detCon)
    : G4VUserPrimaryGeneratorAction(), fConfig(cfg), fDetCon(detCon) {

    fGun = std::make_unique<G4ParticleGun>(1);

    static const std::map<std::string, std::string> particleMap = {
        {"gamma",    "gamma"},
        {"neutron",  "neutron"},
        {"electron", "e-"},
        {"positron", "e+"},
        {"proton",   "proton"},
        {"muon",     "mu-"},
        {"muon+",    "mu+"},
        {"pion",     "pi-"},
        {"alpha",    "alpha"},
        {"triton",   "triton"},
    };

    auto it = particleMap.find(cfg.particle);
    if (it == particleMap.end())
        throw std::runtime_error("Unknown particle: " + cfg.particle);

    G4ParticleDefinition* particle =
        G4ParticleTable::GetParticleTable()->FindParticle(it->second);
    if (!particle)
        throw std::runtime_error("G4 particle not found: " + it->second);

    fGun->SetParticleDefinition(particle);
    fGun->SetParticleEnergy(cfg.energy);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));

    // Gun position
    G4double gunX = 0.0, gunY = 0.0;
    G4double gunZ = -10.0 * cm;
    bool useDetZ = (cfg.mode == SimMode::kFullExperiment  ||
                    cfg.mode == SimMode::kSr90Calibration ||
                    cfg.mode == SimMode::kSr90NoMM        ||
                    cfg.mode == SimMode::kLSCalib          ||
                    cfg.mode == SimMode::kBackScintCalib);
    if (useDetZ && fDetCon)
        gunZ = fDetCon->GetHe3GasCenterZ();

    if (cfg.mode == SimMode::kP2Wedge) {
        // World x/y are gerber coordinates (apex = beam axis at 0,0). The beam
        // passes through the aim point at the drift mid-plane and is tilted by
        // theta from the wedge normal, so an angle scan re-illuminates the same
        // pads instead of walking across the plane (P0.3).
        const G4double th  = cfg.p2_gun_theta_deg * deg;
        const G4double ph  = cfg.p2_gun_phi_deg   * deg;
        const G4ThreeVector dir(std::sin(th)*std::cos(ph),
                                std::sin(th)*std::sin(ph),
                                std::cos(th));

        const G4double aimZ = fDetCon ? fDetCon->GetP2DriftCenterZ() : 0.0;
        const G4ThreeVector aim(cfg.p2_gun_x_mm * mm,
                                cfg.p2_gun_y_mm * mm,
                                aimZ);

        // Back-track along the beam to the gun. The standoff must clear the
        // front window bulge, which reaches GetP2FrontZ upstream of the window
        // plane; at large theta the beam also has to clear it laterally, so
        // require the gun to sit upstream of the bulge apex with margin.
        G4double standoff = cfg.p2_gun_standoff_mm * mm;
        if (fDetCon) {
            const G4double frontZ  = fDetCon->GetP2FrontZ();   // negative
            const G4double needed  = (aimZ - frontZ + 10.0*mm) / std::cos(th);
            if (standoff < needed) {
                G4cout << "PrimaryGeneratorAction: raising gun standoff "
                       << standoff/mm << " -> " << needed/mm
                       << " mm to clear the front window at theta = "
                       << cfg.p2_gun_theta_deg << " deg" << G4endl;
                standoff = needed;
            }
        }

        const G4ThreeVector pos = aim - standoff * dir;
        gunX = pos.x(); gunY = pos.y(); gunZ = pos.z();
        fGun->SetParticleMomentumDirection(dir);

        fP2Mode = true; fP2Aim = aim; fP2Dir = dir; fP2Standoff = standoff;

        if (cfg.p2_beam_spread_mm <= 0.0) {
            // Warn if the fixed aim point lands on a pad boundary. A pencil
            // beam parked on a boundary makes "which pad has the most
            // charge" a coin flip and inflates pad multiplicity, and it does
            // so invisibly -- MX17 lost a first result to exactly this
            // (RESPONSE_SIM_PLAN §7). P2's own historical default, r=355 mm,
            // sat 35 um from a ring boundary.
            const double r = std::hypot(cfg.p2_gun_x_mm, cfg.p2_gun_y_mm);
            const double ringPitch = 11.4290;                 // mm, measured
            const double rInner    = 120.714;                 // first ring centre
            const double frac = std::fabs(std::fmod(r - rInner + 0.5*ringPitch,
                                                    ringPitch) / ringPitch - 0.5);
            if (frac > 0.35)
                G4cout << "PrimaryGeneratorAction: WARNING - the aim point "
                       << "r = " << r << " mm sits " << 100.0*frac
                       << "% of a ring pitch from the nearest pad centre, "
                       << "i.e. near a RADIAL PAD BOUNDARY. Pad multiplicity "
                       << "and any positional observable will be biased. "
                       << "Use --beam-spread (>= one ring pitch, 11.43 mm) "
                       << "for those observables." << G4endl;
        }
    }

    fGun->SetParticlePosition(G4ThreeVector(gunX, gunY, gunZ));

    // Load Sr-90/Y-90 spectrum if requested
    if (!cfg.spectrum_file.empty()) {
        LoadSpectrum(cfg.spectrum_file);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
void PrimaryGeneratorAction::LoadSpectrum(const std::string& filepath) {
    std::ifstream f(filepath);
    if (!f.is_open())
        throw std::runtime_error("Cannot open spectrum file: " + filepath);

    std::vector<double> energies, weights;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        // Skip header lines (contain non-numeric first token)
        std::istringstream ss(line);
        double e, w;
        if (!(ss >> e >> w)) continue;
        if (e < 0 || w < 0) continue;
        energies.push_back(e);
        weights.push_back(w);
    }

    if (energies.size() < 2)
        throw std::runtime_error("Spectrum file too short: " + filepath);

    // Build normalised CDF
    double total = std::accumulate(weights.begin(), weights.end(), 0.0);
    fSpecEnergies = energies;
    fSpecCDF.resize(weights.size());
    double cumsum = 0.0;
    for (size_t i = 0; i < weights.size(); ++i) {
        cumsum += weights[i] / total;
        fSpecCDF[i] = cumsum;
    }
    fSpecCDF.back() = 1.0;  // ensure exact 1 at end
    fUseSpectrum = true;

    G4cout << "PrimaryGeneratorAction: Loaded spectrum from " << filepath
           << "  (" << fSpecEnergies.size() << " points, "
           << "E_max=" << fSpecEnergies.back() << " MeV)" << G4endl;
}

// ─────────────────────────────────────────────────────────────────────────────
double PrimaryGeneratorAction::SampleSpectrum() const {
    double r = G4UniformRand();
    // Inverse CDF via binary search
    auto it = std::lower_bound(fSpecCDF.begin(), fSpecCDF.end(), r);
    size_t i = std::distance(fSpecCDF.begin(), it);
    if (i >= fSpecEnergies.size()) i = fSpecEnergies.size() - 1;
    return fSpecEnergies[i];
}

// ─────────────────────────────────────────────────────────────────────────────
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    if (fUseSpectrum) {
        double E = SampleSpectrum();
        fGun->SetParticleEnergy(E * MeV);
    }

    // --beam-spread: scatter the impact point over a disc in the wedge plane
    // so pad-level observables average over the pad cell instead of being
    // read off one fixed point. Sampled per event, in the plane transverse
    // to the beam so it stays a disc at any --gun-theta.
    if (fP2Mode && fConfig.p2_beam_spread_mm > 0.0) {
        const G4double R = fConfig.p2_beam_spread_mm * mm;
        const G4double rho = R * std::sqrt(G4UniformRand());   // uniform in area
        const G4double psi = CLHEP::twopi * G4UniformRand();
        // Orthonormal basis transverse to the beam.
        G4ThreeVector u = fP2Dir.orthogonal().unit();
        G4ThreeVector v = fP2Dir.cross(u).unit();
        const G4ThreeVector aim = fP2Aim + rho*(std::cos(psi)*u + std::sin(psi)*v);
        fGun->SetParticlePosition(aim - fP2Standoff * fP2Dir);
    }
    fGun->GeneratePrimaryVertex(event);
}
