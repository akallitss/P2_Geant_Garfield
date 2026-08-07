// PrimaryGeneratorAction.cc
// Vacuum mode          : pencil beam at z = -10 cm, along +z.
// Full / sr90 modes    : beam from detector stack front, along +z.
// kLSCalib / kBackScintCalib : bare electron gun at detector front, optionally
//                              sampling the Sr-90/Y-90 beta spectrum from a CSV file.

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "P2PadMap.hh"

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

        if (cfg.p2_beam_spread_mm <= 0.0) WarnIfAimPointOffPad(cfg);
    }

    fGun->SetParticlePosition(G4ThreeVector(gunX, gunY, gunZ));

    // Load Sr-90/Y-90 spectrum if requested
    if (!cfg.spectrum_file.empty()) {
        LoadSpectrum(cfg.spectrum_file);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Is the fixed aim point actually on a pad, in BOTH coordinates?
//
// A pencil beam parked on a pad boundary makes "which pad has the most
// charge" a coin flip and inflates pad multiplicity, invisibly. MX17 lost a
// first result to exactly that (RESPONSE_SIM_PLAN §7).
//
// The pad plane is polar, so there are TWO ways to land on a boundary and
// fixing one does not fix the other. The historical default (r = 355 mm) sat
// 35 um from a *radial* boundary; moving it to a ring centre then put it at
// phi = 30.0000 deg, which is 61 um from the centre of an *azimuthal* gap --
// only 25 % of events deposited in pad copper at all, against 100 % at a
// true pad centre. Hence an exact 2-D check against the real artwork
// (include/P2PadMap.hh, extracted from the gerbers) rather than a radial
// heuristic against the mapping files.
void PrimaryGeneratorAction::WarnIfAimPointOffPad(const SimConfig& cfg) {
    const double x = cfg.p2_gun_x_mm, y = cfg.p2_gun_y_mm;
    const double r = std::hypot(x, y);
    double phi = std::atan2(y, x);

    const P2::PadRing* ring = nullptr;
    int ringIdx = -1;
    for (int i = 0; i < P2::kNPadRings; ++i) {
        if (r >= P2::kPadRings[i].rIn && r <= P2::kPadRings[i].rOut) {
            ring = &P2::kPadRings[i]; ringIdx = i; break;
        }
    }
    if (!ring) {
        G4cout << "PrimaryGeneratorAction: WARNING - the aim point r = " << r
               << " mm is NOT inside any pad ring (pad field spans "
               << P2::kPadFieldRIn << " .. " << P2::kPadFieldROut
               << " mm). The beam is not hitting instrumented area."
               << G4endl;
        return;
    }

    // Nearest pad centre in phi, and how far off it we are.
    const double fidx = (phi - ring->phiFirst) / ring->dPhiPitch;
    const double k    = std::round(fidx);
    const double dphi = std::fabs(phi - (ring->phiFirst + k * ring->dPhiPitch));
    const bool onCopperPhi = (dphi <= 0.5 * ring->dPhiPad);
    const bool inRing      = (k >= 0 && k < ring->n);

    // Radial margin to the ring's copper edge.
    const double rMid = 0.5 * (ring->rIn + ring->rOut);
    const double dr   = std::fabs(r - rMid);

    if (!inRing) {
        G4cout << "PrimaryGeneratorAction: WARNING - the aim point phi = "
               << phi/deg << " deg is outside ring " << ringIdx
               << "'s populated pads." << G4endl;
        return;
    }
    if (!onCopperPhi) {
        G4cout << "PrimaryGeneratorAction: WARNING - the aim point sits in an "
               << "AZIMUTHAL INTER-PAD GAP: " << dphi * r * 1000.0
               << " um from the nearest pad centre in phi (pad half-width "
               << 0.5 * ring->dPhiPad * r * 1000.0 << " um) on ring "
               << ringIdx << ". Pad multiplicity, charge sharing and every "
               << "positional observable will be biased. Use --beam-spread, "
               << "or aim at a pad centre." << G4endl;
        return;
    }
    if (dr > 0.35 * (ring->rOut - ring->rIn)) {
        G4cout << "PrimaryGeneratorAction: WARNING - the aim point is "
               << dr * 1000.0 << " um from the radial centre of ring "
               << ringIdx << " (half-height "
               << 0.5 * (ring->rOut - ring->rIn) * 1000.0
               << " um), i.e. close to a RADIAL pad boundary. Use "
               << "--beam-spread for pad-level observables." << G4endl;
        return;
    }

    if (cfg.verbose)
        G4cout << "PrimaryGeneratorAction: aim point is on pad copper, ring "
               << ringIdx << ", " << dr*1000.0 << " um from the radial centre "
               << "and " << dphi*r*1000.0 << " um from the azimuthal centre."
               << G4endl;
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
