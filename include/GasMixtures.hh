#pragma once
// GasMixtures.hh
//
// Single source of truth for the counting-gas mixtures.
//
// Why this file exists: composition used to live in DetectorConstruction
// (volume fractions -> G4Material) while the W-value lived in a *separate*
// hardcoded map in SteppingAction, keyed by the same string. Two places had
// to agree about what "NeIso" means, and nothing checked that they did. The
// gas-composition bug fixed on 2026-08-05 (research/GAS_FIX_NOTES.md) was
// exactly that class of error. Adding the seven campaign mixtures (P0.2)
// would have doubled the exposure, so composition, density and W are now
// derived from one table.
//
// Everything downstream reads this: DetectorConstruction builds the
// G4Material, SteppingAction takes W from the same record, and the run
// metadata writes the resolved composition into the ROOT file so a run can
// be identified without the submit script.

#include <string>
#include <vector>

namespace gas {

// ── Component properties ────────────────────────────────────────────────
//
// One record per pure gas that appears in any mixture.
//
//   molarMass  [g/mol]
//   nElectrons  molecular electron count, used to weight W by stopping power
//   wValue     [eV]  mean energy per ion pair, pure gas
struct Component {
    std::string name;
    double      molarMass;
    int         nElectrons;
    double      wValue;
};

// ── A mixture: components with VOLUME fractions ─────────────────────────
//
// Volume fraction is how gas systems are actually mixed and how every
// mixture in the literature is quoted. Mass fractions (what
// G4Material::AddMaterial wants) and density are derived, never typed in.
struct Fraction {
    std::string component;
    double      volFrac;
};

struct Mixture {
    std::string           name;
    std::string           label;      // human-readable, e.g. "Ne/iC4H10 90/10"
    std::vector<Fraction> parts;
    bool                  campaign;   // in the §7.4 campaign list
    std::string           note;       // caveats carried into the run banner
};

// ── Lookup ──────────────────────────────────────────────────────────────

// Throws std::runtime_error naming the available mixtures if not found.
const Mixture& Find(const std::string& name);
bool           Exists(const std::string& name);

const std::vector<Mixture>&  All();
const std::vector<Component>& Components();
const Component&              FindComponent(const std::string& name);

// Formatted list for --list-gases and for the "unknown gas" error.
std::string ListMixtures();

// ── Derived quantities ──────────────────────────────────────────────────
//
// Ideal gas at 20 C, 1 atm; partial densities add, so
//     rho_mix  = sum_i f_i * rho_i
//     f_mass_i = f_i * rho_i / rho_mix
double ComponentDensity(const std::string& component);   // g/cm3
double MixtureDensity(const Mixture& mix);               // g/cm3
double MassFraction(const Mixture& mix, std::size_t i);

// Mixture W-value [eV], stopping-power weighted:
//
//     W_mix = sum_i (f_i Z_i) / sum_i (f_i Z_i / W_i)
//
// The weight is the fraction of energy the mixture deposits in each species,
// which at a common pressure scales with the molecular electron count Z_i,
// not with the volume fraction alone. See the header comment in
// GasMixtures.cc for why this replaces the plain harmonic volume weighting
// named in SIM_CAMPAIGN_PLAN P0.2, and for the (larger) Penning caveat.
//
// wCF4 overrides the CF4 component W-value so the documented 35-52 eV spread
// can be bracketed as a systematic (P0.2); pass <= 0 for the table default.
double MixtureWValue(const Mixture& mix, double wCF4 = -1.0);
double MixtureWValue(const std::string& name, double wCF4 = -1.0);

}  // namespace gas
