// GasMixtures.cc
//
// ── On the W-values in this table ───────────────────────────────────────
//
// W (mean energy per ion pair) enters only `nPrimary = floor(edep/W)` in
// SteppingAction. OUTPUT_FORMAT.md already warns that this gives the right
// mean but Bernoulli rather than Fano fluctuations, and that `edep` is the
// observable to trust until Stage B does the conversion properly. Two
// further caveats, both larger than the arithmetic below:
//
//  1. PENNING TRANSFER is not modelled. In neon mixtures an excited Ne atom
//     can ionize a quencher molecule whose ionization potential lies below
//     the Ne excitation energy, so the *effective* W is materially lower
//     than any weighted average of pure-gas values -- by tens of percent in
//     Ne/iC4H10. Every Ne mixture W below is therefore an UPPER BOUND, and
//     comparing Ar with Ne on nPrimary at fixed edep is biased in argon's
//     favour. This is why the campaign plan defers the real number to the
//     Magboltz/Garfield stage with an explicit Penning transfer rate (P0.10).
//
//  2. The pure-gas W values are inherited from the pre-existing table in
//     SteppingAction.cc and have NOT been re-sourced against PDG. Isobutane
//     in particular is carried here at 26.0 eV while the commonly quoted
//     value is nearer 23.4 eV. Left unchanged deliberately: silently moving
//     it would shift every previously produced ArIso number with no record.
//     Flagged as a to-do in SIM_CAMPAIGN_PLAN P0.2.
//
// ── Weighting ───────────────────────────────────────────────────────────
//
// P0.2 as written called for a "harmonic vol-weighted" mixture W, i.e.
// W = 1 / sum_i(f_i / W_i). That weights each species by how much of the
// VOLUME it occupies, but what actually matters is how much of the ENERGY it
// absorbs, which at a common pressure scales with the molecular electron
// count Z_i. The stopping-power weighted form used here,
//
//     W_mix = sum_i (f_i Z_i) / sum_i (f_i Z_i / W_i),
//
// differs from the volume form by ~7 % for Ne/iC4H10 90/10 (32.8 vs 35.0 eV)
// because isobutane carries 34 electrons against neon's 10, so 10 % by
// volume is 27 % of the stopping power. Both sit well inside the Penning
// uncertainty above, but there is no reason to use the less defensible one.

#include "GasMixtures.hh"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace gas {

namespace {

// name, molarMass [g/mol], molecular electrons, W [eV]
const std::vector<Component> kComponents = {
    {"Ar",      39.948,   18, 26.4},
    {"Ne",      20.1797,  10, 36.4},
    {"He",       4.002602, 2, 41.3},
    {"CO2",     44.0095,  22, 33.0},
    {"CF4",     88.0043,  42, 34.0},   // 35-52 eV spread; override with --w-cf4
    {"CH4",     16.0425,  10, 27.3},
    {"C2H6",    30.069,   18, 26.0},
    {"iC4H10",  58.1222,  34, 26.0},   // see caveat 2 above
    {"N2",      28.0134,  14, 34.8},
};

// ── Mixtures ────────────────────────────────────────────────────────────
//
// The seven campaign gases (SIM_CAMPAIGN_PLAN §7.4, task P0.2) plus the
// mixtures inherited from MX17, which stay so the older modes and any
// previously produced sample can still be reproduced.
//
// Naming: where a pair of gases appears at more than one ratio the ratio is
// part of the name. `NeIso` (95/5) is inherited and is NOT one of the
// campaign points -- the campaign uses NeIso9010 and NeIso8020. Nothing
// resolves a bare "NeIso" to a campaign gas by accident.
const std::vector<Mixture> kMixtures = {
    // ── campaign gases ──────────────────────────────────────────────────
    {"ArIso", "Ar/iC4H10 95/5",
     {{"Ar", 0.95}, {"iC4H10", 0.05}}, true,
     "argon co-baseline"},

    {"ArCO2Iso9352", "Ar/CO2/iC4H10 93/5/2",
     {{"Ar", 0.93}, {"CO2", 0.05}, {"iC4H10", 0.02}}, true,
     "non-flammable argon ternary"},

    {"NeIso9010", "Ne/iC4H10 90/10",
     {{"Ne", 0.90}, {"iC4H10", 0.10}}, true,
     "W is an upper bound: Penning transfer not modelled"},

    {"NeIso8020", "Ne/iC4H10 80/20",
     {{"Ne", 0.80}, {"iC4H10", 0.20}}, true,
     "W is an upper bound: Penning transfer not modelled; 20% iC4H10 is above "
     "the flammability limit, see GAS_NOTES.md"},

    {"NeCO2Iso9532", "Ne/CO2/iC4H10 95/3/2",
     {{"Ne", 0.95}, {"CO2", 0.03}, {"iC4H10", 0.02}}, true,
     "ATLAS NSW mixture; worst case on timing (sigma_t ~25 ns) and on the "
     "zero-cluster floor"},

    {"NeIsoCO29055", "Ne/iC4H10/CO2 90/5/5",
     {{"Ne", 0.90}, {"iC4H10", 0.05}, {"CO2", 0.05}}, true,
     "W is an upper bound: Penning transfer not modelled"},

    {"NeCH4937", "Ne/CH4 93/7",
     {{"Ne", 0.93}, {"CH4", 0.07}}, true,
     "assumed Penning rate; no measured 150 um bulk gain curve"},

    {"NeC2H69010", "Ne/C2H6 90/10",
     {{"Ne", 0.90}, {"C2H6", 0.10}}, true,
     "assumed Penning rate; no measured 150 um bulk gain curve"},

    // ── inherited from MX17 ─────────────────────────────────────────────
    {"NeIso", "Ne/iC4H10 95/5",
     {{"Ne", 0.95}, {"iC4H10", 0.05}}, false,
     "inherited, not a campaign point (campaign uses NeIso9010/NeIso8020)"},

    {"ArCF4", "Ar/CF4 90/10",
     {{"Ar", 0.90}, {"CF4", 0.10}}, false, "inherited"},

    {"HeEth", "He/C2H6 96.5/3.5",
     {{"He", 0.965}, {"C2H6", 0.035}}, false, "inherited"},

    {"ArCO2", "Ar/CO2 70/30",
     {{"Ar", 0.70}, {"CO2", 0.30}}, false, "inherited"},

    {"NeCF4", "Ne/CF4 90/10",
     {{"Ne", 0.90}, {"CF4", 0.10}}, false, "inherited"},

    {"ArCF4Iso", "Ar/CF4/iC4H10 88/10/2",
     {{"Ar", 0.88}, {"CF4", 0.10}, {"iC4H10", 0.02}}, false, "inherited"},

    {"ArCF4CO2", "Ar/CF4/CO2 45/40/15",
     {{"Ar", 0.45}, {"CF4", 0.40}, {"CO2", 0.15}}, false, "inherited"},

    // ── pure gases ──────────────────────────────────────────────────────
    {"PureAr",     "Ar 100",      {{"Ar", 1.0}},     false, ""},
    {"PureNe",     "Ne 100",      {{"Ne", 1.0}},     false, ""},
    {"PureHe",     "He 100",      {{"He", 1.0}},     false, ""},
    {"PureCO2",    "CO2 100",     {{"CO2", 1.0}},    false, ""},
    {"PureCF4",    "CF4 100",     {{"CF4", 1.0}},    false, ""},
    {"PureCH4",    "CH4 100",     {{"CH4", 1.0}},    false, ""},
    {"PureEthane", "C2H6 100",    {{"C2H6", 1.0}},   false, ""},
    {"PureIso",    "iC4H10 100",  {{"iC4H10", 1.0}}, false, ""},
};

}  // namespace

// ============================================================
const std::vector<Component>& Components() { return kComponents; }
const std::vector<Mixture>&   All()        { return kMixtures; }

const Component& FindComponent(const std::string& name) {
    for (const auto& c : kComponents)
        if (c.name == name) return c;
    throw std::runtime_error("GasMixtures: unknown component '" + name + "'");
}

bool Exists(const std::string& name) {
    for (const auto& m : kMixtures)
        if (m.name == name) return true;
    return false;
}

const Mixture& Find(const std::string& name) {
    for (const auto& m : kMixtures)
        if (m.name == name) return m;
    throw std::runtime_error("Unknown gas '" + name + "'.\n" + ListMixtures());
}

std::string ListMixtures() {
    std::ostringstream os;
    os << "Available gases (* = campaign gas, SIM_CAMPAIGN_PLAN 7.4):\n";
    for (const auto& m : kMixtures) {
        os << "  " << (m.campaign ? "* " : "  ") << m.name;
        for (std::size_t i = m.name.size(); i < 16; ++i) os << ' ';
        os << m.label;
        if (!m.note.empty()) os << "   [" << m.note << "]";
        os << '\n';
    }
    return os.str();
}

// ============================================================
double ComponentDensity(const std::string& component) {
    // Ideal gas at the declared conditions, 20 C and 1 atm -- NOT 0 C STP.
    // Getting this wrong made every gas 7.3 % too dense until 2026-08-05.
    const double R = 8.314462618;          // J/(mol K)
    const double T = 293.15, P = 101325.0; // K, Pa
    return P * (FindComponent(component).molarMass * 1e-3) / (R * T) * 1e-3;
}

double MixtureDensity(const Mixture& mix) {
    double fsum = 0.0, rho = 0.0;
    for (const auto& p : mix.parts) {
        fsum += p.volFrac;
        rho  += p.volFrac * ComponentDensity(p.component);
    }
    if (std::abs(fsum - 1.0) > 1e-6)
        throw std::runtime_error("Gas " + mix.name +
                                 ": volume fractions must sum to 1");
    return rho;
}

double MassFraction(const Mixture& mix, std::size_t i) {
    return mix.parts[i].volFrac * ComponentDensity(mix.parts[i].component) /
           MixtureDensity(mix);
}

double MixtureWValue(const Mixture& mix, double wCF4) {
    double num = 0.0, den = 0.0;
    for (const auto& p : mix.parts) {
        const Component& c = FindComponent(p.component);
        const double w = (c.name == "CF4" && wCF4 > 0.0) ? wCF4 : c.wValue;
        const double weight = p.volFrac * c.nElectrons;
        num += weight;
        den += weight / w;
    }
    return num / den;
}

double MixtureWValue(const std::string& name, double wCF4) {
    return MixtureWValue(Find(name), wCF4);
}

}  // namespace gas
