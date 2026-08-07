// RunMeta.cc

#include "RunMeta.hh"
#include "GasMixtures.hh"

#include <cstdint>
#include <iomanip>
#include <sstream>

#ifndef P2_GIT_HASH
#define P2_GIT_HASH "unknown"
#endif
#ifndef P2_GIT_DIRTY
#define P2_GIT_DIRTY 0
#endif

namespace meta {

const char* GitHash() { return P2_GIT_HASH; }
bool        GitDirty() { return P2_GIT_DIRTY != 0; }

std::string GeometryDigest(const SimConfig& cfg) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(6);

    // Mode and gas define what was built at all.
    os << "mode=" << static_cast<int>(cfg.mode)
       << ";gas=" << cfg.gas;

    // Composition, not just the name: if a mixture's definition is ever
    // edited, runs made before and after must not hash the same.
    const gas::Mixture& mix = gas::Find(cfg.gas);
    for (const auto& p : mix.parts)
        os << ';' << p.component << '=' << p.volFrac;
    os << ";W=" << gas::MixtureWValue(mix, cfg.w_cf4_eV);

    // P2 stack. Every field DetectorConstruction reads when placing a volume.
    os << ";drift="    << cfg.p2_drift_mm
       << ";amp="      << cfg.p2_amp_um
       << ";meshWire=" << cfg.p2_mesh_wire_um
       << ";meshOpen=" << cfg.p2_mesh_open_um
       << ";frontGap=" << cfg.p2_front_gap_mm
       << ";cathGap="  << cfg.p2_cath_gap_mm
       << ";backGap="  << cfg.p2_back_gap_mm
       << ";bulgeF="   << cfg.p2_bulge_front_mm
       << ";bulgeB="   << cfg.p2_bulge_back_mm
       << ";window="   << cfg.p2_window_um
       << ";cathMy="   << cfg.p2_cath_mylar_um
       << ";cathAl="   << cfg.p2_cath_al_um
       << ";fCu="      << cfg.p2_fcu_coverage
       << ";bCu="      << cfg.p2_bcu_coverage;

    // Beam configuration is not geometry, but a pencil-beam run and a
    // beam-spread run are not interchangeable for any pad-level observable,
    // so they must not share a hash.
    os << ";gunX="   << cfg.p2_gun_x_mm
       << ";gunY="   << cfg.p2_gun_y_mm
       << ";theta="  << cfg.p2_gun_theta_deg
       << ";phi="    << cfg.p2_gun_phi_deg
       << ";spread=" << cfg.p2_beam_spread_mm;

    // Non-P2 modes.
    os << ";al="   << cfg.alThickness_mm
       << ";cfrp=" << cfg.cfrpThickness_mm;

    return os.str();
}

std::string GeometryHash(const SimConfig& cfg) {
    const std::string s = GeometryDigest(cfg);
    std::uint64_t h = 1469598103934665603ULL;          // FNV-1a offset basis
    for (unsigned char c : s) {
        h ^= c;
        h *= 1099511628211ULL;                          // FNV prime
    }
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << h;
    return os.str();
}

}  // namespace meta
