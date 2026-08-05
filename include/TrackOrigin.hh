#pragma once
// TrackOrigin.hh — provenance record carried down each track tree.
//
// The question this exists to answer: for every bit of ionization scored in
// the gas, *which interaction, in which layer, started the chain that
// produced it*.  One-level parentID cannot answer it — a delta ray from a
// photoelectron has creatorProcess == "eIoni" and a parent that is itself a
// secondary, so the link back to the photoelectric conversion is lost.
//
// The rule (implemented in TrackingAction):
//
//   * a primary track gets origin = "primary";
//   * a track whose parent is CHARGED inherits its parent's origin verbatim
//     (deltas and bremsstrahlung of an electron stay attributed to whatever
//     made that electron);
//   * a track whose parent is NEUTRAL starts a new origin, tagged with its
//     own creator process and the volume it was created in.
//
// So the origin is "the interaction at which a neutral ancestor became
// charged" — exactly the conversion we want to bin photon runs by — and for
// electron runs it stays "primary" all the way down the shower.
//
// See docs/research/PHOTON_DISCRIMINATION_NOTES.md §6.2.

#include <string>

struct TrackOrigin {
    // ── The origin interaction ────────────────────────────────────────────
    std::string process    = "unknown";  // phot, compt, conv, Rayl, primary, ...
    std::string volume     = "";         // logical volume it happened in
    double      x = 0, y = 0, z = 0;     // where [mm]
    double      energy     = 0.0;        // KE of the track born there [MeV]
    int         pdg        = 0;          // its PDG code
    int         ancestorID = 0;          // its trackID

    // ── The neutral parent that converted (photon runs: the X-ray) ────────
    // Lets the analysis separate a direct conversion of the incident photon
    // from re-absorption of a fluorescence photon born in the copper: for
    // Cu-K the parent was created at 8.05 keV inside PCB_Cu_F.
    double      parentEnergy      = 0.0; // parent's KE at ITS creation [MeV]
    std::string parentBirthVolume = "";  // where the parent was created

    // ── Bookkeeping for the inheritance rule ──────────────────────────────
    double      ownCharge = 0.0;         // charge of the track holding this
                                         // record (not of the origin)
};
