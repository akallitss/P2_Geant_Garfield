#pragma once
// RunMeta.hh
//
// Provenance that has to travel INSIDE the output file.
//
// SIM_CAMPAIGN_PLAN §9 requires one manifest row per run point carrying
// runID, mode, gas, gap, particle, E, theta, position, N, seed, code git-hash
// and a geometry-parameter hash, and says "no un-manifested runs". Until now
// none of that was in the ROOT file: the submit script was the only thing
// that knew, and the THROWN event count -- which every per-incident-photon
// probability divides by, and which differs from GetEntries() by ~300x under
// --skip-empty -- existed only in the stdout log. A file separated from its
// job log was unnormalizable and, strictly, unidentifiable.
//
// So each worker file now carries a one-row RunMeta tree. `hadd` concatenates
// them, which is the wanted behaviour: the merged file holds one row per
// worker, thrown counts sum, and a mismatch in git hash or geometry hash
// between rows is visible instead of silent.

#include "SimConfig.hh"
#include <string>

namespace meta {

// Source revision the binary was built from, "unknown" if git was
// unavailable at configure time. Dirty means uncommitted changes were
// present -- for a production run that is itself the finding.
const char* GitHash();
bool        GitDirty();

// Every geometry-affecting field of SimConfig, formatted as key=value. This
// is the text that gets hashed; it is also written to the file so a hash
// mismatch can be diagnosed instead of merely detected.
std::string GeometryDigest(const SimConfig& cfg);

// FNV-1a over GeometryDigest, hex. Two runs with the same hash were produced
// by the same geometry; different hashes must not be merged.
std::string GeometryHash(const SimConfig& cfg);

}  // namespace meta
