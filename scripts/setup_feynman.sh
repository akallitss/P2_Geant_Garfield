#!/usr/bin/env bash
# setup_feynman.sh — environment on Feynman (CEA/IRFU cluster, Rocky 8, SLURM).
#
#   source scripts/setup_feynman.sh
#
# Feynman mounts /cvmfs/sft.cern.ch (public, no CERN account needed). LCG_106
# has an el8 build with exactly the production stack of the lxplus runs:
# Geant4 11.2.1, ROOT 6.32.02 and Garfield++ d3bbb29 (what mm_sim --heed was
# linked against for rnd2026_heed). So mm_sim here is the production binary,
# not an approximation of it.
#
# The transparency/gain programs need the pinned Garfield (927e5c21 + our
# patch). Build it once against the same view:
#   source scripts/setup_feynman.sh && bash scripts/setup_garfield.sh --build
# (installs into $P2_GARFIELD_INSTALL, default $P2_WORK/garfield/install).
#
# Sets:
#   P2_WORK   your work area (default: $HOME/p2; home quota is only 6 GB,
#             so point it at a work/scratch allocation when you have one)
#   P2_DATA   data bundle root (default: the shared group copy
#             /feynman/work/projets/idedipdetecs/p2_data, else $P2_WORK/data),
#             same layout as /eos/user/a/akallits/p2_sim. Read-only.
#   P2_OUT    where your outputs go (default: $P2_WORK/out)
#   GARFIELD_INSTALL  LCG_106 Garfield, used by build.sh for mm_sim --heed

export P2_LCG_VIEW="${P2_LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el8-gcc11-opt}"
if [ ! -f "$P2_LCG_VIEW/setup.sh" ]; then
    echo "[setup_feynman] $P2_LCG_VIEW not found: is /cvmfs/sft.cern.ch mounted on this node?" >&2
    return 1 2>/dev/null || exit 1
fi
source "$P2_LCG_VIEW/setup.sh"

_p2_here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export P2_REPO="$_p2_here"
export P2_SETUP="$_p2_here/scripts/setup_feynman.sh"   # sourced by batch/local job wrappers
export P2_WORK="${P2_WORK:-$HOME/p2}"
# Shared, read-only reference data (transferred 2026-11-01); else your own copy.
_p2_shared=/feynman/work/projets/idedipdetecs/p2_data
if [ -z "$P2_DATA" ]; then
    if [ -r "$_p2_shared/gas_tables" ]; then P2_DATA="$_p2_shared"; else P2_DATA="$P2_WORK/data"; fi
fi
export P2_DATA
export P2_OUT="${P2_OUT:-$P2_WORK/out}"     # your outputs: never write into a shared P2_DATA
export P2_GARFIELD_SRC="${P2_GARFIELD_SRC:-$P2_WORK/garfield/src}"
export P2_GARFIELD_INSTALL="${P2_GARFIELD_INSTALL:-$P2_WORK/garfield/install}"

# mm_sim links the LCG_106 Garfield (the production one); the pinned build in
# P2_GARFIELD_INSTALL is only for scripts/garfield/ programs, which find it
# themselves.
export GARFIELD_INSTALL="/cvmfs/sft.cern.ch/lcg/releases/LCG_106/Garfield++/d3bbb299/x86_64-el8-gcc11-opt"
export HEED_DATABASE="$GARFIELD_INSTALL/share/Heed/database"

echo "[setup_feynman] LCG_106 el8: Geant4 $(geant4-config --version)  ROOT $(root-config --version)"
echo "[setup_feynman] P2_DATA=$P2_DATA (read)  P2_OUT=$P2_OUT (write)"
[ -d "$P2_GARFIELD_INSTALL/include/Garfield" ] || \
    echo "[setup_feynman] pinned Garfield not built yet: bash scripts/setup_garfield.sh --build"
