#!/usr/bin/env bash
# setup_local.sh — environment for running the P2 chain on your own machine
# (no CERN account, no CVMFS/EOS/HTCondor). The lxplus counterpart is
# setup_lxplus.sh. Full instructions: docs/SETUP.md.
#
#   conda activate p2sim                 # from environment.yml
#   source scripts/setup_local.sh
#
# Sets:
#   P2_DATA      root of the data bundle (default ~/p2_data). Same layout as
#                /eos/user/a/akallits/p2_sim on EOS, so every path in the docs
#                works with /eos/user/a/akallits/p2_sim -> $P2_DATA.
#   P2_REPO      this checkout.
#   GARFIELD_INSTALL, LD_LIBRARY_PATH, HEED_DATABASE  (via setup_garfield.sh,
#                if Garfield++ is built; build it once with
#                bash scripts/setup_garfield.sh --build)

_p2_here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export P2_REPO="$_p2_here"
export P2_SETUP="$_p2_here/scripts/setup_local.sh"     # sourced by local job wrappers
export P2_DATA="${P2_DATA:-$HOME/p2_data}"
export P2_OUT="${P2_OUT:-$HOME/p2_out}"     # your outputs (kept apart from the data bundle)

_p2_ok=1
for _t in geant4-config root-config cmake; do
    if ! command -v "$_t" >/dev/null 2>&1; then
        echo "[setup_local] $_t not found: activate the conda env first (conda activate p2sim)" >&2
        _p2_ok=0
    fi
done
if [ $_p2_ok = 0 ]; then return 1 2>/dev/null || exit 1; fi

# CMake must find the env's libraries before the system's (e.g. a system
# libexpat older than the one Geant4 was built against fails the configure).
if [ -n "$CONDA_PREFIX" ]; then
    case ":$CMAKE_PREFIX_PATH:" in
        *":$CONDA_PREFIX:"*) ;;
        *) export CMAKE_PREFIX_PATH="$CONDA_PREFIX${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}" ;;
    esac
fi

# conda's activation script sets the Geant4 data variables; check one.
if [ -z "$G4LEDATA" ] || [ ! -d "$G4LEDATA" ]; then
    echo "[setup_local] WARNING: G4LEDATA not set or missing; Geant4 data sets not found." >&2
    echo "[setup_local]          Re-activate the env, or run: geant4-config --datasets" >&2
fi

# Garfield++ (optional: needed for mm_sim --heed and the transparency/gain programs)
if [ -d "${P2_GARFIELD_INSTALL:-$HOME/garfield/install}/include/Garfield" ]; then
    export P2_GARFIELD_INSTALL="${P2_GARFIELD_INSTALL:-$HOME/garfield/install}"
    source "$_p2_here/scripts/setup_garfield.sh"
else
    echo "[setup_local] Garfield++ not built yet (needed for --heed, the production physics)."
    echo "[setup_local]   bash scripts/setup_garfield.sh --build      (~15 min, once)"
fi

[ -d "$P2_DATA" ] || echo "[setup_local] P2_DATA=$P2_DATA does not exist yet: see docs/SETUP.md §3."

echo "[setup_local] Geant4 $(geant4-config --version)  ROOT $(root-config --version)"
echo "[setup_local] P2_DATA=$P2_DATA"
echo "[setup_local] Next: bash scripts/build.sh"
