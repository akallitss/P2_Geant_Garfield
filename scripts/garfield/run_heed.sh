#!/usr/bin/env bash
# run_heed.sh <gas> <gap_mm> <theta_deg> <n> <out> [particle] [E_kin_eV] [seed]
# Builds (if needed) and runs scripts/garfield/heed_mip.cc with the pinned
# Garfield build (see run_transparency.sh).
set -e
export P2_GARFIELD_INSTALL="${P2_GARFIELD_INSTALL:-/afs/cern.ch/work/a/akallits/garfield/install}"
here="$(cd "$(dirname "$0")/../.." && pwd)"
[ -f "$here/scripts/setup_garfield.sh" ] || here="${P2_REPO:-$HOME/P2_Geant_Garfield}"
source "$here/scripts/setup_garfield.sh"
mkdir -p "$(dirname "$5")"
src="$here/scripts/garfield/heed_mip.cc"
exe="$P2_GARFIELD_INSTALL/../bin/heed_mip"
if [ ! -x "$exe" ] || [ "$src" -nt "$exe" ]; then
    mkdir -p "$(dirname "$exe")"
    g++ -O2 -std=c++20 -o "$exe" "$src" -I"$P2_GARFIELD_INSTALL/include" \
        $(root-config --cflags --libs) -L"$P2_GARFIELD_INSTALL/lib64" -lGarfield -lGarfieldRandom \
        -Wl,-rpath,"$P2_GARFIELD_INSTALL/lib64"
fi
"$exe" "$@"
