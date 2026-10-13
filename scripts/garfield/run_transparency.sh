#!/usr/bin/env bash
# run_transparency.sh <mode> <cell> <gasfile|""> <out> [n] [seed]
#
# Builds (if needed) and runs scripts/garfield/mesh_transparency.cc with the pinned Garfield build
# (P2_GARFIELD_INSTALL, default the lxplus work-area install).
#
#   bash scripts/garfield/run_transparency.sh transp \
#     /eos/user/a/akallits/p2_sim/comsol/p2_sps/h4mm \
#     /eos/user/a/akallits/p2_sim/gas_tables/ArIso.gas \
#     /eos/user/a/akallits/p2_sim/transparency/h4mm_ArIso 2000
set -e
export P2_GARFIELD_INSTALL="${P2_GARFIELD_INSTALL:-/afs/cern.ch/work/a/akallits/garfield/install}"
here="$(cd "$(dirname "$0")/../.." && pwd)"
# Under HTCondor the script is copied to the job scratch: use the clone.
[ -f "$here/scripts/setup_garfield.sh" ] || here="${P2_REPO:-$HOME/P2_Geant_Garfield}"
source "$here/scripts/setup_garfield.sh"
mkdir -p "$(dirname "$4")"
src="$here/scripts/garfield/mesh_transparency.cc"
exe="$P2_GARFIELD_INSTALL/../bin/mesh_transparency"
if [ ! -x "$exe" ] || [ "$src" -nt "$exe" ]; then
    mkdir -p "$(dirname "$exe")"
    g++ -O2 -std=c++20 -o "$exe" "$src" -I"$P2_GARFIELD_INSTALL/include" \
        $(root-config --cflags --libs) -L"$P2_GARFIELD_INSTALL/lib64" -lGarfield -lGarfieldRandom \
        -Wl,-rpath,"$P2_GARFIELD_INSTALL/lib64"
fi
"$exe" "$1" "$2" "$3" "$4" "${5:-2000}" "${6:-1}"
