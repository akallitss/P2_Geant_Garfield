#!/usr/bin/env bash
# run_gain.sh <cell> <gasfile> <out> <r1,r2,...> [n]
#
# Builds (if needed) and runs scripts/garfield/mesh_gain.cc with the pinned
# Garfield build (P2_GARFIELD_INSTALL, default the lxplus work-area install).
#
#   bash scripts/garfield/run_gain.sh \
#     /eos/user/a/akallits/p2_sim/comsol/p2_sps/h2mm_sweep_gap150/E625 \
#     /eos/user/a/akallits/p2_sim/gas_tables/ArCO2Iso9352.gas \
#     /eos/user/a/akallits/p2_sim/gain/test 0,0.2,0.4 10
set -e
export P2_GARFIELD_INSTALL="${P2_GARFIELD_INSTALL:-/afs/cern.ch/work/a/akallits/garfield/install}"
here="$(cd "$(dirname "$0")/../.." && pwd)"
# Under HTCondor the script is copied to the job scratch: use the clone.
[ -f "$here/scripts/setup_garfield.sh" ] || here="${P2_REPO:-$HOME/P2_Geant_Garfield}"
source "$here/scripts/setup_garfield.sh"
mkdir -p "$(dirname "$3")"
src="$here/scripts/garfield/mesh_gain.cc"
exe="$P2_GARFIELD_INSTALL/../bin/mesh_gain"
if [ ! -x "$exe" ] || [ "$src" -nt "$exe" ]; then
    mkdir -p "$(dirname "$exe")"
    # build to a private name and rename: parallel jobs must not race
    g++ -O2 -std=c++20 -o "$exe.$(hostname).$$" "$src" -I"$P2_GARFIELD_INSTALL/include" \
        $(root-config --cflags --libs) -L"$P2_GARFIELD_INSTALL/lib64" -lGarfield -lGarfieldRandom \
        -Wl,-rpath,"$P2_GARFIELD_INSTALL/lib64"
    mv -f "$exe.$(hostname).$$" "$exe"
fi
"$exe" "$1" "$2" "$3" "$4" "${5:-20}"
