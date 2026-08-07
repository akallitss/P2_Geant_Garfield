#!/usr/bin/env bash
# setup_garfield.sh
#
# Puts Garfield++ (and the ROOT it was built against) on the environment for
# the P2 campaign's Garfield-dependent work: Magboltz gas tables (P0.10),
# the mesh-transparency field solve (P0.13) and the induced-current model
# (P0.17).
#
#   source scripts/setup_garfield.sh            # set up
#   bash   scripts/setup_garfield.sh --build    # build the pin into ~/garfield
#   bash   scripts/setup_garfield.sh --check    # report what would be used
#
# ─────────────────────────────────────────────────────────────────────────────
# WHICH GARFIELD DO YOU ACTUALLY NEED? Read this before building anything.
#
# There are two answers and the difference is ~15 minutes of your time.
#
# **For Magboltz gas tables (P0.10) the CVMFS Garfield is FINE.** Magboltz is
# vendored inside Garfield at version 11.19 (January 2024) and has not moved:
# between the LCG_108 Garfield and the current pin, `Magboltz/magboltz.f`
# differs only by the fixed-form continuation marker in column 6 (354 lines of
# `/` -> `&`) plus one missing comma in a FORMAT *print* statement. No cross
# section, no transport change. Garfield's built-in Penning table is likewise
# unchanged — MX17 re-probed it and every Penning rate reproduced exactly.
# So if what you need is v_d, D_L, D_T and Townsend/attachment coefficients to
# replace Stage B's placeholders, stop here and use the LCG view.
#
# **For field solves and induction you need a recent build.** LCG_108 ships
# Garfield 664 commits behind and LCG_109 281 behind, and the APIs all still
# exist — which is exactly why this bites silently rather than failing. What
# landed between March and August 2026:
#
#   * the neBEM OpenMP race in the SVD inversion — wrong field solves on a
#     multi-core box, no error;
#   * interface-crossing checks, i.e. electrons no longer tunnel through mesh
#     wires. That *is* the mesh-transparency observable (P0.13); on an old
#     build the number is wrong and looks fine;
#   * `AvalancheMicroscopic::GetIons()`, needed for the ion tail in P0.17;
#   * the `Examples/ResistiveMicromegas` example;
#   * the FFT-convolution fix and arbitrary-PSD noise generators;
#   * the regression test suite itself.
#
# Provenance: MX17_Geant/design/RESPONSE_SIM_PLAN.md §5a, which is where the
# pin and the Magboltz diff were established.
# ─────────────────────────────────────────────────────────────────────────────

P2_GARFIELD_PIN="${P2_GARFIELD_PIN:-927e5c21}"
P2_LCG_VIEW="${P2_LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc14-opt}"

# An install someone else already made. Optional, and deliberately NOT the
# default: it is another user's AFS area, so it depends on their ACLs, their
# account continuing to exist, and AFS continuing to exist. Set it if you have
# been given read access; otherwise --build gives you your own in ~15 min.
#
#   Known shared install (Dylan, lxplus):
#     /afs/cern.ch/user/d/dneff/work/garfield_install/lcg109-927e5c21
#   To grant access to it, the OWNER runs (recursively, per directory):
#     fs setacl -dir <path> -acl <cern-username> read
P2_SHARED_INSTALL="${P2_SHARED_INSTALL:-}"

_mode="setup"
case "$1" in
    --build) _mode="build" ;;
    --check) _mode="check" ;;
esac

# ── Resolve ──────────────────────────────────────────────────────────────────
_gf=""
_how=""
if   [ -n "$P2_GARFIELD_INSTALL" ] && [ -d "$P2_GARFIELD_INSTALL/include/Garfield" ]; then
    _gf="$P2_GARFIELD_INSTALL"; _how="P2_GARFIELD_INSTALL"
elif [ -d "$PWD/garfield/include/Garfield" ]; then
    _gf="$PWD/garfield"; _how="condor scratch tarball"
elif [ -n "$P2_SHARED_INSTALL" ] && [ -r "$P2_SHARED_INSTALL/include/Garfield" ]; then
    _gf="$P2_SHARED_INSTALL"; _how="shared install"
elif [ -d "$HOME/garfield/install/include/Garfield" ]; then
    _gf="$HOME/garfield/install"; _how="personal build"
fi

_have_lcg=0
[ -f "$P2_LCG_VIEW/setup.sh" ] && _have_lcg=1

if [ "$_mode" = "check" ]; then
    echo "P2 Garfield environment check"
    echo "  host            : $(hostname -s)"
    echo "  LCG view        : $([ $_have_lcg = 1 ] && echo "$P2_LCG_VIEW" || echo 'NOT FOUND (not on CVMFS?)')"
    echo "  pin wanted      : $P2_GARFIELD_PIN"
    if [ -n "$_gf" ]; then
        echo "  Garfield        : $_gf   [$_how]"
    else
        echo "  Garfield        : NONE FOUND"
        echo "                    For Magboltz gas tables (P0.10) you do not need one —"
        echo "                    source the LCG view and use its Garfield."
        echo "                    For field solves / induction: bash scripts/setup_garfield.sh --build"
    fi
    if [ -n "$P2_SHARED_INSTALL" ] && [ ! -r "$P2_SHARED_INSTALL/include/Garfield" ]; then
        echo "  NOTE: P2_SHARED_INSTALL is set but not readable — likely an AFS"
        echo "        permission issue. Ask the owner for 'fs setacl ... read',"
        echo "        or just build your own (--build)."
    fi
    return 0 2>/dev/null || exit 0
fi

# ── Build ────────────────────────────────────────────────────────────────────
if [ "$_mode" = "build" ]; then
    set -e
    _src="${P2_GARFIELD_SRC:-$HOME/garfield/src}"
    _inst="${P2_GARFIELD_INSTALL:-$HOME/garfield/install}"
    echo "==> Building Garfield++ @ $P2_GARFIELD_PIN"
    echo "    source : $_src"
    echo "    install: $_inst"
    if [ $_have_lcg = 1 ]; then
        echo "==> Sourcing $P2_LCG_VIEW for ROOT/compiler"
        source "$P2_LCG_VIEW/setup.sh"
    else
        echo "!!  No LCG view; relying on whatever ROOT/cmake are already set up."
    fi
    mkdir -p "$(dirname "$_src")"
    if [ ! -d "$_src/.git" ]; then
        git clone https://gitlab.cern.ch/garfield/garfieldpp.git "$_src"
    fi
    git -C "$_src" fetch --all --tags
    git -C "$_src" checkout "$P2_GARFIELD_PIN"
    cmake -S "$_src" -B "$_src/../build" -DCMAKE_INSTALL_PREFIX="$_inst" \
          -DCMAKE_BUILD_TYPE=Release
    cmake --build "$_src/../build" -j"$(nproc)"
    cmake --install "$_src/../build"
    echo ""
    echo "==> Done. Now:  source scripts/setup_garfield.sh"
    echo "    Upstream examples (incl. ResistiveMicromegas) live in $_src/Examples"
    exit 0
fi

# ── Setup ────────────────────────────────────────────────────────────────────
if [ $_have_lcg = 1 ]; then
    source "$P2_LCG_VIEW/setup.sh"
fi

if [ -z "$_gf" ]; then
    echo "[setup_garfield] No pinned Garfield++ install found." >&2
    if [ $_have_lcg = 1 ]; then
        echo "[setup_garfield] The LCG view IS set up, and its Garfield is" >&2
        echo "[setup_garfield] adequate for MAGBOLTZ GAS TABLES (P0.10) — see the" >&2
        echo "[setup_garfield] header of this script for why, and for what it is" >&2
        echo "[setup_garfield] NOT adequate for (field solves, induction)." >&2
        return 0 2>/dev/null || exit 0
    fi
    echo "[setup_garfield] Build one:  bash scripts/setup_garfield.sh --build" >&2
    return 1 2>/dev/null || exit 1
fi

export GARFIELD_INSTALL="$_gf"
export LD_LIBRARY_PATH="$_gf/lib64:$_gf/lib:$LD_LIBRARY_PATH"
export HEED_DATABASE="${HEED_DATABASE:-$_gf/share/Heed/database}"
for _pd in "$_gf"/lib64/python*/site-packages "$_gf"/lib/python*/site-packages; do
    [ -d "$_pd" ] && export PYTHONPATH="$_pd:$PYTHONPATH"
done

echo "[setup_garfield] host=$(hostname -s)  Garfield=$_gf  [$_how]"
[ -n "$(command -v root-config)" ] && \
    echo "[setup_garfield] ROOT $(root-config --version)"
