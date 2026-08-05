#!/usr/bin/env bash
# benchmark_lxplus.sh
# Measures per-event CPU cost and output size of mm_sim, and fits the
# fixed startup cost. Results and their interpretation: docs/BENCHMARKS.md.
#
# Usage (from the repo root, on lxplus, after building):
#   bash scripts/benchmark_lxplus.sh            # full suite (~15 min)
#   bash scripts/benchmark_lxplus.sh quick      # just the cost-model fit
#
# IMPORTANT: lxplus login nodes are shared and typically sit at load 9-15,
# so wall-clock timings there are meaningless. Everything below reports
# **user CPU time**, which is only mildly affected by contention. Do not run
# two copies of this at once (it invalidates both).

set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${REPO}/build/mm_sim"
WORK="${WORK:-/tmp/$(whoami)_mmbench}"

if [[ ! -x "$BIN" ]]; then
    echo "ERROR: $BIN not found. Run: source scripts/setup_lxplus.sh && bash scripts/build.sh"
    exit 1
fi

rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
echo "# host=$(hostname) nproc=$(nproc) load=$(uptime | sed 's/.*load average: //')"

# --- cost-model scan: fixed startup + marginal per-event cost -------------
echo
echo "## cost model (single thread; fit CPU = init + marginal*N)"
echo "particle|E_MeV|N|cpu_user_s"
fit_run() {
    local P=$1 E=$2 N=$3; shift 3
    local R
    R=$( { /usr/bin/time -f "%U" nice -n 5 "$BIN" -p "$P" -e "$E" -n "$N" \
           -t 1 -o fit "$@" >/dev/null; } 2>&1 | tail -1 )
    echo "$P|$E|$N|$R"
}
for n in 200 1000 5000 20000; do fit_run electron 119 $n;             done
for n in 5000 50000 200000;   do fit_run gamma 0.060 $n --skip-empty; done

[[ "${1:-}" == "quick" ]] && exit 0

# --- throughput / output size across the campaign's run points ------------
echo
echo "## run points (single thread unless noted)"
echo "label|particle|E_MeV|N|threads|cpu_user_s|cpu_sys_s|wall_s|bytes|written"
run() {
    local L=$1 P=$2 E=$3 N=$4 T=$5; shift 5
    local R SZ WR
    R=$( { /usr/bin/time -f "%U %S %e" nice -n 5 "$BIN" -p "$P" -e "$E" \
           -n "$N" -t "$T" -o "b_$L" "$@" > "log_$L.txt"; } 2>&1 | tail -1 )
    SZ=$(du -cb b_${L}_t*.root 2>/dev/null | tail -1 | cut -f1)
    WR=$(grep -o "written=[0-9]*" "log_$L.txt" | cut -d= -f2 | paste -sd+ | bc 2>/dev/null)
    [[ -z "$WR" ]] && WR=$N
    echo "$L|$P|$E|$N|$T|$(echo "$R" | tr ' ' '|')|$SZ|$WR"
}
run g8       gamma    0.008  100000 1 --skip-empty
run g60      gamma    0.060  100000 1 --skip-empty
run g100     gamma    0.100  100000 1 --skip-empty
run g150     gamma    0.150  100000 1 --skip-empty
run g60_ne   gamma    0.060  100000 1 --skip-empty -g NeIso
run g60_nosk gamma    0.060  100000 1
run e30      electron 30     2000   1
run e119     electron 119    2000   1
run e155     electron 155    2000   1
run e119_ne  electron 119    2000   1 -g NeIso
run mu200    muon     200000 2000   1
run g60_t8   gamma    0.060  400000 8 --skip-empty
run e119_t8  electron 119    8000   8

echo "# load at end: $(uptime | sed 's/.*load average: //')"
echo "# outputs left in $WORK"
