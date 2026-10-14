#!/usr/bin/env bash
# make_gaps.sh G1 G2 ... -- one COMSOL class per drift gap [mm], corrected
# geometry (mesh bottom on the 150 um pillars, meshHeight 164.85 um), mesh
# -450 V and the drift electrode a fixed DV above it (default 250 V, as at
# SPS), so the drift field is DV / gap. Outputs:
#   /eos/user/a/akallits/p2_sim/comsol/p2_sps/gaps_gap150/h<G>mm/
# Writes sweep/SolveGap_h<G>mm.class; submission is done by the DAG
# (scripts/garfield/gaps.dag), so the Garfield jobs follow each solve.
#   cd scripts/comsol && bash make_gaps.sh 4 5 6 7 8
set -e
C=/cvmfs/projects.cern.ch/engtools/comsol/comsol64/multiphysics/bin/comsol
DV="${DV:-250}"
mkdir -p sweep logs
for G in "$@"; do
    VC=$(echo "-450 - $DV" | bc -l)
    cls="SolveGap_h${G}mm"
    sed "s/@E@/$G/g; s/@VC@/$VC/g; s/@MH@/164.85/g; s/@PITCH@/68/g; s/@VM@/450/g; \
         s#@OUT@#gaps_gap150/h${G}mm#g; s/@H@/$G/g; s/SolveSweep_E$G/$cls/g" \
        SolveSweep.java.in > "sweep/$cls.java"
    grep -q "class $cls" "sweep/$cls.java" || { echo "class name not set"; exit 1; }
    (cd sweep && $C compile "$cls.java")
    mkdir -p "/eos/user/a/akallits/p2_sim/comsol/p2_sps/gaps_gap150/h${G}mm"
done
