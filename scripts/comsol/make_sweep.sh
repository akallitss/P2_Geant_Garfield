#!/usr/bin/env bash
# make_sweep.sh E1 E2 ... -- one COMSOL class per drift field [V/cm] for the
# h2mm cell (VCathode = -450 V - E * 0.2 cm) and one HTCondor job each.
# MH = wire centre plane [um] (default 150, the colleague's value: the mesh
# bottom is then at ~135 um); MH=164.85 puts the bottom of the calendered
# mesh (flattened at +-1.5 rCut = +-14.85 um) on the 150 um pillars.
# NOSUBMIT=1 only compiles (submission by a DAG).
# TAG is appended to the output folder (h2mm_sweep<TAG>/E<E>).
# PITCH = wire spacing [um] (default 68, measured P2), VM = |mesh voltage| [V]
# (default 450), RW = wire radius [um] (default 9; the calendering cut and the
# thickness scale with it). Nikolopoulos et al. 2011 geometry for the validation:
#   PITCH=63.5 MH=142.85 VM=563 TAG=_nik bash make_sweep.sh ...
# The meshes are re-generated per solve and differ slightly (4.00 M nodes,
# not node-for-node identical), so fields are solved, not superposed.
#   cd scripts/comsol && bash make_sweep.sh 150 300 450 625 900 1250 1800 2500
set -e
MH="${MH:-150}"
PITCH="${PITCH:-68}"
VM="${VM:-450}"
RW="${RW:-9}"
TAG="${TAG:-}"
C=/cvmfs/projects.cern.ch/engtools/comsol/comsol64/multiphysics/bin/comsol
mkdir -p sweep logs
for E in "$@"; do
    VC=$(echo "-$VM - $E * 0.2" | bc -l | sed 's/0*$//; s/\.$//')
    cls="SolveSweep${TAG}_E$E"
    sed "s/@E@/$E/g; s/@VC@/$VC/g; s/@MH@/$MH/g; s/@TAG@/$TAG/g; s/@PITCH@/$PITCH/g; s/@VM@/$VM/g; s/@RW@/$RW/g; s#@OUT@#h2mm_sweep$TAG/E$E#g; s/@H@/2/g; s/SolveSweep_E$E/$cls/g" \
        SolveSweep.java.in > "sweep/$cls.java"
    (cd sweep && $C compile "$cls.java")
    mkdir -p "/eos/user/a/akallits/p2_sim/comsol/p2_sps/h2mm_sweep$TAG/E$E"
    [ -n "$NOSUBMIT" ] || condor_submit -a "cls=sweep/$cls" mesh_cell.sub
done
