#!/usr/bin/env bash
# HTCondor wrapper for SolveP2MeshCell (COMSOL 6.4 from CVMFS, CERN licence).
# Temporary and recovery files go to the job scratch: the AFS home is too
# small for COMSOL's model copies. Outputs are written to EOS by the class.
set -e
C=/cvmfs/projects.cern.ch/engtools/comsol/comsol64/multiphysics/bin/comsol
mkdir -p /eos/user/a/akallits/p2_sim/comsol/p2_sps/h4mm /eos/user/a/akallits/p2_sim/comsol/p2_sps/h2mm
cd "$_CONDOR_SCRATCH_DIR"   # HTCondor has already copied the class here
$C batch -np 8 -nosave -tmpdir "$_CONDOR_SCRATCH_DIR" \
    -recoverydir "$_CONDOR_SCRATCH_DIR" \
    -inputfile "$(basename "$1")" -batchlog solve.log
cat solve.log
