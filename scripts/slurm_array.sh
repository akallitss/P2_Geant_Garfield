#!/usr/bin/env bash
# slurm_array.sh — run line $SLURM_ARRAY_TASK_ID of a command file.
# The SLURM counterpart of the HTCondor .sub/.dag files: one line = one job.
#
#   N=$(wc -l < cmds.txt)
#   mkdir -p logs
#   sbatch --array=1-$N%200 -p htc -o 'logs/%A_%a.out' scripts/slurm_array.sh cmds.txt
#
# submit_condor_p2.py --slurm writes the command file and calls sbatch for
# you. For Garfield sweeps, write one run_transparency.sh / run_gain.sh call
# per line (arguments as in the scripts/garfield/*.sub files).
#SBATCH --job-name=p2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
# 2 GB is enough for Geant4. Garfield transparency/gain jobs peak at ~3 GB
# (the COMSOL field map): submit those with --mem-per-cpu=4000.
set -eo pipefail
cmdfile="$1"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$cmdfile")
if [ -z "$line" ]; then
    echo "no line ${SLURM_ARRAY_TASK_ID} in $cmdfile" >&2
    exit 2
fi
echo "[slurm_array] $(hostname) task ${SLURM_ARRAY_TASK_ID}: $line"
eval "$line"
