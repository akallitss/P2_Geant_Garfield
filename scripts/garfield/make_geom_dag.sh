#!/usr/bin/env bash
# make_geom_dag.sh -- geometry study (2026-10-14): is the SPS / cosmic drift-scan
# discrepancy a mesh / amplification-gap effect? Compiles one COMSOL class per
# cell (scripts/comsol/make_sweep.sh, NOSUBMIT) and writes geom.dag: each
# COMSOL solve, then the Garfield jobs on that cell.
#
# Variations of the corrected cell (mesh bottom on the pillars, mesh -450 V,
# 2 mm cell), each over the drift fields $FIELDS, transparency for
# Ar/CO2/iso 93/5/2 (SPS) and Ar/iso 95/5 (cosmic):
#   _g120 _g135 _g165   amplification gap 120 / 135 / 165 um (base 150)
#   _rw8 _rw10          wire diameter 16 / 20 um (base 18), spacing 68 um
#   _p66 _p73           wire spacing 66 / 73 um = opening 48 / 55 um (base 50)
#   _vm415              base geometry at mesh -415 V (cosmic det1), Ar/iso only
# Gain: at 625 V/cm, gaps 120 / 135 / 150 / 165 um x mesh 430 / 450 / 470 V,
# both gases (mesh_gain.cc), plus every variation's 625 V/cm cell.
# Base cells at 450 V already exist (h2mm_sweep_gap150).
# The bottom of the calendered mesh sits on the pillars: MH = gap + 1.65 rWire.
#
#   bash scripts/garfield/make_geom_dag.sh      (on lxplus, from the repo root)
#   condor_submit_dag -outfile_dir condor_transp scripts/garfield/geom.dag
set -e
cd "$(dirname "$0")/../.."
FIELDS="${FIELDS:-150 300 450 625 900 1250 1800}"
CS=/eos/user/a/akallits/p2_sim/comsol/p2_sps
TR=/eos/user/a/akallits/p2_sim/transparency/geom
GN=/eos/user/a/akallits/p2_sim/gain/geom
mkdir -p condor_transp "$TR" "$GN"
dag=scripts/garfield/geom.dag
echo "# geometry study DAG, written by make_geom_dag.sh $(date -I)" > $dag
seed=100

mh() { echo "$1 + 1.65 * $2" | bc -l | sed 's/0*$//; s/\.$//'; }

# cell TAG MH PITCH VM RW "fields" "gases" gain(0/1)
cell() {
    local tag=$1 m=$2 p=$3 vm=$4 rw=$5 fields=$6 gases=$7 gain=$8
    (cd scripts/comsol && NOSUBMIT=1 MH=$m PITCH=$p VM=$vm RW=$rw TAG=$tag bash make_sweep.sh $fields >/dev/null)
    for E in $fields; do
        local c="C${tag}_E$E"
        echo "JOB $c mesh_cell.sub DIR scripts/comsol" >> $dag
        echo "VARS $c cls=\"sweep/SolveSweep${tag}_E$E\"" >> $dag
        for g in $gases; do
            local n="T${tag}_E${E}_$g"
            seed=$((seed + 1))
            echo "JOB $n scripts/garfield/transparency_geom.sub" >> $dag
            echo "VARS $n cell=\"$CS/h2mm_sweep$tag/E$E\" gas=\"$g\" out=\"$TR/${tag#_}_E${E}_$g\" name=\"${tag#_}_E${E}_$g\" seed=\"$seed\"" >> $dag
            echo "PARENT $c CHILD $n" >> $dag
            if [ "$gain" = 1 ] && [ "$E" = 625 ]; then
                n="N${tag}_E${E}_$g"
                echo "JOB $n scripts/garfield/gain_geom.sub" >> $dag
                echo "VARS $n cell=\"$CS/h2mm_sweep$tag/E$E\" gas=\"$g\" out=\"$GN/${tag#_}_$g\" name=\"${tag#_}_$g\"" >> $dag
                echo "PARENT $c CHILD $n" >> $dag
            fi
        done
    done
}

GASES="ArCO2Iso9352 ArIso"
for g in 120 135 165; do cell _g$g "$(mh $g 9)" 68 450 9 "$FIELDS" "$GASES" 1; done
cell _rw8  "$(mh 150 8)"  68 450 8  "$FIELDS" "$GASES" 1
cell _rw10 "$(mh 150 10)" 68 450 10 "$FIELDS" "$GASES" 1
cell _p66  "$(mh 150 9)"  66 450 9  "$FIELDS" "$GASES" 1
cell _p73  "$(mh 150 9)"  73 450 9  "$FIELDS" "$GASES" 1
cell _vm415 "$(mh 150 9)" 68 415 9  "$FIELDS" "ArIso" 1
# gain slope: mesh 430 / 470 V at 625 V/cm for every gap (transparency too)
for vm in 430 470; do
    for g in 120 135 150 165; do cell _g${g}_vm$vm "$(mh $g 9)" 68 $vm 9 625 "$GASES" 1; done
done
# base cell at 450 V (already solved): gain only
for g in $GASES; do
    echo "JOB N_g150_E625_$g scripts/garfield/gain_geom.sub" >> $dag
    echo "VARS N_g150_E625_$g cell=\"$CS/h2mm_sweep_gap150/E625\" gas=\"$g\" out=\"$GN/g150_$g\" name=\"g150_$g\"" >> $dag
done
echo "wrote $dag: $(grep -c '^JOB C' $dag) COMSOL, $(grep -c '^JOB T' $dag) transparency, $(grep -c '^JOB N' $dag) gain jobs"
