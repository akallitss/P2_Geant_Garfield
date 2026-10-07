#!/usr/bin/env bash
# submit_gas_tables.sh -- one HTCondor job per gas: Magboltz table (P0.10).
#
#   bash scripts/garfield/submit_gas_tables.sh [outdir] [ncoll]
#
# Default outdir /eos/user/a/akallits/p2_sim/gas_tables, ncoll 10.
# Garfield comes from the LCG view (enough for Magboltz, see
# scripts/setup_garfield.sh). ~1-2 h per gas on one core.
set -e
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${1:-/eos/user/a/akallits/p2_sim/gas_tables}"
NCOLL="${2:-10}"
GASES="ArIso ArIso9010 NeIso NeIso8515 ArCO2Iso9352 ArCF4Iso"
JOBDIR="$REPO/condor_gas"
mkdir -p "$OUT" "$JOBDIR/logs"

cat > "$JOBDIR/run_gas.sh" <<EOF
#!/usr/bin/env bash
set -e
source "$REPO/scripts/setup_garfield.sh"
cd "\$_CONDOR_SCRATCH_DIR"
root -l -b -q "$REPO/scripts/garfield/make_gas_table.C(\"\$1\", \"\$1\", $NCOLL)"
cp "\$1.gas" "$OUT/"
EOF
chmod +x "$JOBDIR/run_gas.sh"

cat > "$JOBDIR/gas.sub" <<EOF
executable   = $JOBDIR/run_gas.sh
arguments    = \$(gas)
output       = $JOBDIR/logs/\$(gas).out
error        = $JOBDIR/logs/\$(gas).err
log          = $JOBDIR/logs/condor.log
+JobFlavour  = "tomorrow"
request_cpus = 1
requirements = (OpSysAndVer =?= "AlmaLinux9")
queue gas in ($GASES)
EOF

echo "gases : $GASES"
echo "out   : $OUT"
condor_submit "$JOBDIR/gas.sub"
