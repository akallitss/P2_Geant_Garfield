# How-to: recipes for working on your own

Each recipe says when you need it, the commands, and how to check the result. Install is in
`SETUP.md`; the physics state is in `HANDOVER.md` and `OPEN_PROBLEMS.md`. Everything here
assumes you have run `source scripts/setup_feynman.sh` (Feynman) or `conda activate p2sim &&
source scripts/setup_local.sh` (laptop) in the current shell.

**Rule of thumb:** anything heavier than one Geant4 or Garfield process goes to SLURM on
Feynman, not to a laptop. Each Garfield transparency/gain process needs ~3 GB of memory (measured peak, 2026-11-01): five at once crashed a 14 GB laptop.

---

## 0. Your own environment, from nothing

Do these once, in this order. Nobody else's account or files are needed afterwards.

1. **Accounts:**
   - A Feynman account (IRFU computing). Ask the admins for a `/feynman/work` or
     `/feynman/scratch` directory too: your home quota (6 GB) is too small.
   - A GitHub account. Ask Alexandra to add you as a collaborator on
     `akallitss/P2_Geant_Garfield`, or to transfer the repo to a group organisation.
2. **SSH:** a key on Feynman and on GitHub (`ssh-keygen`, then add `~/.ssh/id_ed25519.pub`
   to GitHub → Settings → SSH keys). Clone with the SSH URL so you can push:
   `git clone git@github.com:akallitss/P2_Geant_Garfield.git`.
3. **Install:** `SETUP.md` §1 (Feynman) or §2 (laptop).
4. **Data:** download the bundle from the CERNBox link into `$P2_DATA` (`SETUP.md` §3),
   including `sps_products.tar`.
5. **Check:** run every check in `SETUP.md` §4. Do this again whenever you rebuild or move
   to a new machine.
6. **Optional:** with Claude Code, start it in the repo. It reads `CLAUDE.md` automatically,
   so the working rules apply to it too.

## 1. Run Geant4

**One point, interactively** (seconds to minutes):
```bash
./build/mm_sim -m p2 -g ArIso -p electron -e 100 -n 1000 -t 1 \
    --heed --no-pixe --pai config --beam-spread 11.43 -o $P2_DATA/test/e100
./build/mm_sim --list-gases     # available mixtures
./build/mm_sim -h               # all options
```

**A scan** (the standard way, which writes a manifest of what you asked for):
```bash
python3 scripts/submit_condor_p2.py --scan rnd2026 --slurm --dry-run      # list the jobs
python3 scripts/submit_condor_p2.py --scan rnd2026 --slurm \
    --outdir $P2_DATA/rnd2026_v2 --jobdir $P2_WORK/jobs/rnd2026_v2
squeue -u $USER                                                            # watch
```
`--only <part of a tag>` selects jobs and `--scale 0.1` shrinks them, which is good for tests.
To add points or a new scan, edit the lists at the top of `scripts/submit_condor_p2.py`
(`RND_GASES`, `RND_DRIFT_GAPS`, `RND_PHOTON_ENERGIES`, …) or add a block in `build_jobs()`.

**Check:** every job log `logs/<tag>.out` ends with `Job done`. Then:
```bash
python3 scripts/collect_results.py --indir $P2_DATA/rnd2026_v2    # merges, summary.csv, warnings
```

## 2. Read a Geant4 output file

```python
import uproot, numpy as np
f = uproot.open("p2_ArIso_3p964mm_electron100MeV_th0_t0.root")
meta = f["RunMeta"].arrays(library="np")
assert meta["heed"][0] == 1 and meta["pixe"][0] == 0         # production physics
n = f["EventTree"]["nPrimDrift"].array(library="np")
thrown = meta["thrown"].sum()                                 # NOT len(n): --skip-empty
print("P(n<10) =", np.sum(n < 10) / thrown)
```
Tree contents: `OUTPUT_FORMAT.md`. Pitfalls: `SETUP.md` §5.

To compare two files of the same point (e.g. after a code change):
`python3 scripts/check_local_vs_reference.py new.root reference.root`.

## 3. Mesh transparency (Garfield++)

**One point:**
```bash
bash scripts/garfield/run_transparency.sh transp \
    $P2_DATA/comsol/p2_sps/h2mm_sweep_gap150/E625 $P2_DATA/gas_tables/ArIso.gas \
    $P2_DATA/transp/E625_ArIso 5000 1 null
```
Arguments are `<mode> <cell dir> <gas table> <output prefix> <n electrons> <seed> <flags>`.
**Always pass `null`.** It prints one `TRANSPARENCY gas=… micro=T +- err … wire=… attached=…` line.
At 5000 electrons it takes ~5 min.

**A sweep on SLURM:** put one command per line in a file, then submit it as an array.
```bash
J=$P2_WORK/jobs/sweep_ArIso; mkdir -p $J/logs; : > $J/cmds.txt
for E in 150 300 450 625 900 1250 1800 2500; do
  echo "source $P2_REPO/scripts/setup_feynman.sh && bash $P2_REPO/scripts/garfield/run_transparency.sh transp \
$P2_DATA/comsol/p2_sps/h2mm_sweep_gap150/E$E $P2_DATA/gas_tables/ArIso.gas $P2_DATA/transp/E${E}_ArIso 5000 1 null" >> $J/cmds.txt
done
sbatch --array=1-$(wc -l < $J/cmds.txt) -p htc --mem-per-cpu=4000 -o "$J/logs/%A_%a.out" scripts/slurm_array.sh $J/cmds.txt
```
Garfield jobs need `--mem-per-cpu=4000` (they peak at ~3 GB). Use that flag, not `--mem`, which
clashes with the script's own setting. One 1000-electron run has ±0.008 statistical error:
use ≥ 5000 electrons, or several seeds, per point.

**Collect the results into a table:**
```bash
python3 scripts/garfield/collect_lines.py TRANSPARENCY $J/logs/*.out -o sweep_ArIso.csv
```
This gives one row per line, with `gas, E_drift, cell_tag, micro, micro_err, wire, attached, …`.
The same works for the finished geometry runs: `$P2_DATA/condor_transp/geom_*.out`.

**Seeds:** for more statistics, run several jobs with different seeds (the 6th argument) and
combine them. The same seed on two different builds does *not* give the same electrons.

## 4. Gain with Penning transfer

```bash
bash scripts/garfield/run_gain.sh $P2_DATA/comsol/p2_sps/h2mm_sweep_gap150/E625 \
    $P2_DATA/gas_tables/ArCO2Iso9352.gas $P2_DATA/gain/E625 0,0.2,0.36,0.5 20
python3 scripts/garfield/collect_lines.py GAIN <logs> -o gain.csv
```
Arguments: `<cell> <gas table> <out> <Penning r values> <n×n start grid>`. The fitted r is ≈ 0.36
(`results/gain/`). Gain changes by e per ~30 V of mesh voltage in the simulation, but per 22–28 V
in data. That mismatch is still open.

## 5. Make a Magboltz gas table

Needed for a new mixture, or a new contamination level. It takes 1–2 h per gas on one core, so
run it as a SLURM job:
```bash
echo "source $P2_REPO/scripts/setup_feynman.sh && cd $P2_REPO && root -l -b -q \
'scripts/garfield/make_gas_table.C(\"ArIso9010\", \"$P2_DATA/gas_tables_new/ArIso9010\", 10)'" > cmds.txt
sbatch --array=1-1 -p htc --mem-per-cpu=4000 -o "gas_%A_%a.out" scripts/slurm_array.sh cmds.txt
```
Then export all tables to the JSON Stage B reads:
```bash
root -l -b -q 'scripts/garfield/export_gas_json.C("'$P2_DATA'/gas_tables_new", "stage_b/gas_tables/magboltz_new.json", "ArIso9010")'
```
**Check:** `python3 scripts/analysis/gas_vd_dt_map.py` plots v_d and D_T for every gas. Compare
the new gas with its neighbours and with the literature before using it.

## 6. Add a new gas mixture

The mixture list is defined in several places, and they must agree (names and fractions):

| file | what |
|---|---|
| `src/GasMixtures.cc` | Geant4 material and W value. The authority: `mm_sim --list-gases` |
| `scripts/garfield/make_gas_table.C` | Magboltz composition |
| `scripts/garfield/mesh_transparency.cc`, `heed_mip.cc` | Garfield composition (`kMixtures`) |
| `stage_b/gas.py` | Stage B gas list |
| `scripts/submit_condor_p2.py` | scan lists (`CAMPAIGN_GASES`, `RND_GASES`) |

Then rebuild (`bash scripts/build.sh`), make its gas table (§5), and run a transparency point.

## 7. Change a geometry number

1. Find the source (drawing, gerber, or a written answer from the detector team) and record it in
   `P2_GEOMETRY.md`. **No source, no change.** Ask instead.
2. Change the default in `include/SimConfig.hh` (or `include/P2Wedge.hh`) and keep
   `scripts/model/p2_model.py` and `P2_MODEL.md` in sync.
3. Rebuild, then run one point. `RunMeta.geometryHash` changes. That is expected, and it means
   older files can't be mixed with new ones.
4. Regenerate the figures if they change: `python3 scripts/model/plot_p2_model.py`.

A **new mesh or gap geometry for Garfield** needs a new COMSOL field map (`SETUP.md` §6).

## 8. Stage B and the VMM on Geant4 output

```bash
python3 -m stage_b.run $P2_DATA/rnd2026_heed/p2_ArIso_3p964mm_electron100MeV_th0_t0.root \
    -o $P2_DATA/stageb/e100_ArIso.root --gain 1e4 --mesh-transparency 0.943 --drift-field 625
python3 -m stage_b.selftest
```
Options: `python3 -m stage_b.run -h`. Every run writes a `StageBMeta` tree with all its
assumptions. Joining Stage B to the VMM emulator (`vmm/vmm_emulator.py`) is the open
deliverable (`OPEN_PROBLEMS.md` §2, item 5). Read `stage_b/README.md` and `vmm/README.md` first.

## 9. Compare with SPS or cosmic data

- **Data products** come from `sps_products.tar` in `$P2_DATA/sps_products/` (Landau MPV vs HV per
  scan and detector, and the corrected `run_config.json` of each run). The analysis scripts find
  them through `scripts/analysis/p2paths.py`.
- **Gas and HV of a run:** read them from the run config, never from memory or the logbook:
  ```bash
  python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(d['gas'], '|', d.get('gas_note',''))" \
      $P2_DATA/sps_products/runs/drift_mesh_scan_1/run_config.json
  ```
- Existing comparisons to start from:
  - `scripts/analysis/sps_drift_scan_vs_sim.py`: transparency vs the drift scan
  - `sps_contamination_fit.py stage_b/gas_tables/magboltz_contam_20261014.json $P2_DATA/transparency/contam`
  - `gain_vs_data.py`
  - `sps_drift_velocity_check.py`
- Drift field = (V_drift − V_mesh) / 0.4 cm for the 4 mm gap.
- Compare shapes and tails, not just means (`../CLAUDE.md`).

## 10. Record a result

1. Put the output table and figure in `docs/results/<topic>/`, with a short `README.md`: what it
   is, the exact command, the git commit, the SLURM job id, and the input data.
2. Add a dated slide or entry to the logbook (`docs/logbook/project/`).
3. Commit on a branch and push it:
   `git switch -c <topic>`, `git add …`, `git commit`, `git push -u origin <topic>`.
   Then open a pull request so someone can see it.
4. Production runs must come from a committed, clean checkout (`RunMeta.gitDirty = 0`).

## 11. Troubleshooting

| symptom | cause and fix |
|---|---|
| `build.sh` prints `--heed disabled` | Garfield not found. Run `source scripts/setup_feynman.sh` (or `setup_local.sh` after building Garfield), then `build.sh` again |
| link error `undefined reference to Garfield::Random::draw` | an old checkout; `CMakeLists.txt` now links `libGarfieldRandom`. Update the repo |
| CMake: `Could NOT find EXPAT … 2.6.1` (laptop) | the system expat was found first. Run `source scripts/setup_local.sh` (it puts conda first) and make sure `expat` is installed in the env |
| `source scripts/setup_feynman.sh` takes minutes | CVMFS is caching on first use. Later shells are faster |
| `Pinned Garfield++ not found` from `run_*.sh` | build it (`bash scripts/setup_garfield.sh --build`), or point `P2_GARFIELD_INSTALL` at it. On Feynman, set `P2_WORK` **before** `source scripts/setup_feynman.sh`, because the install path is derived from it. The scripts refuse to fall back to another Garfield on purpose |
| `FileNotFoundError … P2_SPS_PRODUCTS` | extract `sps_products.tar` into `$P2_DATA`, or set the variable |
| `Disk quota exceeded` on Feynman | you are writing to home. Set `P2_WORK`/`P2_DATA` to your work area |
| `sbatch: --mem, --mem-per-cpu … mutually exclusive` | use `--mem-per-cpu=…` |
| SLURM job failed | `logs/slurm_<id>_<task>.out` and `logs/<tag>.out`. The command line is printed at the top |
| laptop freezes or reboots | too many Garfield processes at once (each needs ~3 GB). Use SLURM |
| `check_local_vs_reference.py` says `DIFFERENT` | read which RunMeta field differs. A different `geometryHash` means the geometry code changed since the reference, so the comparison is not like-for-like |
