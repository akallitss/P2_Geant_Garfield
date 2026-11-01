# Setup: running the P2 chain

No CERN account is needed. Use **Feynman** (CEA/IRFU cluster) for production runs, and your
laptop if you like for development and plots. The only thing neither can do is a new COMSOL
field solve (§6).

## 1. Feynman (recommended)

Feynman mounts `/cvmfs/sft.cern.ch`. Its LCG_106 release has an el8 build with **exactly the
production stack**: Geant4 11.2.1, ROOT 6.32.02 and Garfield++ d3bbb29, the same as `mm_sim` used on
lxplus for `rnd2026_heed`. So you build the same program, not an approximation of it.

**Storage first.** The home quota is 6 GB, which is too small. Ask the Feynman admins for a
`/feynman/work` (or `/feynman/scratch`) directory, or use the group project area, and point
`P2_WORK` at it. Everything below goes there.

```bash
ssh feynman
export P2_WORK=/feynman/work/<...>/p2            # put this line in your ~/.bashrc
mkdir -p $P2_WORK && cd $P2_WORK
git clone https://github.com/akallitss/P2_Geant_Garfield.git
cd P2_Geant_Garfield
source scripts/setup_feynman.sh                  # every new shell; the first time can take ~3 min (CVMFS)
bash scripts/build.sh                            # -> build/mm_sim, "--heed enabled"
bash scripts/setup_garfield.sh --build           # pinned Garfield for transparency/gain, ~20 min, once
```

**Running jobs (SLURM).** `submit_condor_p2.py --slurm` submits a scan as one SLURM array on
the `htc` partition (28 nodes × 128 cores):

```bash
python3 scripts/submit_condor_p2.py --scan rnd2026 --slurm --dry-run     # look first
python3 scripts/submit_condor_p2.py --scan rnd2026 --slurm \
    --outdir $P2_DATA/my_rnd2026 --jobdir $P2_WORK/jobs/my_rnd2026
squeue -u $USER
```

For Garfield sweeps, write one `run_transparency.sh` / `run_gain.sh` command per line into a
file and submit it with `scripts/slurm_array.sh` (usage in its header). The arguments for
each study are in the old HTCondor files `scripts/garfield/*.sub`.

## 2. Your laptop (optional)

You need conda, mamba or micromamba (micromamba is one binary and needs no admin rights).

```bash
git clone https://github.com/akallitss/P2_Geant_Garfield.git && cd P2_Geant_Garfield
conda env create -f environment.yml          # Geant4 11.2.2, ROOT 6.32, compilers, Python
conda activate p2sim
bash scripts/setup_garfield.sh --build       # ~15 min, once, -> ~/garfield/install
source scripts/setup_local.sh                # every new shell
bash scripts/build.sh                        # must print "--heed enabled"
```

conda-forge has Geant4 11.2.2 rather than 11.2.1. The check in §4 shows this makes no difference.
`--local N` replaces `--slurm` and runs the jobs on N local processes. `--only <tag substring>`
and `--scale 0.1` pick and shrink jobs. The whole 106-job `rnd2026` batch is ~20 CPU-hours.

## 3. Data

**On Feynman the data is already there**, unpacked and checksum-verified (2026-11-01), readable by
the `idedipdetecs` group:

```bash
export P2_DATA=/feynman/work/projets/idedipdetecs/p2_data    # put it in ~/.bashrc, before setup_feynman.sh
ls $P2_DATA/_bundle      # MANIFEST.txt and SHA256SUMS of what was transferred
```

Treat it as **read-only reference data**: it is shared with the group. Write your own outputs to
your work area (e.g. `--outdir $P2_WORK/out/...`), not into `$P2_DATA`. If you can't read it, ask to
be added to the `idedipdetecs` group.

**On a laptop**, copy what you need from Feynman (`rsync -a feynman:$P2_DATA/<dir> ~/p2_data/`).
The bundle was made on lxplus with `scripts/make_data_bundle.sh` (EOS data, plus
`sps_products.tar` from `scripts/make_sps_products_bundle.py`). Extract tarballs with
`for f in *.tar; do tar -xf $f -C $P2_DATA; done` after `sha256sum -c SHA256SUMS`.

`$P2_DATA` has the layout of `/eos/user/a/akallits/p2_sim`, so any `/eos/user/a/akallits/p2_sim/X`
in the docs is `$P2_DATA/X`.

| tarball | size | what |
|---|---|---|
| gas_tables, gas_tables_contam | 2 MB | Magboltz `.gas` files |
| rnd2026_heed | 2.4 GB | production Geant4 sample |
| cell_gap150 | 6.4 GB | COMSOL field maps of the final mesh cell, 8 drift fields |
| transparency | 0.2 GB | all Garfield transparency outputs |
| condor_transp | 45 MB | geometry-run logs (`TRANSPARENCY`/`GAIN` lines), the input to your first task |
| sps_products | 1 MB | SPS Landau MPV vs HV per scan and detector, spark rate vs mesh HV (`spark_qa_*.csv`), corrected `run_config.json` per run (credentials stripped), the Fe55 gain table |
| optional | ~19 GB | validation cell (`cell_nik`), cosmic cell, gap-scan cells, legacy sample, COMSOL model |

The SPS and cosmic raw data are not included. The simulation work only needs the products in
`sps_products.tar`, which the analysis scripts find through `scripts/analysis/p2paths.py`.

## 4. Check your setup before changing anything

```bash
python3 -m stage_b.selftest                       # all pass
python3 vmm/vmm_emulator.py                       # smoke test OK

# Geant4 + Heed: rerun one production point, compare with the lxplus file
python3 scripts/submit_condor_p2.py --scan rnd2026 --local 1 \
    --only ArIso_3p964mm_electron100MeV_th0 \
    --outdir $P2_OUT/check --jobdir $P2_OUT/check_jobs
python3 scripts/check_local_vs_reference.py \
    $P2_OUT/check/p2_ArIso_3p964mm_electron100MeV_th0_t0.root \
    $P2_DATA/rnd2026_heed/p2_ArIso_3p964mm_electron100MeV_th0_t0.root   # -> "same physics"

# Garfield: Ar/iso 95/5 at 625 V/cm, published 0.943
bash scripts/garfield/run_transparency.sh transp \
    $P2_DATA/comsol/p2_sps/h2mm_sweep_gap150/E625 \
    $P2_DATA/gas_tables/ArIso.gas $P2_OUT/check/E625_ArIso 1000 1 null
```

(On Feynman, run the Geant4 check with `--slurm` instead of `--local 1`, and the Garfield check as a
SLURM job with `--mem-per-cpu=4000` (`HOWTO.md` §3). It needs ~3 GB and takes about a minute.)

| check | reference (lxplus) | laptop (2026-11-01) | Feynman (2026-11-01, SLURM) |
|---|---|---|---|
| Geant4 + Heed, 10k e⁻ | `rnd2026_heed` file | same geometry and physics; P(n<10) 1.96 vs 1.89 %, KS p = 0.13 | same geometry and physics; P(n<10) 1.97 vs 1.89 %, KS p = 0.49 |
| transparency, Ar/iso 625 V/cm | 0.943 published (3000 e⁻); **0.934 ± 0.002** rerun with 20 seeds × 1000 | 0.944 ± 0.007 (1000 e⁻) | **0.9345 ± 0.0017** (20 seeds × 1000); same seeds on lxplus differ by −0.0005, p = 0.88 |

## 5. Pitfalls when running and reading outputs

- **Use `--physics heed`** (the default) for anything about primaries or efficiency, and
  **`null`** in the Garfield transparency flags. See `../CLAUDE.md`.
- **One output file per thread.** `-t 4` gives four files. Use `hadd`, or `collect_results.py`.
  `--heed` only works single-threaded (`-t 1`).
- **Never divide by `GetEntries()`.** With `--skip-empty` (always on for photons, which are
  ≥ 99.9 % empty) the trees hold only events with a cluster. The denominator is
  `RunMeta.thrown`, summed over files (`OUTPUT_FORMAT.md` §0).
- **Don't mix files with different `RunMeta.geometryHash`**, and check `gitDirty` = 0 on
  production files.
- **Pad-level observables need `--beam-spread 11.43`** (one ring pitch, the default in
  `submit_condor_p2.py`). A pencil beam sits at one spot on the pad pattern and biases them.
- **Job size:** each job pays ~6.5 s startup. Use ≥ 10⁶ photons or ≥ 10⁴ electrons per job.
  Rates are ~1000–2000 γ/s and ~550 e⁻/s per core.
- **Don't edit a Garfield `.cc` while jobs that use it are queued.** `run_*.sh` recompile
  when the source is newer than the binary.
- If you change the physics list, check that the Cu-K (8.05 keV) fluorescence is still produced.
  It is the main wall background, and a production cut used to suppress it.
- Geometry changes: keep `include/P2Wedge.hh`, `include/SimConfig.hh`,
  `scripts/model/p2_model.py` and `P2_MODEL.md` in sync. Every number needs a source
  (`../CLAUDE.md`).
- Conventions: x/y in the gerber frame (beam axis at 0,0), z = 0 at the front window, +z
  downstream; edep in eV, `-e` and ke in MeV.

## 6. What you cannot do here

- **New COMSOL field solves.** They need a COMSOL licence, and the CERN one only works on lxplus.
  The existing maps cover the final cell at 8 drift fields, plus the variants in the optional
  bundle. For a new geometry, use a Saclay COMSOL licence (ask the colleague who built the
  model; it is `comsol_model.tar`, with parameters explained in `../scripts/comsol/README.md`),
  or ask someone with lxplus.
- **The SPS raw data** (on CERN EOS, 137 GB per scan). Ask for extracted products.

## 7. lxplus (only if you get a CERN account)

`source scripts/setup_lxplus.sh && bash scripts/build.sh`. HTCondor is the default of
`submit_condor_p2.py` (no `--slurm`/`--local`), and the Garfield `.sub`/`.dag` files in
`scripts/garfield/` and `scripts/comsol/` are the batch recipes. The old lxplus notes are in
`archive/RUNBOOK.md`.

A single 1000-electron transparency run has a statistical error of ±0.008, so expect
0.92–0.95 for the Garfield check; for anything finer, run several seeds. The published
0.943 was itself a 2σ-high single run (`OPEN_PROBLEMS.md` §2, item 7).
