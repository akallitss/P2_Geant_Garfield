#!/usr/bin/env python3
"""
submit_condor_p2.py — HTCondor submission for the P2 wedge campaign.

Why this is a separate script rather than a flag on submit_condor.py: that
one is the inherited vacuum/Al-shielding scan and speaks the MX17 modes. It
passes no `-m`, and because mm_sim's default mode is `p2` it was silently
running P2-wedge geometry at default gaps under vacuum-mode tags. It now
passes `-m vacuum` explicitly; P2 campaign runs come through here.

What this script knows that the old one did not:
  * `-m p2`, and every P2 geometry flag (--drift-gap, --amp-gap, --back-gap,
    --bulge-*) so a gap scan is expressible;
  * the beam: --gun-x/y, --gun-theta/--gun-phi for the angle scan,
    and --beam-spread, which is REQUIRED for pad-level observables;
  * --skip-empty, without which a 10^7-photon point writes 10^7 rows of
    which ~10^4 are non-empty;
  * the campaign naming scheme (SIM_CAMPAIGN_PLAN §9);
  * a manifest row per run point, written at submit time.

The manifest is belt-and-braces. Since 2026-08-07 every output file carries
its own RunMeta tree (thrown count, git hash, geometry hash, full beam and
geometry configuration), so the manifest is no longer the only record of what
a file is -- collect_results.py reads RunMeta and does not need this file.
Keep writing it anyway: it is the record of what was *intended*, and
diffing intent against RunMeta is how you notice that a third of a campaign
died on the worker nodes.

Usage:
    python3 scripts/submit_condor_p2.py --dry-run
    python3 scripts/submit_condor_p2.py --scan photon --outdir /eos/user/d/dneff/p2
    python3 scripts/submit_condor_p2.py --scan angle --gases ArIso NeIso9010
    python3 scripts/submit_condor_p2.py --scan rnd2026 --outdir /eos/user/a/akallits/p2_sim/rnd2026
"""

import argparse
import csv
import json
import os
import random
import stat
import subprocess
import sys
import textwrap
from pathlib import Path

# ============================================================
# CAMPAIGN GRID — SIM_CAMPAIGN_PLAN §4
# ============================================================

# The §7.4 campaign gases. `mm_sim --list-gases` is the authority; these must
# be a subset of it, and the script checks that before submitting.
CAMPAIGN_GASES = [
    "ArIso", "ArCO2Iso9352",
    "NeIso9010", "NeIso8020", "NeCO2Iso9532", "NeIsoCO29055",
    "NeCH4937", "NeC2H69010",
]

# §4.1 electron scan [MeV]; dense around the elastic lines.
ELECTRON_ENERGIES = [1, 2, 5, 10, 20, 30, 50, 75, 100, 119, 130, 144, 155,
                     175, 200]

# §4.2 photon scan [MeV]. 3-10 keV covers the fluorescence/K-edge structure
# (Ar K 3.2, Cu K 8.05, Fe K 6.4); 50-150 keV is the dominant background band.
PHOTON_ENERGIES = [3e-3, 5e-3, 8e-3, 10e-3, 15e-3, 20e-3, 30e-3, 40e-3,
                   50e-3, 60e-3, 80e-3, 100e-3, 120e-3, 150e-3, 200e-3,
                   300e-3, 500e-3]

# §4.1 angles [deg from the wedge normal]. Full set at the baseline energies.
ANGLES_FULL   = [0, 10, 20, 30, 40]
ANGLES_COARSE = [0, 10]
ANGLE_FULL_ENERGIES = [30, 119, 155]

# §6 gap scan [mm].
DRIFT_GAPS = [1.0, 2.0, 3.0, 4.0]

# §4.3 test-beam anchor: 150 GeV muons, SPS July 2026 (Alexandra 2026-10-06).
MUON_ENERGY_MEV = 150000.0

# ── Final R&D phase, October 2026 (`--scan rnd2026`) ─────────────────────────
# Gain recovery: larger isobutane fraction and larger drift gap, Ar vs Ne,
# plus the e/gamma question. 3.964 mm = drift gas of the 4 mm frame
# (P2_Frame_V2.0.stp on the Dynamask); 5 mm has no frame drawing yet, the
# frame is stretched by --drift-gap.
RND_GASES = ["ArIso", "ArIso9010", "NeIso", "NeIso8515"]
RND_DRIFT_GAPS = [3.964, 5.0]
RND_PHOTON_ENERGIES = [10e-3, 15e-3, 20e-3, 30e-3, 40e-3, 50e-3, 60e-3,
                       80e-3, 100e-3, 150e-3]
RND_ELECTRON_MEV = 100.0             # signal electrons (Alexandra 2026-10-06)
RND_ELECTRON_ANGLES = [0, 20, 40]
RND_TB_GASES = ["ArCO2Iso9352", "ArCF4Iso"]   # SPS July 2026, 4 mm frame

# Statistics. Photon points need far more throws because P(interaction) is
# ~10^-3 and the rarest reported quantity needs >=100 events (§9).
NEVENTS = {
    "electron": 10_000,
    "gamma":    1_000_000,
    "muon":     100_000,
}
NEVENTS_BASELINE_BOOST = 10      # at the baseline energies

# One ring pitch. Below this the impact point does not average over the pad
# cell and pad-level observables are biased -- see SimConfig.hh and
# MX17_Geant/design/RESPONSE_SIM_PLAN.md §7.
RING_PITCH_MM = 11.4290

JOB_FLAVOUR = "workday"          # 8 h; espresso=20min, longlunch=2h


# ============================================================
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scan", default="electron",
                   choices=["electron", "photon", "angle", "gap", "muon", "all",
                            "rnd2026"],
                   help="Which block of the campaign to submit (default: electron)")
    p.add_argument("--gases", nargs="+", default=None,
                   help="Default: the §7.4 list, or RND_GASES for rnd2026")
    p.add_argument("--outdir", default="./p2_results",
                   help="Where ROOT output lands (use EOS for real campaigns)")
    p.add_argument("--jobdir", default="./condor_p2",
                   help="Where wrapper/submit/logs/manifest go")
    p.add_argument("--nevents", type=int, default=None,
                   help="Override events per job")
    p.add_argument("--drift-gap", type=float, default=3.964,
                   help="Drift gas [mm] for non-gap scans (default: 3.964, "
                        "the 4 mm frame)")
    p.add_argument("--beam-spread", type=float, default=RING_PITCH_MM,
                   help=f"Impact-point spread [mm] (default: {RING_PITCH_MM}, "
                        f"one ring pitch). 0 = pencil beam; only use 0 if you "
                        f"are NOT reading a pad-level observable")
    p.add_argument("--exe", default=None)
    p.add_argument("--setup", default=None,
                   help="Setup script to source on the worker "
                        "(default: scripts/setup_lxplus.sh next to this file)")
    p.add_argument("--flavour", default=JOB_FLAVOUR)
    p.add_argument("--dry-run", action="store_true",
                   help="Write the submit file and manifest, do not condor_submit")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def find_exe(explicit=None):
    if explicit:
        return str(Path(explicit).resolve())
    here = Path(__file__).parent
    for c in (here.parent / "build" / "mm_sim",
              here.parent / "install" / "bin" / "mm_sim",
              Path(os.environ.get("MM_SIM_EXE", "/nonexistent"))):
        if c.is_file():
            return str(c.resolve())
    return None


def known_gases(exe):
    """Ask the binary which gases exist, rather than trusting this file.

    Catches the case where the campaign list here has drifted from
    GasMixtures.hh -- otherwise the mismatch surfaces as every job in a block
    exiting non-zero on the worker, hours later.
    """
    try:
        out = subprocess.run([exe, "--list-gases"], capture_output=True,
                             text=True, timeout=60).stdout
    except Exception:
        return None
    names = set()
    for line in out.splitlines():
        parts = line.strip().lstrip("* ").split()
        if parts and not line.startswith("Available"):
            names.add(parts[0])
    return names or None


# ============================================================
def make_tag(gas, particle, energy_mev, gap_mm, theta_deg):
    """Campaign naming, SIM_CAMPAIGN_PLAN §9.

    p2_<gas>_<gap>mm_<part><E><unit>_th<deg>
    """
    if energy_mev >= 1000:
        e_str = f"{energy_mev/1000:.6g}GeV"
    elif energy_mev >= 1:
        e_str = f"{energy_mev:.6g}MeV"
    else:
        e_str = f"{energy_mev*1000:.6g}keV"
    e_str = e_str.replace(".", "p")
    g_str = f"{gap_mm:g}mm".replace(".", "p")
    t_str = f"th{theta_deg:g}".replace(".", "p")
    return f"p2_{gas}_{g_str}_{particle}{e_str}_{t_str}"


def build_jobs(args):
    """Expand the requested scan into (gas, particle, E, gap, theta, N) rows."""
    jobs = []

    def n_for(particle, energy_mev, baseline=False):
        if args.nevents:
            return args.nevents
        n = NEVENTS[particle]
        return n * NEVENTS_BASELINE_BOOST if baseline else n

    if args.scan == "rnd2026":
        for gas in args.gases:
            for gap in RND_DRIFT_GAPS:
                for e in RND_PHOTON_ENERGIES:
                    jobs.append((gas, "gamma", e, gap, 0.0, n_for("gamma", e)))
                for th in RND_ELECTRON_ANGLES:
                    jobs.append((gas, "electron", RND_ELECTRON_MEV, gap, float(th),
                                 n_for("electron", RND_ELECTRON_MEV)))
        for gas in RND_TB_GASES:
            jobs.append((gas, "muon", MUON_ENERGY_MEV, RND_DRIFT_GAPS[0], 0.0,
                         n_for("muon", MUON_ENERGY_MEV)))
        return jobs

    want = ({"electron", "photon", "angle", "gap", "muon"}
            if args.scan == "all" else {args.scan})

    for gas in args.gases:
        if "electron" in want:
            for e in ELECTRON_ENERGIES:
                baseline = e in ANGLE_FULL_ENERGIES
                jobs.append((gas, "electron", float(e), args.drift_gap, 0.0,
                             n_for("electron", e, baseline)))
        if "photon" in want:
            for e in PHOTON_ENERGIES:
                jobs.append((gas, "gamma", float(e), args.drift_gap, 0.0,
                             n_for("gamma", e)))
        if "angle" in want:
            for e in ELECTRON_ENERGIES:
                angles = (ANGLES_FULL if e in ANGLE_FULL_ENERGIES
                          else ANGLES_COARSE)
                for th in angles:
                    if th == 0.0:
                        continue            # already covered by the electron block
                    jobs.append((gas, "electron", float(e), args.drift_gap,
                                 float(th), n_for("electron", e,
                                                  e in ANGLE_FULL_ENERGIES)))
        if "gap" in want:
            for gap in DRIFT_GAPS:
                for e in ANGLE_FULL_ENERGIES:
                    jobs.append((gas, "electron", float(e), gap, 0.0,
                                 n_for("electron", e, True)))
        if "muon" in want:
            jobs.append((gas, "muon", MUON_ENERGY_MEV, args.drift_gap, 0.0,
                         n_for("muon", MUON_ENERGY_MEV)))

    # Deduplicate: the electron and gap blocks overlap at gap == args.drift_gap.
    seen, out = set(), []
    for j in jobs:
        k = j[:5]
        if k not in seen:
            seen.add(k)
            out.append(j)
    return out


# ============================================================
def write_wrapper(job_dir: Path, exe: str, setup: str) -> Path:
    wrapper = job_dir / "run_job_p2.sh"
    wrapper.write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -e
        source "{setup}"

        GAS="$1"; PARTICLE="$2"; ENERGY="$3"; NEVENTS="$4"
        OUTFILE="$5"; SEED="$6"; GAP="$7"; THETA="$8"; SPREAD="$9"
        SKIPEMPTY="${{10}}"

        echo "=== mm_sim p2 job ==="
        echo "  Node      : $(hostname)"
        echo "  Gas       : $GAS"
        echo "  Particle  : $PARTICLE  @ $ENERGY MeV"
        echo "  Events    : $NEVENTS"
        echo "  Drift gap : $GAP mm"
        echo "  Theta     : $THETA deg   Spread: $SPREAD mm"
        echo "  Output    : $OUTFILE"
        echo "====================="

        EXTRA=""
        if [ "$SKIPEMPTY" = "1" ]; then EXTRA="--skip-empty"; fi

        # -m p2 is explicit even though it is the default. The whole reason
        # this script exists is that an implicit mode silently produced the
        # wrong thing in the inherited scripts.
        "{exe}" \\
            -m p2               \\
            -g "$GAS"           \\
            -p "$PARTICLE"      \\
            -e "$ENERGY"        \\
            -n "$NEVENTS"       \\
            -o "$OUTFILE"       \\
            -s "$SEED"          \\
            -t 1                \\
            --drift-gap "$GAP"  \\
            --gun-theta "$THETA" \\
            --beam-spread "$SPREAD" \\
            $EXTRA

        echo "Job done: $(date)"
    """))
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


def write_submit(job_dir: Path, wrapper: Path, jobs, outdir: Path,
                 flavour: str, spread: float, seed0: int) -> Path:
    submit_file = job_dir / "p2_campaign.sub"
    log_dir = job_dir / "logs"
    log_dir.mkdir(exist_ok=True, parents=True)

    lines = [
        f"executable            = {wrapper}",
        f"output                = {log_dir}/$(tag).out",
        f"error                 = {log_dir}/$(tag).err",
        f"log                   = {log_dir}/condor.log",
        f'+JobFlavour           = "{flavour}"',
        "request_cpus          = 1",
        "request_memory        = 2048",
        "request_disk          = 4096",
        "should_transfer_files = YES",
        "when_to_transfer_output = ON_EXIT",
        'requirements          = (OpSysAndVer =?= "AlmaLinux9")',
        "",
        "arguments             = $(gas) $(particle) $(energy) $(nevents) "
        "$(outfile) $(seed) $(gap) $(theta) $(spread) $(skipempty)",
        "",
        "queue gas,particle,energy,nevents,outfile,seed,gap,theta,spread,"
        "skipempty,tag from (",
    ]

    rng = random.Random(seed0)
    rows = []
    for (gas, particle, energy, gap, theta, nev) in jobs:
        tag = make_tag(gas, particle, energy, gap, theta)
        outfile = outdir / tag
        seed = rng.randint(1, 2**31 - 1)
        skipempty = 1 if particle == "gamma" else 0
        lines.append(f"  {gas}, {particle}, {energy:.10g}, {nev}, {outfile}, "
                     f"{seed}, {gap:g}, {theta:g}, {spread:g}, {skipempty}, {tag}")
        rows.append(dict(tag=tag, gas=gas, particle=particle,
                         energy_MeV=energy, nevents=nev, drift_gap_mm=gap,
                         theta_deg=theta, beam_spread_mm=spread,
                         skip_empty=skipempty, seed=seed,
                         outfile=str(outfile)))
    lines.append(")")
    submit_file.write_text("\n".join(lines) + "\n")
    return submit_file, rows


def write_manifest(job_dir: Path, rows, args, exe):
    """Intent, recorded at submit time. Compare against RunMeta afterwards."""
    csv_path = job_dir / "manifest.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    git = "unknown"
    try:
        git = subprocess.run(["git", "rev-parse", "--short=12", "HEAD"],
                             cwd=Path(__file__).parent.parent,
                             capture_output=True, text=True).stdout.strip()
    except Exception:
        pass
    (job_dir / "manifest.json").write_text(json.dumps(
        dict(scan=args.scan, gases=args.gases, exe=exe, git=git,
             drift_gap_mm=args.drift_gap, beam_spread_mm=args.beam_spread,
             n_jobs=len(rows)), indent=2))
    return csv_path


# ============================================================
def main():
    args = parse_args()
    if args.gases is None:
        args.gases = RND_GASES if args.scan == "rnd2026" else CAMPAIGN_GASES

    exe = find_exe(args.exe)
    if not exe:
        sys.exit("ERROR: mm_sim not found. Build it, or pass --exe.")

    setup = args.setup or str((Path(__file__).parent / "setup_lxplus.sh").resolve())
    if not Path(setup).is_file():
        sys.exit(f"ERROR: setup script not found: {setup}")

    avail = known_gases(exe)
    if avail:
        check = args.gases + (RND_TB_GASES if args.scan == "rnd2026" else [])
        bad = [g for g in check if g not in avail]
        if bad:
            sys.exit(f"ERROR: unknown gas(es) {bad}.\n"
                     f"Run `{exe} --list-gases` for the current table.")

    jobs = build_jobs(args)
    if not jobs:
        sys.exit("ERROR: the requested scan expanded to zero jobs.")

    if args.beam_spread <= 0:
        print("WARNING: --beam-spread 0 is a pencil beam. Pad multiplicity, "
              "charge sharing and any positional observable will be biased by "
              "wherever the aim point happens to sit relative to the pad "
              "boundaries. Only do this deliberately.")

    outdir = Path(args.outdir).resolve()
    job_dir = Path(args.jobdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    job_dir.mkdir(parents=True, exist_ok=True)

    wrapper = write_wrapper(job_dir, exe, setup)
    submit_file, rows = write_submit(job_dir, wrapper, jobs, outdir,
                                     args.flavour, args.beam_spread, args.seed)
    manifest = write_manifest(job_dir, rows, args, exe)

    total_events = sum(r["nevents"] for r in rows)
    print(f"  Scan         : {args.scan}")
    print(f"  Gases        : {len(args.gases)}  {args.gases}")
    print(f"  Jobs         : {len(rows)}")
    print(f"  Total events : {total_events:,}")
    print(f"  Beam spread  : {args.beam_spread} mm")
    print(f"  Executable   : {exe}")
    print(f"  Output       : {outdir}")
    print(f"  Submit file  : {submit_file}")
    print(f"  Manifest     : {manifest}")

    if args.dry_run:
        print("\n  --dry-run: not submitting. First few jobs:")
        for r in rows[:10]:
            print(f"    {r['tag']}  ({r['nevents']:,} events)")
        if len(rows) > 10:
            print(f"    ... and {len(rows)-10} more")
        return

    subprocess.run(["condor_submit", str(submit_file)], check=True)
    print(f"\n  Submitted. Then:  python3 scripts/collect_results.py --indir {outdir}")


if __name__ == "__main__":
    main()
