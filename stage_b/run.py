#!/usr/bin/env python3
"""run.py — Stage B entry point.

    python3 -m stage_b.run <stageA.root> -o padhits.root

Reads a Stage A file, digitizes its gas clusters, writes a `PadHitTree`
(one row per (event, pad) that collected charge) plus a `StageBMeta` tree
carrying every parameter the run was made under.

The metadata is not optional decoration. Stage B has four placeholders in it
(gas transport, mesh transparency, induced current, pillars), so a number
produced here is only interpretable alongside the assumptions that produced
it. StageBMeta records them, and mirrors Stage A's geometry and git hashes so
a Stage B product can always be traced back to the Stage A file it came from.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# Allow both `python3 -m stage_b.run` and `python3 stage_b/run.py`.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from stage_b.clusters import read_clusters
    from stage_b.digitize import DigitizerConfig, digitize
    from stage_b.frame import read_frame
    from stage_b.gas import (JSONTable, MagboltzTable, PlaceholderTable,
                              SPS_DRIFT_FIELD_V_PER_CM)
    from stage_b.padplane import PadPlane
else:
    from .clusters import read_clusters
    from .digitize import DigitizerConfig, digitize
    from .frame import read_frame
    from .gas import (JSONTable, MagboltzTable, PlaceholderTable,
                      SPS_DRIFT_FIELD_V_PER_CM)
    from .padplane import PadPlane


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("infile", help="Stage A ROOT file (one worker file, or hadd-merged)")
    p.add_argument("-o", "--outfile", default=None,
                   help="Output ROOT file (default: <infile>_padhits.root)")
    p.add_argument("--gain", type=float, default=1.0e4)
    p.add_argument("--polya-theta", type=float, default=2.0)
    p.add_argument("--mesh-transparency", type=float, default=1.0,
                   help="Constant eps (P0.13: a field map would be better)")
    p.add_argument("--transport", default="magboltz",
                   choices=["magboltz", "placeholder"],
                   help="Gas transport: Magboltz tables (default, dry gas, "
                        "stage_b/gas_tables/) or the old PLACEHOLDER values")
    p.add_argument("--drift-field", type=float, default=SPS_DRIFT_FIELD_V_PER_CM,
                   help="Drift field [V/cm] at which the transport is taken "
                        f"(default {SPS_DRIFT_FIELD_V_PER_CM:g}: SPS July 2026, "
                        "drift 700 V - mesh 450 V over 4 mm)")
    p.add_argument("--gas-table", default=None,
                   help="A hand-made JSON transport table (JSONTable); "
                        "overrides --transport")
    p.add_argument("--v-scale", type=float, default=1.0,
                   help="Multiply drift velocity, to state a wet-gas "
                        "assumption explicitly (MX17: water dominates v_d)")
    p.add_argument("--mapping-revision", default="B", choices=["A", "B"])
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--max-events", type=int, default=None)
    p.add_argument("--quiet", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    infile = Path(args.infile)
    outfile = Path(args.outfile) if args.outfile else \
        infile.with_name(infile.stem + "_padhits.root")

    frame = read_frame(infile)
    if args.gas_table:
        table = JSONTable(args.gas_table, args.v_scale)
    elif args.transport == "magboltz":
        table = MagboltzTable(v_scale=args.v_scale)
    else:
        table = PlaceholderTable(args.v_scale)
    transport = table.for_gas(frame.gas, args.drift_field)
    padplane = PadPlane.load(args.mapping_revision)
    cfg = DigitizerConfig(gain=args.gain, polya_theta=args.polya_theta,
                          mesh_transparency=args.mesh_transparency,
                          mapping_revision=args.mapping_revision,
                          seed=args.seed)

    clusters = read_clusters(infile, args.max_events)
    hits = digitize(clusters, frame, transport, padplane, cfg)

    if not args.quiet:
        print(f"  Stage A       : {infile.name}")
        print(f"  Gas           : {frame.gas} ({frame.gas_label})")
        print(f"  Drift gap     : {frame.drift_gap_mm} mm   "
              f"amp {frame.amp_gap_um} um")
        print(f"  Transport     : v_d = {transport.v_drift_mm_per_ns*1000:.1f} um/ns"
              f"   [{transport.provenance}]")
        print(f"  Events        : {hits.n_events}   thrown {frame.thrown}")
        print(f"  Electrons     : {hits.n_electrons_total}"
              f"   attached {hits.n_attached}"
              f"   mesh-lost {hits.n_mesh_lost}"
              f"   off-plane {hits.n_offplane}")
        print(f"  Pad hits      : {len(hits.pad)}")
        if hits.n_events:
            per_ev = len(hits.pad) / hits.n_events
            print(f"  Pads / event  : {per_ev:.3f}")
        if frame.beam_spread_mm <= 0:
            print("  NOTE: Stage A used a PENCIL BEAM (beamSpread_mm = 0). "
                  "Pad multiplicity and charge sharing from this file reflect "
                  "one impact point, not the pad cell.")
        if frame.git_dirty:
            print("  WARNING: Stage A file came from a DIRTY working tree.")
        if hits.n_electrons_total:
            off = hits.n_offplane / hits.n_electrons_total
            if off > cfg.warn_offplane_frac:
                print(f"  WARNING: {off:.1%} of electrons landed outside the "
                      f"instrumented pad area.")

    write_output(outfile, hits, frame, transport, cfg, args.drift_field)
    if not args.quiet:
        print(f"  Wrote         : {outfile}")
    return hits


def write_output(outfile: Path, hits, frame, transport, cfg, drift_field):
    import uproot as up
    with up.recreate(outfile) as f:
        f["PadHitTree"] = {
            "eventID":   hits.event.astype(np.int32),
            "pad":       hits.pad.astype(np.int32),
            "charge_e":  hits.charge_e.astype(np.float64),
            "t_first_ns": hits.t_first_ns.astype(np.float64),
            "nElectrons": hits.n_electrons.astype(np.int32),
        }
        prov = dict(hits.provenance)
        prov_json = json.dumps(prov, default=str)
        f["StageBMeta"] = {
            "gas":            np.array([frame.gas]),
            "geometryHash":   np.array([frame.geometry_hash]),
            "gitHash":        np.array([frame.git_hash]),
            "gitDirty":       np.array([int(frame.git_dirty)], np.int32),
            "thrown":         np.array([frame.thrown], np.int64),
            "nEvents":        np.array([hits.n_events], np.int64),
            "driftGap_mm":    np.array([frame.drift_gap_mm]),
            "ampGap_um":      np.array([frame.amp_gap_um]),
            "beamSpread_mm":  np.array([frame.beam_spread_mm]),
            "gain":           np.array([cfg.gain]),
            "polyaTheta":     np.array([cfg.polya_theta]),
            "meshTransparency": np.array([cfg.mesh_transparency]),
            "vDrift_mm_per_ns": np.array([transport.v_drift_mm_per_ns]),
            "sigmaL_mm_per_sqrt_mm": np.array([transport.sigma_L_mm_per_sqrt_mm]),
            "sigmaT_mm_per_sqrt_mm": np.array([transport.sigma_T_mm_per_sqrt_mm]),
            "attachLength_mm": np.array([transport.attach_length_mm]),
            "driftField_V_per_cm": np.array([float(drift_field)]),
            "vScale":         np.array([transport.v_scale]),
            "mappingRevision": np.array([cfg.mapping_revision]),
            "seed":           np.array([cfg.seed], np.int64),
            "nElectronsTotal": np.array([hits.n_electrons_total], np.int64),
            "nAttached":      np.array([hits.n_attached], np.int64),
            "nMeshLost":      np.array([hits.n_mesh_lost], np.int64),
            "nOffPlane":      np.array([hits.n_offplane], np.int64),
            "transportProvenance": np.array([transport.provenance]),
            "provenanceJSON": np.array([prov_json]),
        }


if __name__ == "__main__":
    main()
