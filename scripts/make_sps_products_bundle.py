#!/usr/bin/env python3
"""make_sps_products_bundle.py — pack the SPS / cosmic-bench data products that
the analysis scripts read (≈1 MB) into sps_products.tar for the data bundle.

Run on a machine that has the products (Alexandra's laptop):

    python3 scripts/make_sps_products_bundle.py \\
        --data /local/home/ak271430/Documents/PostDocSaclay/data -o sps_products.tar

then put sps_products.tar next to the other bundle tarballs
(/eos/user/a/akallits/p2_sim/handover_bundle) and re-run
`sha256sum *.tar > SHA256SUMS` there.

Contents (paths inside the tar, relative to $P2_DATA):
  sps_products/analysis/<P2_IN|P2_MID|P2_OUT>/<scan>/.../beam_mpv_vs_hv_spark_vetoed.csv
  sps_products/analysis/<det>/<scan>/.../26_hv_spark_qa/spark_qa_<det>.csv  (spark rate vs mesh HV)
  sps_products/runs/<run>/run_config.json   (corrected gas + gas_note; HV per sub-run)
  cosmic_products/det2/.../fe55_gain_vs_hv_spark_vetoed.csv

Every "username"/"password" field is removed from the run configs (they hold
the HV controller's login), so the bundle carries no credentials.
"""
import argparse
import io
import json
import tarfile
from pathlib import Path

SECRET_KEYS = {"username", "password", "passwd", "token", "api_key", "apikey", "secret"}


def scrub(obj):
    if isinstance(obj, dict):
        return {k: scrub(v) for k, v in obj.items() if k.lower() not in SECRET_KEYS}
    if isinstance(obj, list):
        return [scrub(v) for v in obj]
    return obj


def add_bytes(tar, name, data):
    info = tarfile.TarInfo(name)
    info.size = len(data)
    tar.addfile(info, io.BytesIO(data))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", required=True,
                   help="folder holding SPS_Beam_Test/ and Cosmic_Bench/")
    p.add_argument("-o", "--out", default="sps_products.tar")
    a = p.parse_args()
    data = Path(a.data)
    ana = data / "SPS_Beam_Test/mpgd26_workspace/products/analysis"
    runs = data / "SPS_Beam_Test/mpgd26_workspace/eos_inventory/runs"
    fe55 = data / "Cosmic_Bench/Analysis/det2"

    n = {"csv": 0, "spark_qa": 0, "run_config": 0, "fe55": 0}
    with tarfile.open(a.out, "w") as tar:
        for f in sorted(ana.rglob("beam_mpv_vs_hv_spark_vetoed.csv")):
            tar.add(f, arcname=f"sps_products/analysis/{f.relative_to(ana)}")
            n["csv"] += 1
        # stage 26_hv_spark_qa: spark rate vs mesh HV per sub-run (HV-headroom curves)
        for f in sorted(ana.rglob("spark_qa_*.csv")):
            tar.add(f, arcname=f"sps_products/analysis/{f.relative_to(ana)}")
            n["spark_qa"] += 1
        for f in sorted(runs.glob("*/run_config.json")):
            clean = json.dumps(scrub(json.load(open(f))), indent=1).encode()
            assert b'"password"' not in clean.lower()
            add_bytes(tar, f"sps_products/runs/{f.parent.name}/run_config.json", clean)
            n["run_config"] += 1
        for f in sorted(fe55.rglob("fe55_gain_vs_hv_spark_vetoed.csv")):
            tar.add(f, arcname=f"cosmic_products/det2/{f.relative_to(fe55)}")
            n["fe55"] += 1
    print(f"wrote {a.out}: {n}")


if __name__ == "__main__":
    main()
