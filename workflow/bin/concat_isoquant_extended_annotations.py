#!/usr/bin/env python3
import polars as pl
import argparse
import glob
import os
import sys

# ---------- GTF handling functions ----------------------

def read_gtf(path: str) -> pl.DataFrame:
    """
    Read gtf, comment lines skipped
    """

    return pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            infer_schema_length=100,
            new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"],
            schema_overrides={"seqname" : str},
            comment_prefix="#",
            ignore_errors=True
        )

def build_tx_to_gene_map(annotation: str,gtf_paths: list[str], out_gtf_path: str):
    if not gtf_paths:
        raise SystemExit("No GTF files provided.")

    gtf = read_gtf(annotation)

    for p in gtf_paths:
        print(f"[gtf] parsing: {p}")
        tmp = read_gtf(p).filter(pl.col("source") == "IsoQuant")
        gtf = pl.concat([gtf, tmp], how="vertical_relaxed")


    merged = gtf.unique()
    print(f"[gtf] unique rows after global unique: {merged.height}")

    merged = merged.sort(by=["seqname", "start", "end", "feature", "strand", "source", "frame"])
    print(f"[gtf] sorted")

    GTF_COLS = ["seqname","source","feature","start","end","score","strand","frame","attributes"]
    merged.select(GTF_COLS).write_csv(out_gtf_path, separator="\t", include_header=False, quote_style='never')
    print(f"[gtf] wrote merged GTF -> {out_gtf_path}")

    return


def main():
    ap = argparse.ArgumentParser(
            description="Concat IsoQuant extended annotations using Polars."
            )

    ap.add_argument("--out_gtf", help="Name for output combined GTF", required=True)
    ap.add_argument("--annotation", help="The reference annotation", required=True)
    ap.add_argument("--gtfs", action="extend", nargs="+", help="GTF files", required=True)
    args = ap.parse_args()

    gtf_paths = []
    gtf_paths.extend(sorted(args.gtfs))
    build_tx_to_gene_map(args.annotation, gtf_paths, args.out_gtf)
    return


if __name__ == "__main__":
    main()
