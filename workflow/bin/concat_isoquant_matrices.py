#!/usr/bin/env python3
import polars as pl
import argparse
import glob
import os
import sys

# ---------- GTF handling functions ----------------------

def extract_tx_gene(df: pl.DataFrame) -> pl.DataFrame:

    pairs = df.select([
        pl.col("attributes").str.extract(r'transcript_id "([^"]+)"', 1).alias("transcript_id"),
        pl.col("attributes").str.extract(r'gene_id "([^"]+)"', 1).alias("gene_id")
        ]).drop_nulls(subset=["transcript_id", "gene_id"]).unique()

    return pairs


def build_tx_to_gene_map(gtf_path: str) -> pl.DataFrame:
    if not gtf_path:
        raise SystemExit("No GTF file provided.")

    gtf = pl.read_csv(gtf_path,
            separator="\t",
            has_header=False,
            infer_schema_length=100,
            new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"],
            schema_overrides={"seqname" : str},
            comment_prefix="#",
            ignore_errors=True)

    tx_gene = extract_tx_gene(gtf)
    print(f"[map] unique transcript_id <-> gene_id pairs: {tx_gene.height}")
    print(tx_gene)

    conflicts = (
            tx_gene.group_by("transcript_id").agg(pl.col("gene_id").n_unique().alias("n_genes"))
                .filter(pl.col("n_genes") > 1)
                .select("transcript_id", "n_genes")
                )

    n_conf = conflicts.height
    if n_conf > 0:
        print(f"\n[ERROR] Detected {n_conf} isoform_id(s) mapping to multiple gene_id values across provided GTFs.")

        # Build a long-form report of conflicting mappings
        detail = (
                tx_gene.join(conflicts.select("transcript_id"), on="transcript_id", how="inner")
                    .select(["transcript_id", "gene_id"])
                    .sort(["transcript_id", "gene_id"])
                    )
        # Show a small preview
        print("[ERROR] Example conflicts (first 20 rows):")
        print(detail.head(20).to_pandas().to_string(index=False))

        sys.exit(1)

    return tx_gene

# ------------- helper functions --------------------------------------

def strip_before_first_underscore(colname: str) -> str:
    """
    Remove everything up to and including the first underscore. 
    If no underscore exists, return unchanged. 
    """
    if "_" in colname:
        return colname.split("_", 1)[1]
    return colname

def fix_sample_headers(df: pl.DataFrame, id_cols: list[str]) -> pl.DataFrame:
    ren = {}
    for c in df.columns:
        if c in id_cols:
            continue
        new = strip_before_first_underscore(c)
        if new != c:
            ren[c] = new
    return df.rename(ren) if ren else df

def concat_counts(dfs: list[pl.DataFrame]) -> pl.DataFrame:
    if not all(df.columns == dfs[0].columns for df in dfs):
        print("Not all matrices have the same columns")
        sys.exit(1)

    df = pl.concat(dfs)

    dupes = df.filter(df["gene_id"].is_duplicated())
    if dupes.height > 0:
        print("Duplicate gene_id entries found:")
        print(dupes)
        sys.exit(1)

    return df

# ----------------------------------------------------------------

def read_counts_matrix(df_path: str) -> pl.DataFrame:
    df = pl.read_csv(df_path, 
            separator="\t",
            infer_schema_length=100,
            comment_prefix="#")

    df = fix_sample_headers(df, ["gene_id"])
    return df

# ----------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
            description="Concat IsoQuant counts using Polars. Works for transcript and gene counts"
            )

    ap.add_argument("--counts_matrices", action="extend", nargs="+", help="List of counts matrices", required=True)
    ap.add_argument("--counts_type", help="String informing what type of counts we recieved (gene or transcript)", required=True)
    ap.add_argument("--out", help="Name for output combined TSV", required=True)
    ap.add_argument("--gtf", help="GTF path")
    args = ap.parse_args()

    dfs = [read_counts_matrix(p) for p in args.counts_matrices]
        
    combined_counts_matrix = concat_counts(dfs)

    if args.counts_type == "transcript":
        gtf_key = build_tx_to_gene_map(args.gtf)

        combined_counts_matrix = combined_counts_matrix.rename({ "gene_id" : "transcript_id" })
        combined_counts_matrix = gtf_key.join(combined_counts_matrix, on="transcript_id", how="right")

        combined_counts_matrix = combined_counts_matrix.select(["transcript_id", "gene_id", pl.all().exclude("transcript_id", "gene_id")])

    combined_counts_matrix.write_csv(args.out, separator='\t', quote_style='never')
    return


if __name__ == "__main__":
    main()
