#!/usr/bin/env python3
import polars as pl
import glob
import os

# read in all the files
file_paths = glob.glob("*.abridged.tsv")

dataframes = []

for file_path in file_paths:
    try:
        filename = os.path.basename(file_path)
        sample = filename.replace("_ctat-LR-fusion.fusion_predictions.abridged.tsv", "")
        df = pl.read_csv(file_path, separator="\t")
        df = df.select("#FusionName", "LeftGene", "LeftLocalBreakpoint", "LeftBreakpoint", "RightGene", "RightLocalBreakpoint", "RightBreakpoint", "SpliceType", "annots", "num_LR").rename({"num_LR": sample})
        dataframes.append(df)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

df_combined = dataframes.pop()

while dataframes:
    df_combined = df_combined.join(dataframes.pop(), on=("#FusionName", "LeftGene", "LeftLocalBreakpoint", "LeftBreakpoint", "RightGene", "RightLocalBreakpoint", "RightBreakpoint", "SpliceType", "annots"), how="full", coalesce=True)

df_combined.write_csv("combined_ctat-LR-fusion.fusion_predictions.abridged.tsv", separator="\t", quote_style='never')
