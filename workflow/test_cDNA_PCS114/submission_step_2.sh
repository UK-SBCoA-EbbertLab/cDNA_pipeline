#!/bin/bash

nextflow ../main.nf --step 2 --path "/scratch/bag222/data/ont_data/test_pychopper/" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./PCS114/" \
    --cdna_kit "PCS114" \
    --is_chm13 "False" \
    --track_reads "True" \
    --mapq "10" \
    --quality_score "9" \
    --prefix "test_PCS114" \
    --is_dRNA "False" \
    --trim_dRNA "False" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" -resume
