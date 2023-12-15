#!/usr/bin/env bash


nextflow ../main.nf --path "/scratch/bag222/data/ont_data/test_dRNA_data/basecalling/" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --out_dir "./test_dRNA_NOtrim/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --step "2" \
    --is_dRNA "True" \
    --contamination_ref "../../references/master_contaminant_reference.fasta" \
    --is_chm13 "False" -resume
