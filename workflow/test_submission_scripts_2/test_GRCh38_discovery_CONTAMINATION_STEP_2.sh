#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/test_data/*.fastq" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --cdna_kit "PCS111" \
    --out_dir "./test_discovery_PCS111/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --mapq "10" \
    --step "2" \
    --contamination_ref "../../references/master_contaminant_reference.fasta" \
    --is_chm13 "False" -resume
