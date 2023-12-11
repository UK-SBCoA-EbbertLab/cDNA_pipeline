#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/quick_Test/*.fastq" \
    --trim_dRNA "True" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --out_dir "./test_discovery_dRNA/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --mapq "10" \
    --step "2" \
    --is_dRNA "True" \
    --contamination_ref "../../references/master_contaminant_reference.fasta" \
    --is_chm13 "False" -resume
