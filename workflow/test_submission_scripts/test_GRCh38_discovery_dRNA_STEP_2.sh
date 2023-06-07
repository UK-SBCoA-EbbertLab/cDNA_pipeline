#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/2023-06-06_brain_directRNA_intronic_reads/*.fastq" \
    --ont_reads_txt "/scratch/bag222/data/ont_data/2023-06-06_brain_directRNA_intronic_reads/*.txt" \
    --ref "../../../../cDNA_pipeline/references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../../../cDNA_pipeline/references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../../../cDNA_pipeline/references/hg38.HouseKeepingGenes.bed" \
    --cdna_kit "PCS111" \
    --out_dir "./test_discovery/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --is_dRNA "True" \
    --mapq "0" \
    --step "2" \
    --is_chm13 "False" -resume
