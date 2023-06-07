#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/test_data/*.fastq" \
    --ref "../../../../cDNA_pipeline/references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../../../cDNA_pipeline/references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../../../cDNA_pipeline/references/hg38.HouseKeepingGenes.bed" \
    --cdna_kit "PCS111" \
    --out_dir "./test_discovery_unstranded/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --mapq "10" \
    --step "2" \
    --is_chm13 "False" -resume
