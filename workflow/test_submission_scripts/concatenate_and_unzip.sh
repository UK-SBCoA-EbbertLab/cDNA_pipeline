#!/usr/bin/env bash


nextflow ../main.nf --path "/scratch/bag222/data/mini_test_unzipping/" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --cdna_kit "PCS114" \
    --out_dir "./test_discovery_PCS114/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --mapq "10" \
    --step "2" \
    --is_chm13 "False" -resume
