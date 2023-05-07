#!/usr/bin/env bash


nextflow ../main.nf --bam "./results/test_discovery/mapping_cDNA/*bam" \
    --bai "./results/test_discovery/mapping_cDNA/*bai" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --out_dir "./test_discovery/" \
    --bambu_track_reads "True" \
    --mapq "20" \
    --step "2" \
    --is_chm13 "False" -resume
