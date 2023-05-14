#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/test_data/*.fastq" \
    --ont_reads_txt "/scratch/bag222/data/ont_data/test_data/*.txt" \
    --ref "../../references/chm13v2.0_ERCC.fa" \
    --annotation "../../references/CHM13.v2.0.gff3" \
    --ercc "../../references/ERCC92.gtf" \
    --cdna_kit "PCS111" \
    --out_dir "./test_discovery/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --mapq "10" \
    --step "2" \
    --is_chm13 "True" -resume
