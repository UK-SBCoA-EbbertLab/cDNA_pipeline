#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.106.gtf" \
    --ercc "None" \
    --out_dir "./GRCh38_test/" \
    --cdna_kit "PCS111" \
    --is_chm13 "False"  -resume
