#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/uky_aged_data/*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/uky_aged_data/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.106_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed/" \
    --out_dir "./GRCh38-106_aged_uky_ont/" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --is_discovery "True" \
    --is_chm13 "False"  -resume
