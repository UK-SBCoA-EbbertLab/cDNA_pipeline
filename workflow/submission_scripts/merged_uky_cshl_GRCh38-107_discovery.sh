#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/merged_uky_and_cshl_aged_data/*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/merged_uky_and_cshl_aged_data/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./merged_uky_cshl_GRCh38-107_discovery/" \
    --is_discovery "True" \
    --is_chm13 "False" -resume -bg
