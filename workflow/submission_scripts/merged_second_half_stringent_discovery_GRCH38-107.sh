#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data_split/*secondHalf.fastq" \
    --ont_reads_txt "/scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data_split/*secondHalf.txt" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./merged_second_half_stringent/" \
    --NDR "0.1" \
    --is_discovery "True" \
    --is_chm13 "False" -resume
