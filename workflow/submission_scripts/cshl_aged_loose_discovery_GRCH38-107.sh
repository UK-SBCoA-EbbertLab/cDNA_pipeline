#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/cshl_aged_data/*.fastq" \
    --ont_reads_txt "/scratch/bag222/data/ont_data/cshl_aged_data/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./cshl_aged_loose/" \
    --NDR "1.00" \
    --is_discovery "True" \
    --is_chm13 "False" -resume