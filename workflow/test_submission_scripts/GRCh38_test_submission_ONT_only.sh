#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/ont_test_data/*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/ont_test_data/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.101_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./GRCh38_cDNA_test_ONT_only/" \
    --is_discovery "False" \
    --is_chm13 "False" -resume -bg
