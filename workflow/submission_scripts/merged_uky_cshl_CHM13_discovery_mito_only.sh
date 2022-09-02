#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/merged_uky_and_cshl_aged_data/*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/merged_uky_and_cshl_aged_data/*.txt" \
    --ref "../../references/chm13v2.0_ERCC_mito.fa" \
    --annotation "../../references/CHM13.v2.0._mito.gff3" \
    --ercc "../../references/ERCC92.gtf" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./merged_uky_cshl_CHM13_discovery_mito_only/" \
    --is_discovery "True" \
    --is_chm13 "True" -resume -bg
