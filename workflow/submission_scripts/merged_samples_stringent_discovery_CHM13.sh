#!/usr/bin/env bash


nextflow ../../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/*.fastq" \
    --ont_reads_txt "/scratch/bag222/data/ont_data/merged_uky_and_cshl_aged_data/*.txt" \
    --ref "../../../references/chm13v2.0_ERCC.fa" \
    --annotation "../../../references/CHM13.v2.0.gff3" \
    --ercc "../../../references/ERCC92.gtf" \
    --multiqc_config "../../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./CHM13_merged_aged_stringent/" \
    --NDR "0.1" \
    --is_discovery "True" \
    --is_chm13 "True" -resume
