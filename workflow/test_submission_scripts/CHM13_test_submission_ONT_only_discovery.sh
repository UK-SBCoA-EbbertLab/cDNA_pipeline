#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/test_data/*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/test_data/*.txt" \
    --ref "../../references/chm13v2.0_ERCC.fa" \
    --annotation "../../references/CHM13.v2.0.gff3" \
    --ercc "../../references/ERCC92.gtf" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --cdna_kit "PCS111" \
    --out_dir "./CHM13_comparison_test_ONT_only_discovery/" \
    --is_discovery "True" \
    --is_chm13 "True" -resume -bg
