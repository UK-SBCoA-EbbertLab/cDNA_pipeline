#!/usr/bin/env bash

nextflow run ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/technical_replicates_data/*.fastq" \
    --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2022_ont_data/technical_replicates_data/*.txt" \
    --ref "../../references/chm13v2.0_ERCC.fa" \
    --annotation "../../references/CHM13.v2.0.gff3" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed/" \
    --ercc "../../references/ERCC92.gtf" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --is_discovery "False" \
    --out_dir "./CHM13_repeats_cdna_comparison/" \
    --is_chm13 "True"  -resume -bg
