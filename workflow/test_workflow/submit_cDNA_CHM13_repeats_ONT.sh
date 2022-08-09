#!/usr/bin/env bash

nextflow ../main.nf --ont_reads_fq "/home/bag222/mnt/repeated_flow_cells_data/*.fastq" \
    --ont_reads_txt "/home/bag222/mnt/repeated_flow_cells_data/*.txt" \
    --ref "../../references/chm13v2.0_ERCC.fa" \
    --annotation "../../references/CHM13.v2.0.gff3" \
    --ercc "../../references/ERCC92.gtf" \
    --out_dir "./CHM13_repeats_ONT/" \
    --is_chm13 "True"  -resume
