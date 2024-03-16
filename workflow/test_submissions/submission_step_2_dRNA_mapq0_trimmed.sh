#!/bin/bash

nextflow ../main.nf --step 2 \
    --path "../../../data/ont_data/test_dRNA_data/basecalling/" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./test_mapq0_dRNA_trimmed/" \
    --is_chm13 "False" \
    --track_reads "True" \
    --mapq "0" \
    --trim_dRNA "False" \
    --quality_score "9" \
    --is_dRNA "True" \
    --trim_dRNA "True"  
