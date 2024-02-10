#!/bin/bash

nextflow ../main.nf --step 2 \
    --bam "./results/test_mapq10/mapping_cDNA/*.bam" \
    --bai "./results/test_mapq10/mapping_cDNA/*.bai" \
    --prefix "mapq20" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./test_mapq20_filtered_bams/" \
    --cdna_kit "PCS111" \
    --is_chm13 "False" \
    --track_reads "True" \
    --mapq "20" \
    --quality_score "9" \
    --is_dRNA "False" \
    --trim_dRNA "False" 
