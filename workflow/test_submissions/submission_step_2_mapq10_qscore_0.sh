#!/bin/bash

nextflow ../main.nf --step 2 \
    --ont_reads_fq "../../../data/ont_data/test_data/*.fastq" \
    --ont_reads_txt "../../../data/ont_data/test_data/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./test_mapq10_qscore_0/" \
    --cdna_kit "PCS111" \
    --is_chm13 "False" \
    --track_reads "True" \
    --mapq "10" \
    --qscore_thresh "0" \
    --is_dRNA "False" \
    --trim_dRNA "False" 
