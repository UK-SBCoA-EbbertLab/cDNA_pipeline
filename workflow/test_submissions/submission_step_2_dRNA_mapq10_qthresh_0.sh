#!/bin/bash

nextflow ../main.nf --step 2 \
    --path "../../../data/ont_data/test_dRNA_data/basecalling/" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./test_mapq10_dRNA_qthresh_0/" \
    --is_chm13 "False" \
    --track_reads "True" \
    --mapq "10" \
    --qscore_thresh "0" \
    --is_dRNA "True" \
    --trim_dRNA "False" 
