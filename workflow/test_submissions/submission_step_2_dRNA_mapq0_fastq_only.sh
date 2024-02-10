#!/bin/bash

nextflow ../main.nf --step 2 \
    --ont_reads_fq "./results/test_mapq0_dRNA/concatenated_fastq_and_sequencing_summary_files/*.fastq" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./test_mapq0_dRNA/" \
    --is_chm13 "False" \
    --track_reads "True" \
    --mapq "0" \
    --quality_score "9" \
    --is_dRNA "True" \
    --trim_dRNA "False" 
