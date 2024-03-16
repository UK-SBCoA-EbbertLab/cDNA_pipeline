#!/bin/bash

nextflow ../main.nf --step 2 \
    --path "../../../data/ont_data/test_dRNA_data/basecalling/" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./test_mapq10_qthresh0_dRNA_trimmed_contamination/" \
    --is_chm13 "False" \
    --track_reads "True" \
    --contamination_ref "../../references/master_contaminant_reference.fasta" \
    --mapq "10" \
    --trim_dRNA "False" \
    --qscore_thresh "0" \
    --is_dRNA "True" \
    --trim_dRNA "True" 
