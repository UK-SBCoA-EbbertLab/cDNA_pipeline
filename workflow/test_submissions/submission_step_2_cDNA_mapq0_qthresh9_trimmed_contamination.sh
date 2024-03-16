#!/bin/bash

nextflow ../main.nf --step 2 \
    --ont_reads_fq "../../../data/ont_data/test_data/*.fastq" \
    --ont_reads_txt "../../../data/ont_data/test_data/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./test_mapq0_qthresh9_cDNA_trimmed_contamination/" \
    --is_chm13 "False" \
    --track_reads "True" \
    --contamination_ref "../../references/master_contaminant_reference.fasta" \
    --mapq "0" \
    --qscore_thresh "9" \
    --is_dRNA "False" 
