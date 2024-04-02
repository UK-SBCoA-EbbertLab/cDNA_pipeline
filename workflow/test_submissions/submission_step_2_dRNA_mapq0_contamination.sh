#!/bin/bash

nextflow ../main.nf --step 2 \
    --path "../../../data/ont_data/test_dRNA_data/basecalling/" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --out_dir "./test_mapq0_dRNA_contaminant/" \
    --is_chm13 "False" \
    --track_reads "True" \
    --contamination_ref "../../references/master_contaminant_reference.fasta" \
    --mapq "0" \
    --quality_score "9" \
    --is_dRNA "True" \
    --trim_dRNA "False" -resume
