#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/R10_V14_cDNA_Test_202310/R10_official_test/*.fastq" \
    --ont_reads_txt "/scratch/bag222/data/ont_data/R10_V14_cDNA_Test_202310/R10_official_test/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --cdna_kit "PCS114" \
    --out_dir "./R10_V14_cDNA_Official_202310/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --mapq "10" \
    --step "2" \
    --contamination_ref "../../references/master_contaminant_reference.fasta" \
    --is_chm13 "False" -resume
