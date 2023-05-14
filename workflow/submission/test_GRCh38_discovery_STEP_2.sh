#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/GTEx_data_maddy/sequence_data/PCR_cDNA/*.fastq" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --cdna_kit "PCS109" \
    --out_dir "./GTEx_GRCh38_version_107_mapq_10_track_reads_multithreaded/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --mapq "10" \
    --step "2" \
    --is_chm13 "False" -resume
