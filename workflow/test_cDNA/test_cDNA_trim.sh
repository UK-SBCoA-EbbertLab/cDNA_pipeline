#!/bin/bash

nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/test_data/*.fastq" \
          --ont_reads_txt "/scratch/bag222/data/ont_data/test_data/*.txt" \
          --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
          --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
          --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
          --out_dir "./GRCh38_cDNA/" \
          --bambu_track_reads "True" \
          --is_dRNA "False" \
          --qscore_thresh "9" \
          --mapq "10" \
          --step "2" \
          --is_chm13 "False" -resume
