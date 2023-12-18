#!/bin/bash

nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/test_dRNA_data/*.fastq" \
          --ont_reads_txt "/scratch/bag222/data/ont_data/test_dRNA_data/*.txt" \
          --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
          --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
          --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
          --out_dir "./GRCh38_dRNA/" \
          --bambu_track_reads "True" \
          --is_dRNA "True" \
          --mapq "10" \
          --step "2" \
          --trim_dRNA "True" \
          --is_chm13 "False" -resume
