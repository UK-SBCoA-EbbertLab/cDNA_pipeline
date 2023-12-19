#!/bin/bash

nextflow ../main.nf --path "/scratch/bag222/data/ont_data/test_dRNA_data/basecalling/" \
          --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
          --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
          --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
          --out_dir "./GRCh38_dRNA/" \
          --bambu_track_reads "True" \
          --is_dRNA "True" \
          --qscore_thresh "9" \
          --mapq "10" \
          --step "2" \
          --trim_dRNA "True" \
          --is_chm13 "False" -resume
