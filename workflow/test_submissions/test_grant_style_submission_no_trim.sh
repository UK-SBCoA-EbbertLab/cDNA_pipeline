#!/bin/bash


nextflow ../main.nf --ont_reads_fq "/pscratch/mteb223_uksr/grant_dRNA_data_for_test/data/*.fastq" \
                --ont_reads_txt "/pscratch/mteb223_uksr/grant_dRNA_data_for_test/data/*.txt" \
                --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
                --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
                --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
                --out_dir "./HEK_test_no_trim/" \
                --bambu_track_reads "True" \
                --is_dRNA "True" \
                --trim_dRNA "False" \
                --mapq "0" \
                --step "2" \
                --is_chm13 "False"

