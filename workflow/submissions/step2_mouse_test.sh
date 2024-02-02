#!/bin/bash

nextflow ../main.nf --step 2 \
  --ont_reads_fq "../../../data/ont_data/big_experiment_data_transfer/Banner_Batch_02_20240113_60hrs/Batch_2_sample_1108_PAS87381_a667ecd3_cc3c2861-60hrs.fastq" \
  --ont_reads_txt "../../../data/ont_data/big_experiment_data_transfer/Banner_Batch_02_20240113_60hrs/Batch_2_sample_1108_PAS87381_a667ecd3_cc3c2861-60hrs.txt" \
  --ref "/project/mteb223_uksr/sequencing_resources/references/Ensembl/mouse_m39_soft_mask/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa" \
  --annotation "/project/mteb223_uksr/sequencing_resources/annotations/Ensembl/mouse_m39_release_110/Mus_musculus.GRCm39.110.gtf" \
  --out_dir "./Mouse_test/" \
  --cdna_kit "PCS114" \
  --is_chm13 "False" \
  --track_reads "False" \
  --mapq "10" \
  --contamination_ref "../../references/master_contaminant_reference.fasta" -resume
