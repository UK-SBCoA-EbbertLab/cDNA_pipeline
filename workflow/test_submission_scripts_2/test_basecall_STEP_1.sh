#!/bin/bash

nextflow ../main.nf --step 1 --basecall_path "/scratch/bag222/data/ont_data/test_dRNA_data/basecalling/" \
    --basecall_speed "hac" \
    --basecall_mods "" -resume
