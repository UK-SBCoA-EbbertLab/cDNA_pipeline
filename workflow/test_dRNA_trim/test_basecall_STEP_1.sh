#!/bin/bash

nextflow ../main.nf --step 1 --basecall_path "/scratch/bag222/data/ont_data/test_dRNA_data/basecalling_2/" \
    --basecall_speed "fast" \
    --basecall_mods "" -resume
