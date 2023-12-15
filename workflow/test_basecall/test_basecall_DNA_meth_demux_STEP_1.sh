#!/bin/bash

nextflow ../main.nf --step 1 --basecall_path "/scratch/bag222/data/test_DNA_methylayion_demux/" \
    --basecall_speed "hac" \
    --basecall_mods "5mCG_5hmCG" \
    --basecall_config "False" \
    --basecall_trim "all" \
    --basecall_compute "cpu" \
    --basecall_demux "True" \
    --outdir "test_basecall_DNA" 
