#!/bin/bash


for dir in ./results/*
do
    dir=$(basename "${dir}")

    rclone sync --transfers 64 --checkers 64 --copy-links --progress \
        "./results/${dir}/" "gemini1-1:/mnt/gemini1-7/mteb223_uksr/bag222/processed_data/2022_transcript_discovery_pipeline_output/${dir}/"

done
