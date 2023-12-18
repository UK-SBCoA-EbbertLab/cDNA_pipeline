#!/usr/bin/env bash


nextflow ../main.nf --bambu_rds "./results/Banner_Batch_1_120823/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --out_dir "./Banner_Batch_1_120823/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --NDR "auto" \
    --multiqc_input "./results/Banner_Batch_1_120823/multiQC_input/**" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --fai "./results/Banner_Batch_1_120823/fai/*.fai" \
    --step 3 \
    --is_chm13 "False" -resume
