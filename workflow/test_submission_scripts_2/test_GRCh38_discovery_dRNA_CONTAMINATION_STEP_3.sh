#!/usr/bin/env bash


nextflow ../main.nf --bambu_rds "./results/test_discovery_dRNA/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --out_dir "./test_discovery_dRNA/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --NDR "auto" \
    --multiqc_input "./results/test_discovery_dRNA/multiQC_input/**" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --fai "./results/test_discovery_dRNA/fai/*.fai" \
    --step 3 \
    --is_chm13 "False" -resume
