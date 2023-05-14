#!/usr/bin/env bash


nextflow ../main.nf --bambu_rds "./results/test_discovery/bambu_prep/*.rds" \
    --ref "../../references/chm13v2.0_ERCC.fa" \
    --annotation "./results/test_discovery/CHM13_gtf/CHM13_v2.0_ERCC.gtf" \
    --out_dir "./test_discovery/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --NDR "auto" \
    --fai "./results/test_discovery/fai/chm13v2.0_ERCC.fa.fai" \
    --multiqc_input "./results/test_discovery/multiQC_input/**" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --fai "./results/test_discovery/fai/chm13v2.0_ERCC.fa.fai" \
    --step 3 \
    --is_chm13 "True" -resume
