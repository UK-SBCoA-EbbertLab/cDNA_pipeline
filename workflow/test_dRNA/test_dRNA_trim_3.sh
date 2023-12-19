#!/bin/bash

nextflow ../main.nf --step 3 \
    --bambu_rds "./results/GRCh38_dRNA/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --fai "./results/GRCh38_dRNA/fai/*.fai" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --is_discovery "True" \
    --track_reads "False" \
    --NDR "auto" \
    --multiqc_input "./results/GRCh38_dRNA/multiQC_input/**" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --out_dir "./GRCh38_dRNA/" \
    --is_chm13 "False"
