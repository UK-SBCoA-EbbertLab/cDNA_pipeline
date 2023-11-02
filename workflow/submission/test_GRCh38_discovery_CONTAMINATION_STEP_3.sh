#!/usr/bin/env bash


nextflow ../main.nf --bambu_rds "./results/R10_V14_cDNA_Official_202310/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --out_dir "./R10_V14_cDNA_Official_202310/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --NDR "auto" \
    --multiqc_input "./results/R10_V14_cDNA_Official_202310/multiQC_input/**" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --fai "./results/R10_V14_cDNA_Official_202310/fai/*.fai" \
    --step 3 \
    --is_chm13 "False" -resume
