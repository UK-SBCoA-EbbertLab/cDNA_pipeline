#!/usr/bin/env bash


nextflow ../main.nf --bambu_rds "./results/GTEx_GRCh38_version_107_mapq_10_track_reads_multithreaded/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
    --out_dir "./GTEx_GRCh38_version_107_mapq_10_track_reads_multithreaded/" \
    --is_discovery "True" \
    --bambu_track_reads "True" \
    --NDR "auto" \
    --multiqc_input "./results/GTEx_GRCh38_version_107_mapq_10_track_reads_multithreaded/multiQC_input/**" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --fai "./results/GTEx_GRCh38_version_107_mapq_10_track_reads_multithreaded/fai/*.fai" \
    --step 3 \
    --is_chm13 "False" -resume
