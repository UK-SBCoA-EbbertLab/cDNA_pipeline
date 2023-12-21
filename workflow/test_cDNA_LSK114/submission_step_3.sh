#!/bin/bash

nextflow ../main.nf --step 3 \
    --bambu_rds "./results/LSK114/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --fai "./results/LSK114/fai/*.fai" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --is_discovery "True" \
    --track_reads "True" \
    --NDR "auto" \
    --multiqc_input "./results/LSK114/multiQC_input/**" \
    --multiqc_config "../../../references/multiqc_config.yaml" \
    --out_dir "./LSK114/" \
    --intermediate_qc "./results/LSK114/intermediate_qc_reports/" \
    --is_chm13 "False"  


