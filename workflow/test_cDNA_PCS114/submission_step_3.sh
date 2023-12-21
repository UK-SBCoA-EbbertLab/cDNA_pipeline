#!/bin/bash

nextflow ../main.nf --step 3 \
    --bambu_rds "./results/PCS114/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --fai "./results/PCS114/fai/*.fai" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --is_discovery "True" \
    --track_reads "True" \
    --NDR "auto" \
    --multiqc_input "./results/PCS114/multiQC_input/**" \
    --multiqc_config "../../../references/multiqc_config.yaml" \
    --out_dir "./PCS114/" \
    --intermediate_qc "./results/PCS114/intermediate_qc_reports/" \
    --is_chm13 "False"  


