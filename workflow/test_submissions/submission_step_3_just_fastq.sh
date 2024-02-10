#!/bin/bash

nextflow ../main.nf --step 3 \
    --bambu_rds "./results/test_mapq10_just_fastqs/bambu_prep/*.rds" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --fai "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" \
    --annotation "../../references/Homo_sapiens.GRCh38.107.gtf" \
    --is_discovery "True" \
    --track_reads "True" \
    --NDR "auto" \
    --multiqc_input "./results/test_mapq10_just_fastqs/multiQC_input/**" \
    --multiqc_config "../../references/multiqc_config.yaml" \
    --out_dir "./test_mapq10_just_fastqs/" \
    --intermediate_qc "./results/test_mapq10_just_fastqs/intermediate_qc_reports/" \
    --is_chm13 "False" 


