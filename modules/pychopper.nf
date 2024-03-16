process PYCHOPPER {

    publishDir "results/${params.out_dir}/multiQC_input/pychopper/", mode: 'copy', overwrite: true, pattern: "*pychopper.stats"
    
    label "large"

    input:
        tuple val(id), path(fastq)
        val(txt)
        val(cdna_kit)
        val(quality_score)

    output:
        val "$id", emit: id
        path "${id}_pychop.fq", emit: fastq
        path "${id}.txt", emit: txt     
        path "$fastq", emit: original_fastq
        path "*pychopper.stats", emit: multiQC
        env(NUM_PASS_READS), emit: num_pass_reads

    script:
    
    """
    if [[ "${txt}" != "None" ]] &&  [[ "${txt}" != "${id}.txt" ]]; then
        cp "${txt}" "./${id}.txt"
    elif [[ "${txt}" == "None" ]]; then 
        touch "./${id}.txt"
    fi

    NUM_PASS_READS=\$(filter_by_mean_base_quality.py "${fastq}" "${quality_score}" "${id}_qscore_${quality_score}.fastq")
    
    ## Pychopper does not have PCS114 primers yes, need to create them ##
    if [[ "${cdna_kit}" == "PCS114" ]]; then
    
        ## Create primer config file ##
        echo "+:MySSP,-MyVNP|-:MyVNP,-MySSP" > primer_config.txt
    
        ## Create custom primers for PCS114 ##
        echo ">MyVNP" > custom_pimers.fas
        echo "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT" >> custom_pimers.fas
        echo ">MySSP" >> custom_pimers.fas
        echo "TTTCTGTTGGTGCTGATATTGCTTT" >> custom_pimers.fas

        ## Run pychopper with the custom primers and primer config ##
        pychopper -m edlib -b custom_pimers.fas -c primer_config.txt \
            -t 50 \
            -Q $quality_score \
            -r "${id}_pychopper_report.pdf" \
            -u "${id}_pychopper.unclassified.fq" \
            -w "${id}_pychopper.rescued.fq" \
            -S "${id}_pychopper.stats" \
            -A "${id}_pychopper.scores" \
            "${id}_qscore_${quality_score}.fastq" "${id}_pychop.fq"

    ## All other kits just use default settings ##
    else

        pychopper -m edlib \
            -t 50 \
            -Q $quality_score \
            -k $cdna_kit \
            -r "${id}_pychopper_report.pdf" \
            -u "${id}_pychopper.unclassified.fq" \
            -w "${id}_pychopper.rescued.fq" \
            -S "${id}_pychopper.stats" \
            -A "${id}_pychopper.scores" \
            "${id}_qscore_${quality_score}.fastq" "${id}_pychop.fq"
    
    fi

    rm "${id}_qscore_${quality_score}.fastq"

    """
}
