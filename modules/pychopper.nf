process PYCHOPPER {

    publishDir "results/${params.out_dir}/multiQC_input/pychopper/", mode: 'copy', overwrite: true, pattern: "*pychopper.stats"

    label "large"

    input:
        tuple val(id), path(fastq)
        val(txt)
        val(cdna_kit)

    output:
        val "$id", emit: id
        path "${id}_pychop.fq", emit: fastq
        val "$txt", emit: txt
        path "$fastq", emit: original_fastq
        path "*pychopper.stats", emit: multiQC

    script:
    
    """

    ## Pychopper does not have PCS114 primers yes, need to create them ##
    if [[ "${cdna_kit}" == "PCS114" ]];
    then 
        ## Create primer config file ##
        echo "+:MySSP,-MyVNP|-:MyVNP,-MySSP" > primer_config.txt
    
        ## Create custom primers (using PCS111 here for test, neew to add PCS114 once we have access ##
        echo ">MyVNP" > custom_pimers.fas
        echo "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT" >> custom_pimers.fas
        echo ">MySSP" >> custom_pimers.fas
        echo "TTTCTGTTGGTGCTGATATTGCTTT" >> custom_pimers.fas

        ## Run pychopper with the custom primers and primer config ##
        pychopper -m edlib -b custom_pimers.fas -c primer_config.txt \
            -t 50 \
            -Q 9 \
            -k $cdna_kit \
            -r "${id}_pychopper_report.pdf" \
            -u "${id}_pychopper.unclassified.fq" \
            -w "${id}_pychopper.rescued.fq" \
            -S "${id}_pychopper.stats" \
            -A "${id}_pychopper.scores" \
            "${fastq}" "${id}_pychop.fq"

    ## All other kits just use default settings ##
    else

        pychopper -t 50 \
            -Q 9 \
            -k $cdna_kit \
            -r "${id}_pychopper_report.pdf" \
            -u "${id}_pychopper.unclassified.fq" \
            -w "${id}_pychopper.rescued.fq" \
            -S "${id}_pychopper.stats" \
            -A "${id}_pychopper.scores" \
            "${fastq}" "${id}_pychop.fq"
    """
}
