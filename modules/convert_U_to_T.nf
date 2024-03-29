process CONVERT_U_TO_T {

    label "medium"

    input:
        tuple val(id), file(fastq)
        val(txt) 
        val(qscore)

    output:
        tuple val("$id"), path("${id}_U_to_T_qscore_${qscore}.fastq"), emit: fastq
        path "${id}.txt", emit: txt
        env(NUM_PASS_READS), emit: num_pass_reads

    script:
    """
       
        if [[ "${txt}" != "None" ]] &&  [[ "${txt}" != "${id}.txt" ]]; then
            cp "${txt}" "./${id}.txt"
        elif [[ "${txt}" == "None" ]]; then
            touch "./${id}.txt"
         fi
 
        ## convert U to T
        convert_U_to_T.py $fastq "${id}_U_to_T.fastq"

        ## Filter by mean base quality threshold
        NUM_PASS_READS=\$(filter_by_mean_base_quality.py "${id}_U_to_T.fastq" "${qscore}" "${id}_U_to_T_qscore_${qscore}.fastq")

        ## Delete intermediate files
        rm "${id}_U_to_T.fastq"

    """
}
