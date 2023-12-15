process CONVERT_U_TO_T {

    publishDir "results/${params.out_dir}/trimmed_dRNA_files/", mode: 'symlink', overwrite: true, pattern: "*"
    
    label "medium"

    input:
        tuple val(id), file(fastq)
        val(txt)

    output:
        tuple val("$id"), path("${id}_U_to_T.fastq"), emit: fastq
        val("$txt"), emit: txt

    script:
    """
        
        ## convert U to T
        convert_U_to_T.py $fastq "${id}_U_to_T.fastq"

    """
}
