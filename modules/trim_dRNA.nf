process TRIM_dRNA {

    publishDir "results/${params.out_dir}/trimmed_dRNA_files/", mode: 'symlink', overwrite: true, pattern: "*.trimmed.fastq"
    
    label "medium"

    input:
        tuple val(id), file(fastq)
        val(txt)

    output:
        tuple val("$id"), path("${id}.trimmed.fastq"), emit: fastq
        val("$txt"), emit: txt

    script:
    """
        
        ## Convert U to T since Dorado trimmer does not support U in fastq file sequences
        convert_U_to_T.py $fastq "${id}_fixed.fastq"

        ## Trim adapters and primers
        dorado trim --emit-fastq "${id}_fixed.fastq" > "${id}.trimmed.fastq" 

        ## Delete intermediate data
        rm "${id}_fixed.fastq"
    
    """
}
