process dRNA_TRIM {

    publishDir "results/${params.out_dir}/concatenated_fastq_and_sequencing_summary_files/", mode: 'copy', overwrite: true, pattern: "*.fastq"
    
    label "medium"

    input:
        tuple val(id), file(fastq)

    output:
        tuple val("$id"), path("${id}.trimmed.fastq")

    script:
    """
        cat $fastq

        dorado trim $fastq
        
    """
}
