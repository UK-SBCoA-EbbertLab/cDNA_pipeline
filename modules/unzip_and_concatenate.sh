process UNZIP_AND_CONCATENATE {

    publishDir "results/${params.out_dir}/concatenated_fastq_and_sequencing_summary_files/", mode: 'copy', overwrite: true, pattern: "*.fastq"
    
    label "medium"

    input:
        tuple val(id), file(reads)

    output:
        tuple val("$id"), path("${id}.fastq")

    script:
    """

        echo "HELLO"

        find -L . -maxdepth 1 -name "*.fastq.gz" -exec gunzip --keep --force {} \\;
        
        echo "MID"

        find . -type f -maxdepth 1 -name "*.fastq" -exec cat {} \\; >> "${id}.fastq"
    
        echo "FINAL"

    """
}

process FIX_SEQUENCING_SUMMARY_NAME {

    publishDir "results/${params.out_dir}/concatenated_fastq_and_sequencing_summary_files/", mode: 'copy', overwrite: true, pattern: "*.txt"

    label "small"

    input: 
        tuple val(id), file(txt)

    output:
        path "${id}.txt"

    script:
    """
    
        mv $txt "${id}.txt"
    
    """

}
