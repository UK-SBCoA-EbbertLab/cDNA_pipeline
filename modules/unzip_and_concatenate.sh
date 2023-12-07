process UNZIP_AND_CONCATENATE {

    label "medium"

    input:
        tuple val(id), file(reads)

    output:
        val "$id", emit: id
        path "${id}.fastq", emit: fastq

    """
    
    find . -type f -name "*.fastq.gz" -exec gunzip {} +

    find . -type f -name "*.fastq" -exec cat {} + > "${id}.fastq"
    
    """
}
