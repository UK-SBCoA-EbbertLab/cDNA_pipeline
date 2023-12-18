process TRIM_dRNA {

    publishDir "results/${params.out_dir}/trimmed_dRNA_files/", mode: 'symlink', overwrite: true, pattern: "*"
    
    label "medium"

    input:
        tuple val(id), file(fastq)
        val(txt)

    output:
        tuple val("$id"), path("${id}.trimmed.fastq"), emit: fastq
        val("$txt"), emit: txt
        path("*adapter_data*.txt"), emit: outty
    
    script:
    """
        ## Trim adapters and primers
        porechop_abi -abi -i "${fastq}" -o "${id}.trimmed.fastq" > "${id}_adapter_data.txt"
       
    """
}
