process FAST5_to_POD5 {

    publishDir "results/${params.out_dir}/fast5_to_pod5/", mode: "copy", overwrite: true
    
    label 'normal'

    input:
        tuple val(id), path(fast5_dir)

    output:
        tuple val("${id}"), path("${id}_pod5s/")


    script:
        """
        
        pod5 convert fast5 "./${fast5_dir}/*.fast5" --output "{$id}_pod5s/" --one-to-one "./${fast5_dir}/"
        
        """

}



process BASECALL {

    label 'gpu'

    input:
        tuple val(id), path(pod5_dir)
        val speed
        val modifications

    output:
        tuple val(id), path('pass/*.fastq'), emit: fastq
        val '*.txt', emit: txt

   script:
        """

        dorado basecaller "${speed},{modifications}" $pod5_dir --emit-fastq --no-trim
        
        """
}


process GATHER_BASECALL {

    publishDir "results/${params.out_dir}/basecall_output/", mode: "copy", overwrite: true

    label 'local'

    input:
        tuple val(id), path(fastq)
        val txt

    output:
        tuple val("$id"), path( '${id}.fastq'), emit: fastq
        val '*.txt', emit: txt

    script:
        """
        
        find . -type f -maxdepth 1 -name "*.fastq" ! -name "${id}.fastq" -exec cat {} \\; >> "${id}.fastq"
        
        """
}
