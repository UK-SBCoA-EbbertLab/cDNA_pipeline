process FAST5_to_POD5 {

    publishDir "results/${params.out_dir}/fast5_to_pod5/", mode: "copy", overwrite: true
    
    label 'large'

    input:
        tuple val(id), path(fast5)

    output:
        tuple val("${id}"), path("*.pod5")


    script:
        """

        pod5 convert fast5 *.fast5 --output . --one-to-one . --threads 50
        
        """

}


process BASECALL {

    label 'medium_large'

    input:
        tuple val(id), path(pod5_dir)
        val speed
        val modifications

    output:
        tuple val("${id}"), path('*.fastq'), emit: fastq
        val '${id}.txt', emit: txt

   script:
        """
        
        dorado basecaller "${speed}" . -x cpu > "${id}.bam"
        
        dorado summary "${id}.bam" > "{id}.txt"
        
        samtools fastq -T "*" "${id}.bam" > "${id}.fastq"
        
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
