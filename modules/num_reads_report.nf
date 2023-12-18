process MAKE_NUM_READS_REPORT {

    publishDir "results/${params.out_dir}/multiQC_input/num_reads_report/", pattern: "*", mode: "copy", overwrite: true

    label 'large'

    input:
        tuple val(id), val(num_all_fastq_reads),  val(mapq), path(json), path(flagstat)

    output:
        path("*"), emit: outty

    script:
        """

        echo "${id}, ${num_all_fastq_reads}, ${mapq}, ${json}, ${flagstat}"

        echo "${id}, ${num_all_fastq_reads}, ${mapq}, ${json}, ${flagstat}" > "${id}.tsv"

        """

}
