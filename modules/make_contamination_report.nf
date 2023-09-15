process MAKE_CONTAMINATION_REPORT {

    publishDir "results/${params.out_dir}/contamination_report/samtools_stat_files/", mode: "copy", pattern: "*.*stat"

    label 'medium'

    input:
        val(id)
        val(num_reads)
        val(num_unmapped_reads)
        val(num_contaminant_reads)

    output:
        path("${id}_test.txt"), emit: stat

    script:
        """

        echo "${id}" >> "${id}_test.txt"
        echo "${num_reads}" >> "${id}_test.txt"
        echo "${num_unmapped_reads}" >> "${id}_test.txt"
        echo "${num_contaminant_reads}" >> "${id}_test.txt"

        
        """

}

