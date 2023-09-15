process MAKE_CONTAMINATION_REPORT {

    publishDir "results/${params.out_dir}/multiQC_input/contamination/", mode: "copy", pattern: "*"

    label 'medium'

    input:
        val(id)
        val(num_reads)
        val(num_unmapped_reads)
        val(num_contaminant_reads)

    output:
        path("*"), emit: outty

    script:
        """

        echo "${id}" >> "${id}_test.txt"
        echo "${num_reads}" >> "${id}_test.txt"
        echo "${num_unmapped_reads}" >> "${id}_test.txt"
        echo "${num_contaminant_reads}" >> "${id}_test.txt"

        
        """

}

