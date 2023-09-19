process MAKE_CONTAMINATION_REPORT_1 {

    label 'small'

    input:
        val(id)
        val(num_reads)
        val(num_unmapped_reads)
        val(num_contaminant_reads)

    output:
        path("*.tsv")

    shell:
        '''
        
        mapped_to_contaminant="!{num_contaminant_reads}"
        unmapped_reads="!{num_unmapped_reads}"
        reads="!{num_reads}"
        ID="!{id}"
        
        unmapped=$(awk -v var1=$mapped_to_contaminant -v var2=$unmapped_reads 'BEGIN { print  ( var2 - var1 ) }') 
        mapped_to_target=$(awk -v var1=$unmapped_reads -v var2=$reads 'BEGIN { print  ( var2 - var1 ) }')

        echo "${ID}\t${mapped_to_target}\t${unmapped}\t${mapped_to_contaminant}" > "${ID}.tsv"
        '''

}

process MAKE_CONTAMINATION_REPORT_2 {

    publishDir "results/${params.out_dir}/multiQC_input/contamination/", mode: "copy", pattern: "*"

    label 'small'

    input:
        path(report)

    output:
        path("*"), emit: outty

    script:
        """
        
        echo "Sample ID\tMapped to Target\tUnmapped\tMapped to Contaminant" >> "Percent_Contaminant_Reads_mqc.tsv"
        cat $report >> "Percent_Contaminant_Reads_mqc.tsv"

        """

}

