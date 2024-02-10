process MULTIQC_GRCh38 {

    publishDir "results/${params.out_dir}/multiQC_output", mode: "copy", overwrite: true

    label 'tiny'

    input:
        path(multiqc_input)
        path(multiqc_config)
    
    output: 
       path "*"

    script:
        """    
        multiqc -c $multiqc_config -n multiQC_report.html .
        """
}
