process MULTIQC_GRCh38 {

    publishDir "results/${params.out_dir}/multiqc/"

    label 'tiny'

    input:
        path(QC_1)
        path(QC_2)
        path(QC_3)
        path(QC_4)
    
    output: 
        path "*"

    script:
        """    
        multiqc -c multiqc_config.yaml -n multiQC_report.html .
        """
}

process MULTIQC_CHM13 {

    publishDir "results/${params.out_dir}/multiqc/"

    label 'tiny'

    input:
        path(QC_1)
        path(QC_2)
        path(QC_3)
    
    output: 
        path "*"

    script:
        """    
        multiqc -c multiqc_config.yaml -n multiQC_report.html .
        """
}
