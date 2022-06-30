process PYCOQC {

    publishDir "results/${params.out_dir}/pycoQC/"

    label 'large'

    input:
        val(id)
        path(seq_summary)
        path(total_bam)
        path(total_bai)

    output:
        path "pycoqc.json", emit: multiqc
        path "*", emit: all_output
    
    script:
        """
        pycoQC -f $seq_summary \
            -a $total_bam \
            -o "./${id}_pycoqc.html" \
            -j "./${id}_pycoqc.json" 
        """ 
}

