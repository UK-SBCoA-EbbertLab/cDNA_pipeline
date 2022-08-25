process PYCOQC {

    publishDir "results/${params.out_dir}/pycoQC/", mode: "copy", overwrite: true

    label 'huge'

    input:
        path(seq_summary)
        path(total_bam)
        path(total_bai)

    output:
        path "*", emit: all_output
    
    script:
        """
        pycoQC -f $seq_summary \
            -v \
            -a $total_bam \
            -o "./pycoqc.html" \
            -j "./pycoqc.json" 
        """ 
}

