process PYCOQC {

    publishDir "results/${params.out_dir}/pycoQC/", mode: "copy", overwrite: true

    label 'huge'

    input:
        val(id)
        path(seq_summary)
        path(total_bam)
        path(total_bai)

    output:
        path "*", emit: all_output

    script:
        """
        fix_sequencing_summary.py $seq_summary "${id}_sequencing_summary_pycoqc.txt"

        pycoQC -f "${id}_sequencing_summary_pycoqc.txt" \
            -v \
            -a $total_bam \
            -o "./${id}_pycoqc.html" \
            -j "./${id}_pycoqc.json"
        """
}

