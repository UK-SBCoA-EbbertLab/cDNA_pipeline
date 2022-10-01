process PYCOQC {

    publishDir "results/${params.out_dir}/QC/pycoqc/", mode: 'copy', overwrite: true, pattern: "*pycoqc*"

    label 'huge'

    input:
        val(id)
        path(seq_summary)
        path(total_bam)
        path(total_bai)

    output:
        path "*pycoqc*", emit: multiQC

    script:
        """
        fix_sequencing_summary.py $seq_summary "${id}_sequencing_summary_pyco.txt"

        pycoQC -f "${id}_sequencing_summary_pyco.txt" \
            -v \
            -a $total_bam \
            -o "./${id}_pycoqc.html" \
            -j "./${id}_pycoqc.json"
        """
}

