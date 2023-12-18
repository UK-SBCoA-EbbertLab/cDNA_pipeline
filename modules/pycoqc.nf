process PYCOQC {

    publishDir "results/${params.out_dir}/multiQC_input/pycoqc/", mode: 'copy', overwrite: true, pattern: "*pycoqc*"

    label 'huge'

    input:
        val(id)
        path(fastq)
        path(seq_summary)
        path(total_bam)
        path(total_bai)
        val(quality_score)

    output:
        path "*pycoqc*", emit: multiQC

    script:
        """
        fix_sequencing_summary_pychopper.py $fastq $seq_summary "${id}_sequencing_summary_pyco.txt"


        pycoQC -f "${id}_sequencing_summary_pyco.txt" \
            -v \
            -a $total_bam \
            --min_pass_qual $quality_score \
            -o "./${id}_pycoqc.html" \
            -j "./${id}_pycoqc.json"
        """
}

process PYCOQC_dRNA {


    publishDir "results/${params.out_dir}/multiQC_input/pycoqc/", mode: 'copy', overwrite: true, pattern: "*pycoqc*"

    label 'huge'

    input:
        val(id)
        path(fastq)
        path(seq_summary)
        path(total_bam)
        path(total_bai)
        val(all_fastq_reads)
        val(mapq)
        path(stats)

    output:
        path "*pycoqc*", emit: multiQC
        tuple val("${id}"), env(num_reads_trimmed), val("${mapq}"), path("${id}_pycoqc.json"), path("${stats}"), emit: num_reads_report

    script:
        """
        num_reads_trimmed=\$(fix_sequencing_summary_porechop.py $fastq $seq_summary "${id}_sequencing_summary_pyco.txt"


        pycoQC -f "${id}_sequencing_summary_pyco.txt" \
            -v \
            -a $total_bam \
            --min_pass_qual 0 \
            -o "./${id}_pycoqc.html" \
            -j "./${id}_pycoqc.json"
        """
}

