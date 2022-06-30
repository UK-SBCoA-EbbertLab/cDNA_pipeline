process SEQ_SUMMARY {

    publishDir "results/${params.out_dir}/pychopper/"

    label "huge"

    input:
        val(id)
        path(fastq)
        path(sequencing_summary)

    output:
        val "$id", emit: id
        path "$fastq", emit: fastq
        path "${id}_sequencing_summary_pychop.txt", emit: txt

    script:
        """
        fix_sequencing_summary_pychopper.py $fastq $sequencing_summary "${id}_sequencing_summary_pychop.txt"
        """
}
