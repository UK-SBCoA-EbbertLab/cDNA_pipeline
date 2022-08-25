process SEQ_SUMMARY {

    label "huge"

    input:
        val(id)
        path(seq_summary)
        path(total_bam)
        path(total_bai)

    output:
        val "$id", emit: id
        path "${id}_sequencing_summary_fixed.txt", emit: txt
        path "$total_bam", emit: bam
        path "$total_bai", emit: bai


    script:
        """
        fix_sequencing_summary.py $seq_summary "${id}_sequencing_summary_fixed.txt"
        """
}
