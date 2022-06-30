process RSEQC {

    publishDir "results/${params.out_dir}/RseQC/"
 
    label "medium_small"

    input:
        val(id)
        path(bam)
        path(bai)

    output:
        path "*"

    script:
        """
        bam_stat.py -i $bam > "${id}_bam_stat"
        read_GC.py -i $bam -o "${id}_read_gc"
        read_NVC.py -i $bam -o "${id}_read_nvc"
        read_quality.py -i $bam -o "${id}_read_quality"
        read_duplication.py -i $bam -o "${id}_read_duplication"
        """
}
