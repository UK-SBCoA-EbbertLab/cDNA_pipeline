process RSEQC {

    publishDir "results/${params.out_dir}/multiQC_input/RseQC/", mode: "copy", overwrite: true, pattern: "${id}*"

 
    label "large"

    input:
        val(id)
        path(bam)
        path(bai)
        path(housekeeping)
        path(annotation)
        val(mapq)

    output:
        path "${id}*", emit: multiQC

    script:
        """

        geneBody_coverage.py -i $bam -r $housekeeping -o "${id}_RSEQC_geneBody_coverage"

        convert_gtf_to_bed12.py "${annotation}"
        
        junction_annotation.py -i "${bam}" -o "${id}_RSEQC_junction_annotation" -r *.bed
       
        junction_saturation.py -i "${bam}" -o "${id}_RSEQC_junction_saturation" -r *.bed

        tin.py -i "${bam}" -r *.bed > "${id}_RSEQC_TIN.txt" 

        bam_stat.py -i "${bam}" -q "${mapq}" > "${id}_RSEQC_bam_stat.txt"

        read_distribution.py -i "${bam}" -r *.bed > "${id}_RSEQC_read_distribution.txt"

        read_GC.py -i "${bam}" -q "${mapq}" -o "${id}_RSEQC_read_GC"       
        

        """
}
