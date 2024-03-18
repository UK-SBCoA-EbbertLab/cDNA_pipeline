process RSEQC_GENE_BODY_COVERAGE {

    publishDir "results/${params.out_dir}/multiQC_input/RseQC/", mode: "copy", overwrite: true, pattern: "${id}*"

 
    label "large"

    input:
        val(id)
        path(bam)
        path(bai)
        path(housekeeping)

    output:
        path "${id}*", emit: multiQC

    script:
        """

        geneBody_coverage.py -i $bam -r $housekeeping -o "${id}_RSEQC_geneBody_coverage"

        """
}

process RSEQC_BAM_STAT {

    publishDir "results/${params.out_dir}/multiQC_input/RseQC/", mode: "copy", overwrite: true, pattern: "${id}*"

 
    label "large"

    input:
        val(id)
        path(bam)
        path(bai)
        val(mapq)

    output:
        path "${id}*", emit: multiQC

    script:
        """

        bam_stat.py -i "${bam}" -q "${mapq}" > "${id}_RSEQC_bam_stat.txt"

        """
}

process RSEQC_READ_GC {

    publishDir "results/${params.out_dir}/multiQC_input/RseQC/", mode: "copy", overwrite: true, pattern: "${id}*"

 
    label "large"

    input:
        val(id)
        path(bam)
        path(bai)
        val(mapq)

    output:
        path "${id}*", emit: multiQC

    script:
        """

        read_GC.py -i "${bam}" -q "${mapq}" -o "${id}_RSEQC_read_GC"       
        
        """
}


process CONVERT_GTF_TO_BED12 {

    publishDir "results/${params.out_dir}/bed_annotation/", mode: "copy", overwrite: true, pattern: "*.bed"

 
    label "large"

    input:
        path(annotation)

    output:
        path "*.bed", emit: bed12

    script:
        """
        
        convert_gtf_to_bed12.py "${annotation}"

        """
}

process RSEQC_JUNCTION_ANNOTATION {

    publishDir "results/${params.out_dir}/multiQC_input/RseQC/", mode: "copy", overwrite: true, pattern: "${id}*"

 
    label "large"

    input:
        val(id)
        path(bam)
        path(bai)
        path(bed)

    output:
        path "${id}*", emit: multiQC

    script:
        """

        junction_annotation.py -i "${bam}" -o "${id}_RSEQC_junction_annotation" -r "${bed}"

        """
}


process RSEQC_JUNCTION_SATURATION {

    publishDir "results/${params.out_dir}/multiQC_input/RseQC/", mode: "copy", overwrite: true, pattern: "${id}*"

 
    label "large"

    input:
        val(id)
        path(bam)
        path(bai)
        path(bed)

    output:
        path "${id}*", emit: multiQC

    script:
        """

        junction_saturation.py -i "${bam}" -o "${id}_RSEQC_junction_saturation" -r "${bed}"

        """
}

process RSEQC_TIN {

    publishDir "results/${params.out_dir}/multiQC_input/RseQC/", mode: "copy", overwrite: true, pattern: "${id}*"

 
    label "large"

    input:
        val(id)
        path(bam)
        path(bai)
        path(bed)

    output:
        path "${id}*", emit: multiQC

    script:
        """

        tin.py -i "${bam}" -r "${bed}" > "${id}_RSEQC_TIN.txt"

        """
}


process RSEQC_READ_DISTRIBUTION {

    publishDir "results/${params.out_dir}/multiQC_input/RseQC/", mode: "copy", overwrite: true, pattern: "${id}*"

 
    label "large"

    input:
        val(id)
        path(bam)
        path(bai)
        path(bed)

    output:
        path "${id}*", emit: multiQC

    script:
        """

        read_distribution.py -i "${bam}" -r "${bed}" > "${id}_RSEQC_read_distribution.txt"

        """

}
