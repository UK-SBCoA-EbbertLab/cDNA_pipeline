process ISOQUANT {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'isoquant'

    input:
	val(mapq)
        path(bams)
        path(bais)
        path(ref)
        path(gtf)
        path(fai)

    output:
        path("IsoQuant/*")

    script:

	
        """
	export HOME=\${PWD}

	isoquant.py \
		--reference $ref \
		--genedb $gtf \
		--complete_genedb \
		--data_type ont \
		--bam $bams \
		--output IsoQuant \
		--prefix mapq_${mapq} \
		--count_exons \
		--threads 16 \
		--normalization_method none \
		--counts_format matrix

        """
}
