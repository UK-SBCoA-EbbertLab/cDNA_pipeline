process NORMALIZE_COUNTS {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'medium'

    input:
	path(counts_file)

    output:
        path("normalization/*")
	path("normalization/*_tmm.tsv"), emit: tmm_normalized
	path("normalization/*_counts.tsv"), emit: raw_counts


    script:

	
        """

	normalize_counts.R ${counts_file} .

        """
}

