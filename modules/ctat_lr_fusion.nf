process CTAT_LR_FUSION {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'ctat_lr_fusion'

    input:
        val(id)
	path(fastq)
        path(ctat_lib_dir)

    output:
        path("ctat_LR_fusion_outdir_${id}/*")
	path("${id}_ctat-LR-fusion.fusion_predictions.abridged.tsv"), emit:abridged_tsv
	path("${id}_ctat-LR-fusion.fusion_predictions.tsv"), emit:predicitons_tsv
	path("${id}_ctat-LR-fusion.fusion_inspector_web.html"), emit: html

    script:

	
        """
	ctat-LR-fusion -T ${fastq} \
		--genome_lib_dir ${ctat_lib_dir} \
		--CPU 8 \
		--vis \
		--output ctat_LR_fusion_outdir_${id}

	mv ctat_LR_fusion_outdir_${id}/ctat-LR-fusion.fusion_predictions.abridged.tsv ${id}_ctat-LR-fusion.fusion_predictions.abridged.tsv
	mv ctat_LR_fusion_outdir_${id}/ctat-LR-fusion.fusion_predictions.tsv ${id}_ctat-LR-fusion.fusion_predictions.tsv
	mv ctat_LR_fusion_outdir_${id}/ctat-LR-fusion.fusion_inspector_web.html ${id}_ctat-LR-fusion.fusion_inspector_web.html

        """
}

process COMBINE_FUSION_OUTPUT_TSVS {

    publishDir "results/${params.out_dir}/ctat_LR_fusion", mode: "copy", overwrite: true

    label 'ctat_lr_fusion'

    input:
	path(tsvs)

    output:
	path("combined_ctat-LR-fusion.fusion_predictions.abridged.tsv")

    script:

	"""
	ctat_lr_fusion_create_matrix.py
	"""
}
	

