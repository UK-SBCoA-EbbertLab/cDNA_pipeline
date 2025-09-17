process CTAT_LR_FUSION {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'ctat_lr_fusion'

    input:
        val(id)
	path(fastq)
        path(ctat_lib_dir)

    output:
        path("ctat_LR_fusion_outdir_${id}/*")
	path("${id}_ctat-LR-fusion.fusion_predictions.abridged.tsv")
	path("${id}_ctat-LR-fusion.fusion_predictions.tsv")
	path("${id}_ctat-LR-fusion.fusion_inspector_web.html")

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
