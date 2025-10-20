process GFFCOMPARE {

    publishDir "results/${params.out_dir}/gffcompare/", mode: "copy", overwrite: true

    label 'small'

    input:
        path(extended_annotation)
        path(reference_annotation)

    output:
        path "*"
	path "*.tracking", emit: tracking_file

    script:
        """
        gffcompare -r $reference_annotation $extended_annotation -o "gffcompare_output"
        """
}
