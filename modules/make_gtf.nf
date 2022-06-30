process MAKE_GTF {

    publishDir "results/${params.out_dir}/annotations/", mode: "copy", overwrite: true

    label 'medium'

    input:
        path annotation

    output:
        path '*.gtf'

    script:
        """
        chm13_gff3_to_gtf.py $annotation "CHM13_v2.0.gtf"
        """
}
