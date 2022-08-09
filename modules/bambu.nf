process BAMBU_PREP {

    label 'large'

    input:
        path(bam)
        path(bai)
        path(ref)
        path(gtf)
        path(fai)

    output:
        path("bambu_prep/*.rds")

    script:
        """
        mkdir -p bambu_prep

        bambu_prep.R $bam $ref $gtf
        """
}

process BAMBU_DISCOVERY {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'medium'

    input:
        path(rc_files)
        path(ref)
        path(gtf)
        path(fai)
        

    output:
        path("./bambu_discovery/extended_annotations.gtf"), emit:gtf
        path("bambu_discovery/*"), emit: outty

    shell:
        '''
        mkdir bambu_discovery

        dummy="!{rc_files}"

        rc_files2="$(tr ' ' ',' <<<$dummy)"
    
        bambu_discovery.R $rc_files2 "!{ref}" "!{gtf}"
        '''
}

