process BAMBU_PREP_DISCOVERY {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'large'

    input:
        path(bam)
        path(bai)
        path(ref)
        path(gtf)
        path(fai)
        val(NDR)

    output:
        path("bambu_prep_discovery/*.rds")

    script:
        """
        mkdir -p bambu_prep_discovery

        bambu_prep_discovery.R $bam $ref $gtf $NDR
        """
}

process BAMBU_PREP_QUANT {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'large'

    input:
        path(bam)
        path(bai)
        path(ref)
        path(gtf)
        path(fai)

    output:
        path("bambu_prep_quant/*.rds")

    script:
        """
        mkdir -p bambu_prep_quant

        bambu_prep_quant.R $bam $ref $gtf
        """
}



process BAMBU_DISCOVERY {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'medium_large'

    input:
        path(rc_files)
        path(ref)
        path(gtf)
        path(fai)
        val(NDR)
        

    output:
        path("./bambu_discovery/extended_annotations.gtf"), emit:gtf
        path("bambu_discovery/*"), emit: outty

    shell:
        '''
        mkdir bambu_discovery

        dummy="!{rc_files}"

        rc_files2="$(tr ' ' ',' <<<$dummy)"
    
        bambu_discovery.R $rc_files2 "!{ref}" "!{gtf}" "!{NDR}"
        '''
}


process BAMBU_QUANT {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'medium'

    input:
        path(rc_files)
        path(ref)
        path(gtf)
        path(fai)
        

    output:
        path("./bambu_quant/extended_annotations.gtf"), emit:gtf
        path("bambu_quant/*"), emit: outty

    shell:
        '''
        mkdir bambu_quant

        dummy="!{rc_files}"

        rc_files2="$(tr ' ' ',' <<<$dummy)"
    
        bambu_quant.R $rc_files2 "!{ref}" "!{gtf}"
        '''
}


