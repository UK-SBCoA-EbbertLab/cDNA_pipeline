process BAMBU_PREP {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'bambu_prep_job'

    input:
        val(id)
        val(mapq)
        path(bam)
        path(bai)
        path(ref)
        path(gtf)
        path(fai)
        val(track_reads)

    output:
        path("bambu_prep/*.rds")

    script:
        """
        mkdir -p bambu_prep

        bambu_prep.R $bam $ref $gtf $track_reads

        mv ./bambu_prep/*.rds "./bambu_prep/${id}_mapq_${mapq}.rds"
        """
}

process BAMBU_DISCOVERY {


    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'bambu_discovery'

    input:
        path(rc_files)
        path(ref)
        path(gtf)
        path(fai)
        val(NDR)
        val(new_gene_and_isoform_prefix)
        

    output:
        path("./bambu_discovery/extended_annotations.gtf"), emit:gtf
        path("bambu_discovery/*"), emit: outty
	path("bambu_discovery/*.tsv"), emit: counts_files

    shell:
        '''
        mkdir bambu_discovery

        dummy="!{rc_files}"

        rc_files2="$(tr ' ' ',' <<<$dummy)"
    
        bambu_discovery_low_memory.R $rc_files2 "!{ref}" "!{gtf}" "!{NDR}" "!{new_gene_and_isoform_prefix}" "bambu_discovery"
        '''
}


process BAMBU_QUANT {

    publishDir "results/${params.out_dir}/bambu_quant/", mode: "copy", overwrite: true

    label 'bambu_prep_job'

    input:
        path(rc_file)
        path(ref)
        path(gtf)
        path(fai)
        val(track_reads)
        

    output:
        path("*CPM_transcript.txt"), emit: all_transcripts_cpm
        path("*counts_gene.txt"), emit: gene_counts
        path("*counts_transcript.txt"), emit: all_transcripts_counts
        path("*.rds"), emit: rds
        path("*fullLengthCounts_transcript.txt"), emit: full_transcripts_counts
        path("*uniqueCounts_transcript.txt"), emit: unique_transcripts_counts
        path("*.txt"), emit: counts_files

    script:
        """

        bambu_quant.R $rc_file $ref $gtf $track_reads
       
        """
}

process MERGE_BAMBU {

    publishDir "results/${params.out_dir}/bambu_final_files/", mode: "copy", overwrite: true

    label "large"

    input:
        path(gtf)
        path(all_transcript_cpm)
        path(gene_counts)
        path(all_transcript_counts)
        path(full_transcript_counts)
        path(unique_transcript_counts)
    
    output:
        path("*.txt"), emit: outty
        path("$gtf"), emit: gtf

    script:
        """
        
        merge_counts_matrices.py $all_transcript_cpm
       
        merge_counts_matrices.py $gene_counts
       
        merge_counts_matrices.py $all_transcript_counts
       
        merge_counts_matrices.py $full_transcript_counts
       
        merge_counts_matrices.py $unique_transcript_counts

        """

}
