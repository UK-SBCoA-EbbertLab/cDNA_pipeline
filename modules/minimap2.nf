process MINIMAP2_cDNA {

    publishDir "results/${params.out_dir}/mapping_cDNA/", pattern: "*.ba*", mode: "copy", overwrite: true

    label 'large'

    input:
        val(id)
        path(fastq)
        file(index)
        val(txt)
        val(num_pass_reads)

    output:
        val("$id"), emit: id
        path("$fastq"), emit: fastq
        path("${id}.bam"), emit: bam
        path("${id}.bam.bai"), emit: bai
        val("$txt"), emit: txt
        env(NUM_READS), emit: num_reads
        val("${num_pass_reads}"), emit: num_pass_reads

    script:
        """
        minimap2 -t 50 -ax splice \
            -uf \
            $index \
            $fastq > "${id}_all.bam" \
     

        samtools sort -@ -12 "${id}_all.bam" -o "${id}.bam" 
        samtools index "${id}.bam"

        NUM_READS=\$(samtools view -F 0x40 "${id}.bam" | cut -f1 | sort | uniq | wc -l)

        rm "${id}_all.bam"
        """

}

process MINIMAP2_dRNA {

    publishDir "results/${params.out_dir}/mapping_dRNA/", pattern: "*.ba*", mode: "copy", overwrite: true

    label 'large'

    input:
        tuple val(id), path(fastq)
        file(index)
        val(txt)
        val(num_pass_reads)

    output:
        val("$id"), emit: id
        path("$fastq"), emit: fastq
        path("${id}.bam"), emit: bam
        path("${id}.bam.bai"), emit: bai
        val("${txt}"), emit: txt
        env(NUM_READS), emit: num_reads
        val("${num_pass_reads}"), emit: num_pass_reads

    script:
        """
        
        minimap2 -t 50 -ax splice \
            -k14 -uf \
            $index \
            $fastq > "${id}_all.bam" \


        samtools sort -@ -12 "${id}_all.bam" -o "${id}.bam"
        samtools index "${id}.bam"

        NUM_READS=\$(samtools view -F 0x40 "${id}.bam" | cut -f1 | sort | uniq | wc -l)

        rm "${id}_all.bam"
        """

}


process FILTER_BAM {

    publishDir "results/${params.out_dir}/bam_filtering/", mode: "copy", pattern: "*filtered_mapq*.*stat", overwrite: true
    publishDir "results/${params.out_dir}/multiQC_input/minimap2/", mode: "copy", pattern: "*unfiltered.*stat", overwrite: true

   
    label 'medium_small'

    input:
        val(id)
        val(mapq)
        path(bam)
        path(bai)
        path(fastq)
        val(txt)
        val(num_pass_reads)

    output:
        val("$id"), emit: id
        path("${id}_filtered_mapq_${mapq}.bam"), emit: bam_filtered
        path("${id}_filtered_mapq_${mapq}.bam.bai"), emit: bai_filtered
        path("$bam"), emit: bam_unfiltered
        path("$bai"), emit: bai_unfiltered
        path("*.*stat"), emit: QC
        path("*unfiltered.flagstat"), emit: flagstat_unfiltered 
        path("*filtered_mapq*.flagstat"), emit: flagstat_filtered
        path("$fastq"), emit: fastq
        val("$txt"), emit: txt
        val("${num_pass_reads}"), emit: num_pass_reads
        

    script:
        """
 
        samtools flagstat $bam > "${id}_unfiltered.flagstat"
        samtools idxstats $bam > "${id}_unfiltered.idxstat"

       
        samtools view -b -q $mapq -F 2304 -@ 12 $bam > 'intermediate.bam'
        samtools sort -@ 12 "intermediate.bam" -o '${id}_filtered_mapq_${mapq}.bam'
        samtools index '${id}_filtered_mapq_${mapq}.bam'
        samtools flagstat "${id}_filtered_mapq_${mapq}.bam" > "${id}_filtered_mapq_${mapq}.flagstat"
        samtools idxstats "${id}_filtered_mapq_${mapq}.bam" > "${id}_filtered_mapq_${mapq}.idxstat"

        rm "intermediate.bam"
        """

}



process FILTER_BAM_ONLY {

    publishDir "results/${params.out_dir}/bam_filtering/", mode: "copy", pattern: "*filtered_mapq*.*stat", overwrite: true


    label 'medium_small'

    input:
        tuple val(id), path(bam)
        val(bai)
        val(mapq)

    output:
        val("$id"), emit: id
        path("${id}_filtered_mapq_${mapq}.bam"), emit: bam
        path("${id}_filtered_mapq_${mapq}.bam.bai"), emit: bai
        path("*filtered_mapq*.*stat"), emit: QC

    script:
        """

        samtools view -b -q $mapq -F 2304 -@ 12 $bam > 'intermediate.bam'
        samtools sort -@ 12 "intermediate.bam" -o '${id}_filtered_mapq_${mapq}.bam'
        samtools index '${id}_filtered_mapq_${mapq}.bam'
        samtools flagstat "${id}_filtered_mapq_${mapq}.bam" > "${id}_filtered_mapq_${mapq}.flagstat"
        samtools idxstats "${id}_filtered_mapq_${mapq}.bam" > "${id}_filtered_mapq_${mapq}.idxstat"

        rm "intermediate.bam"
        """

}
