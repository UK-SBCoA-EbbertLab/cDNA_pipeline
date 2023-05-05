process MINIMAP2_cDNA {

    publishDir "results/${params.out_dir}/mapping_cDNA/", pattern: "*sorted*.ba*", mode: "copy", overwrite: true
    publishDir "results/${params.out_dir}/multiQC_input/minimap2/", pattern: "*sorted.*stat", mode: "copy", overwrite: true

    label 'large'

    input:
        val(id)
        path(fastq)
        path(index)
        path(txt)

    output:
        val("$id"), emit: id
        path("$fastq"), emit: fastq
        path("${id}_all_sorted.bam"), emit: bam
        path("${id}_all_sorted.bam.bai"), emit: bai
        path("*all*stat"), emit: QC_out
        path("$txt"), emit: txt

    script:
        """
        minimap2 -t 16 -ax splice \
            -uf \
            $index \
            $fastq > "${id}_all.bam" \
     

        samtools sort -@ -12 "${id}_all.bam" -o "${id}_all_sorted.bam" 
        samtools index "${id}_all_sorted.bam"
        samtools flagstat "${id}_all_sorted.bam" > "${id}_all_sorted.flagstat"
        samtools idxstats "${id}_all_sorted.bam" > "${id}_all_sorted.idxstat"
    
        """

}

process FILTER_BAM {

publishDir "results/${params.out_dir}/QC/filtering_bam/", pattern: "*sorted.*stat"
    
    label 'medium_small'

    input:
        val(id)
        val(mapq)
        path(bam)
        path(bai)

    output:
        val("$id"), emit: id
        path("${id}_filtered_mapq_${mapq}_sorted.bam"), emit: bam
        path("${id}_filtered_mapq_${mapq}_sorted.bam.bai"), emit: bai
        path("*sorted.*stat"), emit: QC

    script:
        """
        
        samtools view -b -q $mapq -F 2304 -@ 12 $bam > 'intermediate.bam'
        samtools sort -@ 12 "intermediate.bam" -o '${id}_filtered_mapq_${mapq}_sorted.bam'
        samtools index '${id}_filtered_mapq_${mapq}_sorted.bam'
        samtools flagstat "${id}_filtered_mapq_${mapq}_sorted.bam" > "${id}_filtered_mapq_${mapq}_sorted.flagstat"
        samtools idxstats "${id}_filtered_mapq_${mapq}_sorted.bam" > "${id}_filtered_mapq_${mapq}_sorted.idxstat"

        rm "intermediate.bam"
        """

}
