process MINIMAP2_cDNA {

    publishDir "results/${params.out_dir}/mapping_cDNA/"

    label 'large'

    input:
        val(id)
        path(fastq)
        path(index)

    output:
        val("$id"), emit: id
        path("${id}_all_sorted.bam"), emit: bam_all
        path("${id}_all_sorted.bam.bai"), emit: bai_all
        path("${id}_mapped_filtered_sorted.bam"), emit: bam_mapped
        path("${id}_mapped_filtered_sorted.bam.bai"), emit: bai_mapped

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
    
        samtools view -b -q 10 -F 2304 -@ 12 '${id}_all.bam' > '${id}_mapped_filtered.bam'
        samtools sort -@ 12 "${id}_mapped_filtered.bam" -o '${id}_mapped_filtered_sorted.bam'
        samtools index '${id}_mapped_filtered_sorted.bam'
        samtools flagstat "${id}_mapped_filtered_sorted.bam" > "${id}_mapped_filtered_sorted.flagstat"
        samtools idxstats "${id}_mapped_filtered_sorted.bam" > "${id}_mapped_filtered_sorted.idxstat"
        
        rm "${id}_all.bam"
        """

}


process MINIMAP2_QC {

    label 'large'

    input:
        val(id)
        path(fastq)
        path(txt)
        path(index)

    output:
        val("$id"), emit: id
        path("$txt"), emit: txt
        path("${id}_all_sorted.bam"), emit: bam
        path("${id}_all_sorted.bam.bai"), emit: bai

    script:
        """
        minimap2 -t 16 -ax splice \
            -uf \
            $index \
            $fastq > "${id}_all.bam" \


        samtools sort -@ -12 "${id}_all.bam" -o "${id}_all_sorted.bam"
        samtools index "${id}_all_sorted.bam"

        rm "${id}_all.bam"
        """
}

