process MAP_CONTAMINATION_cDNA {

    publishDir "results/${params.out_dir}/contamination_report/aligment/${id}", mode: "copy", pattern: "*"

    label 'bambu_prep_job'

    input:
        val(id)
        tuple path(bam), path(index)
        path(bai)
        val(num_reads)

    output:
        val("$id"), emit: id
        val("$num_reads"), emit: num_reads
        val("!NUM_UNMAPPED_READS"), emit: num_unmapped_reads
        val("!NUM_CONTAMINANT_READS"), emit: num_contaminant_reads
        path("${id}*"), emit: outty

    script:
        """

        samtools view -h -b -f 4 "${bam}" > "${id}_unmapped.bam"
        samtools fastq "${id}_unmapped.bam" > "${id}_unmapped_reads.fastq"

        NUM_UNMAPPED_READS=\$(samtools view -F 0x40 "${id}_unmapped.bam" | cut -f1 | sort | uniq | wc -l)
        
        minimap2 -t 12 -ax splice \
            --split-prefix /tmp/tmp_name \
            -uf \
            $index \
            "${id}_unmapped_reads.fastq" > "${id}_contaminants_unsorted.bam" 

        samtools view -b -@ 12 -F 256 "${id}_contaminants_unsorted.bam" > "${id}_contaminants_unsorted_primary.bam"
        samtools sort -@ -12 "${id}_contaminants_unsorted_primary.bam" -o "${id}_contaminants.bam"
        samtools index "${id}_contaminants.bam"

        samtools idxstat "${id}_contaminants.bam" > "${id}_contaminants.idxstat"

        NUM_CONTAMINANT_READS=\$(samtools view -F 0x40 "${id}_contaminants.bam" | cut -f1 | sort | uniq | wc -l)

        rm "${id}_contaminants_unsorted.bam" "${id}_unmapped.bam"
        """

}

process MAP_CONTAMINATION_dRNA {

    publishDir "results/${params.out_dir}/contamination_report/aligment/${id}", mode: "copy", pattern: "*"

    label 'bambu_prep_job'

    input:
        val(id)
        tuple path(bam), path(index)
        path(bai)
        val(num_reads)

    output:
        val("$id"), emit: id
        val("$num_reads"), emit: num_reads
        val("!NUM_UNMAPPED_READS"), emit: num_unmapped_reads
        val("!NUM_CONTAMINANT_READS"), emit: num_contaminant_reads
        path("${id}*"), emit: outty

    script:
        """

        samtools view -h -b -f 4 "${bam}" > "${id}_unmapped.bam"
        samtools fastq "${id}_unmapped.bam" > "${id}_unmapped_reads.fastq"

        NUM_UNMAPPED_READS=\$(samtools view -F 0x40 "${id}_unmapped.bam" | cut -f1 | sort | uniq | wc -l)
        
        minimap2 -t 12 -ax splice \
            --split-prefix /tmp/tmp_name \
            -k14 -uf \
            $index \
            "${id}_unmapped_reads.fastq" > "${id}_contaminants_unsorted.bam" 

        samtools view -b -@ 12 -F 256 "${id}_contaminants_unsorted.bam" > "${id}_contaminants_unsorted_primary.bam"
        samtools sort -@ -12 "${id}_contaminants_unsorted_primary.bam" -o "${id}_contaminants.bam"
        samtools index "${id}_contaminants.bam"

        samtools idxstat "${id}_contaminants.bam" > "${id}_contaminants.idxstat"

        NUM_CONTAMINANT_READS=\$(samtools view -F 0x40 "${id}_contaminants.bam" | cut -f1 | sort | uniq | wc -l)

        rm "${id}_contaminants_unsorted.bam" "${id}_unmapped.bam"
        """

}
