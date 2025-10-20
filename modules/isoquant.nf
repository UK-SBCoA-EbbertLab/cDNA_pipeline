process ISOQUANT {

    publishDir "results/${params.out_dir}/", mode: "copy", overwrite: true

    label 'isoquant'

    input:
	tuple val(chr), path(bams), path(bais)
        path(ref)
        path(gtf)

    output:
        path("IsoQuant_${chr}/*")
	path("IsoQuant_${chr}/**/*.discovered_transcript_grouped_counts.tsv"), emit: discovery_raw_transcript
        path("IsoQuant_${chr}/**/*.transcript_grouped_counts.tsv"), emit: quantification_raw_transcript
        path("IsoQuant_${chr}/**/*.discovered_gene_grouped_counts.tsv"), emit: discovery_raw_gene
        path("IsoQuant_${chr}/**/*.gene_grouped_counts.tsv"), emit: quantification_raw_gene
        path("IsoQuant_${chr}/**/*.extended_annotation.gtf"), emit: extended_annotation

    script:

	
        """
	export HOME=\${PWD}

	isoquant.py \
		--reference $ref \
		--genedb $gtf \
		--complete_genedb \
		--data_type ont \
		--bam $bams \
		--output IsoQuant \
		--prefix chr_${chr} \
		--count_exons \
		--threads 16 \
		--normalization_method none \
		--counts_format matrix

        """
}

process SEPARATE_CHROMOSOMES {

    label 'medium_large'

    input:
	path(ref)
	tuple val(name), path(bams_and_bais)

    output:
	path("*.bam"), emit: chr_bams
	path("*.bai"), emit: chr_bais

    script:
    """
	chroms=\$(seq 1 22)
	chroms="\${chroms} MT"

	#split the bam files, add chromosome to the beginning of the name
	for chr in \${chroms}; do
	    echo "[INFO] Extracting \${chr}"
	    samtools view -b ${name}.bam \${chr} > \${chr}_${name}.bam
	    samtools index \${chr}_${name}.bam
	done

	echo "[INFO] Extracting all other chromosomes and contigs"
	other=\$(samtools idxstats ${name}.bam | awk -vOFS='\t' '\$1 !~ /^([0-9]+|MT|chr[0-9]+|chrMT|\\*)\$/' | cut -f1)
	if [[ -n "\${other}" ]]; then
	    samtools view -b ${name}.bam \${other} > "others_${name}.bam"
	    samtools index "others_${name}.bam"
	else 
	    echo "[INFO] Missing some contigs."
	fi

    """
}

process CONCAT_CHROM_MATRICES {

    publishDir "results/${params.out_dir}/IsoQuant_merged", mode: "copy", overwrite: true

    label 'medium_large'

    input:
	path(discovery_transcripts)
	path(discovery_gene)
	path(quantification_transcripts)
	path(quantification_gene)
	path(annotation)
	path(gtfs)

    output:
	path("*.tsv"), emit: counts_files
	path("*discovered_transcript_grouped_counts.tsv"), emit: combined_discovery_transcript
	path("*discovered_gene_grouped_counts.tsv"), emit: combined_discovery_gene
	path("*transcript_grouped_counts.tsv"), emit: combined_quant_transcript
	path("*gene_grouped_counts.tsv"), emit: combined_quant_gene
	path("*extended_annotation.gtf"), emit: combined_ext_annotation

    script:
    """
	concat_isoquant_extended_annotations.py --out_gtf "isoquant_combined.extended_annotation.gtf" --gtfs ${gtfs} --annotation ${annotation}
	concat_isoquant_matrices.py --counts_matrices ${discovery_gene} --counts_type "gene" --out "isoquant_combined.discovered_gene_grouped_counts.tsv"
	concat_isoquant_matrices.py --counts_matrices ${quantification_gene} --counts_type "gene" --out "isoquant_combined.gene_grouped_counts.tsv"
	concat_isoquant_matrices.py --counts_matrices ${discovery_transcripts} --counts_type "transcript" --out "isoquant_combined.discovered_transcript_grouped_counts.tsv" --gtf "isoquant_combined.extended_annotation.gtf"
	concat_isoquant_matrices.py --counts_matrices ${quantification_transcripts} --counts_type "transcript" --out "isoquant_combined.transcript_grouped_counts.tsv" --gtf "isoquant_combined.extended_annotation.gtf"

    """
}

//process COMPARE_BAMBU_AND_ISOQUANT {

//    publishDir "results/${params.out_dir}/bambu_isoquant_comparison", mode: "copy", overwrite: true

//    label 'bambu_isoquant_compare'

//    input:
//	path(bambu_counts)
//	path(isoquant_counts)
//	path(gff_compare)

//    output:
//	path(*)

//    script:
//    """
//	compare_bambu_isoquant.R ${bambu_counts} ${isoquant_counts} ${gff_compare}
//    """

//}
