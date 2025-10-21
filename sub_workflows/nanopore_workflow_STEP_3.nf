// Import Modules
include {BAMBU_DISCOVERY; BAMBU_QUANT; MERGE_BAMBU} from '../modules/bambu'
include {GFFCOMPARE ; GFFCOMPARE as GFFCOMPARE_QUERY} from '../modules/gffcompare'
include {MAKE_TRANSCRIPTOME} from '../modules/make_transcriptome'
include {MULTIQC_GRCh38} from '../modules/multiqc'
include {MAKE_CONTAMINATION_REPORT_2} from '../modules/make_contamination_report.nf'
include {MERGE_QC_REPORT} from '../modules/num_reads_report.nf'
include {ISOQUANT; SEPARATE_CHROMOSOMES ; CONCAT_CHROM_MATRICES} from '../modules/isoquant.nf'
include {NORMALIZE_COUNTS; NORMALIZE_COUNTS as NORMALIZE_BAMBU_COUNTS} from '../modules/counts_normalization.nf'

workflow NANOPORE_STEP_3 {

    take:
        ref
        fai
        annotation
        NDR
        track_reads
        bambu_rds
	filtered_bam
        multiqc_input
        multiqc_config
        contamination
        num_reads
        read_length 
        quality_thresholds
        new_gene_and_isoform_prefix
	quantification_tool

    main:

	contamination.view() 
      
 
        if ((params.multiqc_input != 'None')) {
	    MERGE_QC_REPORT(num_reads.collect(), read_length.collect(), quality_thresholds.collect())
        } 

	MAKE_CONTAMINATION_REPORT_2(contamination.collect())
	MULTIQC_GRCh38(multiqc_input.concat(MAKE_CONTAMINATION_REPORT_2.out.flatten(), MERGE_QC_REPORT.out.flatten()).collect(), multiqc_config)


	if (quantification_tool != "isoquant") {
            if (params.is_discovery == true)
            {

                BAMBU_DISCOVERY(bambu_rds.collect(), ref, annotation, fai, NDR, new_gene_and_isoform_prefix)
            
                BAMBU_QUANT(bambu_rds, ref, BAMBU_DISCOVERY.out.gtf, fai, track_reads)
            
                MERGE_BAMBU(BAMBU_DISCOVERY.out.gtf, BAMBU_QUANT.out.all_transcripts_cpm.collect(), BAMBU_QUANT.out.gene_counts.collect(),
                        BAMBU_QUANT.out.all_transcripts_counts.collect(), BAMBU_QUANT.out.full_transcripts_counts.collect(),
                        BAMBU_QUANT.out.unique_transcripts_counts.collect())

                new_annotation = BAMBU_DISCOVERY.out.gtf
                GFFCOMPARE(new_annotation, annotation)
                MAKE_TRANSCRIPTOME(ref, fai, new_annotation)

		MERGE_BAMBU.out.outty
			.set { bambu_counts }

            }

            else
            {
            
                BAMBU_QUANT(bambu_rds, ref, annotation, fai, track_reads)
		BAMBU_QUANT.out.counts_files
			.set { bambu_counts }

                MERGE_BAMBU(annotation, BAMBU_QUANT.out.all_transcripts_cpm.collect(), BAMBU_QUANT.out.gene_counts.collect(),
                        BAMBU_QUANT.out.all_transcripts_counts.collect(), BAMBU_QUANT.out.full_transcripts_counts.collect(),
                        BAMBU_QUANT.out.unique_transcripts_counts.collect())
		
                MERGE_BAMBU.out.outty
			.set { bambu_counts }
            }

	    NORMALIZE_BAMBU_COUNTS(bambu_counts)
	}
	
	if (quantification_tool != "bambu")
	{
	    
	    SEPARATE_CHROMOSOMES(ref, filtered_bam)
	
	    SEPARATE_CHROMOSOMES.out.chr_bams
		.collect()
		.flatten()
		.map { bam ->
		    def name = bam.getBaseName()
		    def (chr, sample) = name.split('_', 2)
		    tuple(chr, bam)
		}
		.groupTuple()
		.set { chrom_bams }

	    SEPARATE_CHROMOSOMES.out.chr_bais
		.collect()
		.flatten()
		.map { bai ->
		    def name = bai.getBaseName()
		    def (chr, sample) = name.split('_', 2)
		    tuple(chr, bai)
		}
		.groupTuple()
		.set { chrom_bais }

	    chrom_bams
		.join(chrom_bais)
	        .set{ chrom_tuple }

	    ISOQUANT(chrom_tuple, ref, annotation)
	    CONCAT_CHROM_MATRICES(ISOQUANT.out.discovery_raw_transcript.collect(), ISOQUANT.out.discovery_raw_gene.collect(), ISOQUANT.out.quantification_raw_transcript.collect(), ISOQUANT.out.quantification_raw_gene.collect(), annotation, ISOQUANT.out.extended_annotation.collect())

	    NORMALIZE_COUNTS(CONCAT_CHROM_MATRICES.out.counts_files)
	}


	if (quantification_tool == "both") {
	    GFFCOMPARE_QUERY(new_annotation, CONCAT_CHROM_MATRICES.out.combined_ext_annotation)
	    //COMAPARE_BAMBU_ISOQUANT(NORMALIZE_BAMBU_COUNTS.out.raw_counts, NORMALZE_COUNTS.out.raw_counts, GFFCOMPARE_QUERY.out.tracking_file)

	}
        
}
