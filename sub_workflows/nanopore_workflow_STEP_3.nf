// Import Modules
include {BAMBU_DISCOVERY; BAMBU_QUANT; MERGE_BAMBU} from '../modules/bambu'
include {GFFCOMPARE} from '../modules/gffcompare'
include {MAKE_TRANSCRIPTOME} from '../modules/make_transcriptome'
include {MULTIQC_GRCh38} from '../modules/multiqc'
include {MAKE_CONTAMINATION_REPORT_2} from '../modules/make_contamination_report.nf'
include {MERGE_QC_REPORT} from '../modules/num_reads_report.nf'

workflow NANOPORE_STEP_3 {

    take:
        ref
        fai
        annotation
        NDR
        track_reads
        bambu_rds
        multiqc_input
        multiqc_config
        contamination
        num_reads
        read_length 
        quality_thresholds
        new_gene_and isoform_prefix

    main:
 
       
	contamination.count { n ->
	    if( n > 0 ) {
	        MAKE_CONTAMINATION_REPORT_2(contamination.collect())
	    } else {
                MAKE_CONTAMINATION_REPORT_2 = 'None'
            }
    	}

        if ((params.multiqc_input != 'None')) {
	    MERGE_QC_REPORT(num_reads.collect(), read_length.collect(), quality_thresholds.collect())
	   
	    MULTIQC_GRCh38(multiqc_input.concat(MAKE_CONTAMINATION_REPORT_2.out.flatten(), MERGE_QC_REPORT.out.flatten()).collect(), multiqc_config)
        } 

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

        }

        else
        {
            
            BAMBU_QUANT(bambu_rds, ref, annotation, fai, track_reads)

            MERGE_BAMBU(annotation, BAMBU_QUANT.out.all_transcripts_cpm.collect(), BAMBU_QUANT.out.gene_counts.collect(),
                        BAMBU_QUANT.out.all_transcripts_counts.collect(), BAMBU_QUANT.out.full_transcripts_counts.collect(),
                        BAMBU_QUANT.out.unique_transcripts_counts.collect())

        }
}
