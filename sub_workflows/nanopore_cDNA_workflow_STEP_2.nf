// Import Modules
include {MAKE_FAI} from '../modules/make_fai'
include {MAKE_INDEX_cDNA ; MAKE_INDEX_cDNA_CONTAMINATION_CHM13} from '../modules/make_index'
include {MAKE_INDEX_cDNA as MAKE_INDEX_CONTAMINANTS} from '../modules/make_index'
include {CHM13_GTF; CHM13_GTF_ERCC} from '../modules/chm13_gff3_to_gtf'
include {PYCHOPPER} from '../modules/pychopper'
include {PYCOQC} from '../modules/pycoqc'
include {MINIMAP2_cDNA; FILTER_BAM} from '../modules/minimap2'
include {RSEQC_GENE_BODY_COVERAGE ; RSEQC_BAM_STAT ; RSEQC_READ_GC ; CONVERT_GTF_TO_BED12 ; RSEQC_JUNCTION_ANNOTATION ; RSEQC_JUNCTION_SATURATION ; RSEQC_TIN ; RSEQC_READ_DISTRIBUTION} from '../modules/rseqc'
include {BAMBU_PREP} from '../modules/bambu'
include {MAP_CONTAMINATION_cDNA} from '../modules/contamination'
include {MAKE_CONTAMINATION_REPORT_1} from '../modules/make_contamination_report.nf'
include {MAKE_QC_REPORT_TRIM} from '../modules/num_reads_report.nf'


workflow NANOPORE_cDNA_STEP_2 {

    take:
        ref
        annotation
        housekeeping
        ont_reads_txt
        ont_reads_fq
        ercc
        cdna_kit
        track_reads
        mapq
        contamination_ref
        quality_score

    main:
        MAKE_FAI(ref)
        MAKE_INDEX_cDNA(ref)
        PYCHOPPER(ont_reads_fq, ont_reads_txt, cdna_kit, quality_score)
        MINIMAP2_cDNA(PYCHOPPER.out.id, PYCHOPPER.out.fastq,  MAKE_INDEX_cDNA.out, PYCHOPPER.out.txt, PYCHOPPER.out.num_pass_reads)
        FILTER_BAM(MINIMAP2_cDNA.out.id, mapq, MINIMAP2_cDNA.out.bam, MINIMAP2_cDNA.out.bai, MINIMAP2_cDNA.out.fastq,  MINIMAP2_cDNA.out.txt, MINIMAP2_cDNA.out.num_pass_reads)
        
        if (params.contamination_ref != "None") {

            MAKE_INDEX_CONTAMINANTS(contamination_ref)

            MAKE_INDEX_cDNA_CONTAMINATION_CHM13()

            BAM_AND_INDEX = MINIMAP2_cDNA.out.bam.combine(MAKE_INDEX_CONTAMINANTS.out).combine(MAKE_INDEX_cDNA_CONTAMINATION_CHM13.out)
            
            MAP_CONTAMINATION_cDNA(MINIMAP2_cDNA.out.id, BAM_AND_INDEX, MINIMAP2_cDNA.out.bai, MINIMAP2_cDNA.out.num_reads)

            MAKE_CONTAMINATION_REPORT_1(MAP_CONTAMINATION_cDNA.out.id, MAP_CONTAMINATION_cDNA.out.num_reads, MAP_CONTAMINATION_cDNA.out.num_unmapped_reads_before_chm13, 
            MAP_CONTAMINATION_cDNA.out.num_unmapped_reads_after_chm13, MAP_CONTAMINATION_cDNA.out.num_contaminant_reads)

        }


        if ((params.ont_reads_txt != "None") || (params.path != "None")) {
            
            PYCOQC(FILTER_BAM.out.id, FILTER_BAM.out.fastq, FILTER_BAM.out.txt, FILTER_BAM.out.bam_unfiltered, FILTER_BAM.out.bai_unfiltered, quality_score, mapq, FILTER_BAM.out.flagstat_unfiltered, 
                    FILTER_BAM.out.flagstat_filtered, FILTER_BAM.out.num_pass_reads)
            
            MAKE_QC_REPORT_TRIM(PYCOQC.out.num_reads_report, quality_score)

        }

        if (params.is_chm13 == true)
        {
            if (params.ercc == "None") 
            { 
                CHM13_GTF(annotation)
                annotation = CHM13_GTF.out.collect()
            }
            
            else 
            {
                CHM13_GTF_ERCC(annotation, ercc)
                annotation = CHM13_GTF_ERCC.out.collect()
            }
        }

        else if ((params.is_chm13 == false) && (params.housekeeping != "None"))
        {
       
            CONVERT_GTF_TO_BED12(annotation)
            RSEQC_GENE_BODY_COVERAGE(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, housekeeping)
            RSEQC_BAM_STAT(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, mapq)
            RSEQC_READ_GC(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, mapq)
            RSEQC_JUNCTION_ANNOTATION(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, CONVERT_GTF_TO_BED12.out.bed12)
            RSEQC_JUNCTION_SATURATION(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, CONVERT_GTF_TO_BED12.out.bed12)
            RSEQC_READ_DISTRIBUTION(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, CONVERT_GTF_TO_BED12.out.bed12)        
            
            if (params.rseqc_tin == true) {

                RSEQC_TIN(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, CONVERT_GTF_TO_BED12.out.bed12)
            
            } 
        
        }
        
        
        BAMBU_PREP(FILTER_BAM.out.id, mapq, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, ref, annotation, MAKE_FAI.out, track_reads)

}
