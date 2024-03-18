// Import Modules
include {MAKE_FAI} from '../modules/make_fai'
include {MAKE_INDEX_dRNA ; MAKE_INDEX_dRNA_CONTAMINATION_CHM13} from '../modules/make_index'
include {MAKE_INDEX_dRNA as MAKE_INDEX_CONTAMINANTS} from '../modules/make_index'
include {CHM13_GTF; CHM13_GTF_ERCC} from '../modules/chm13_gff3_to_gtf'
include {PYCOQC_dRNA} from '../modules/pycoqc'
include {MINIMAP2_dRNA; FILTER_BAM} from '../modules/minimap2'
include {RSEQC_GENE_BODY_COVERAGE ; RSEQC_BAM_STAT ; RSEQC_READ_GC ; CONVERT_GTF_TO_BED12 ; RSEQC_JUNCTION_ANNOTATION ; RSEQC_JUNCTION_SATURATION ; RSEQC_TIN ; RSEQC_READ_DISTRIBUTION} from '../modules/rseqc'
include {BAMBU_PREP} from '../modules/bambu'
include {MAP_CONTAMINATION_dRNA} from '../modules/contamination'
include {MAKE_CONTAMINATION_REPORT_1} from '../modules/make_contamination_report.nf'
include {TRIM_dRNA} from '../modules/trim_dRNA.nf'
include {CONVERT_U_TO_T} from '../modules/convert_U_to_T.nf'
include {MAKE_QC_REPORT_TRIM ; MAKE_QC_REPORT_NO_TRIM} from '../modules/num_reads_report.nf'

workflow NANOPORE_dRNA_STEP_2 {

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
        MAKE_INDEX_dRNA(ref)
        CONVERT_U_TO_T(ont_reads_fq, ont_reads_txt, quality_score)
        

        if (params.trim_dRNA == true) { 
            
            TRIM_dRNA(CONVERT_U_TO_T.out.fastq, CONVERT_U_TO_T.out.txt, CONVERT_U_TO_T.out.num_pass_reads)
            MINIMAP2_dRNA(TRIM_dRNA.out.fastq,  MAKE_INDEX_dRNA.out, TRIM_dRNA.out.txt, TRIM_dRNA.out.num_pass_reads)

        } else {
    
            MINIMAP2_dRNA(CONVERT_U_TO_T.out.fastq,  MAKE_INDEX_dRNA.out, CONVERT_U_TO_T.out.txt, CONVERT_U_TO_T.out.num_pass_reads)
    
        }
       

        FILTER_BAM(MINIMAP2_dRNA.out.id, mapq, MINIMAP2_dRNA.out.bam, MINIMAP2_dRNA.out.bai, MINIMAP2_dRNA.out.fastq,  MINIMAP2_dRNA.out.txt, MINIMAP2_dRNA.out.num_pass_reads)
        

        if (params.contamination_ref != "None") {
      
            MAKE_INDEX_CONTAMINANTS(contamination_ref)

            MAKE_INDEX_dRNA_CONTAMINATION_CHM13()
                  
            BAM_AND_INDEX = MINIMAP2_dRNA.out.bam.combine(MAKE_INDEX_CONTAMINANTS.out).combine(MAKE_INDEX_dRNA_CONTAMINATION_CHM13.out)
                              
            MAP_CONTAMINATION_dRNA(MINIMAP2_dRNA.out.id, BAM_AND_INDEX, MINIMAP2_dRNA.out.bai, MINIMAP2_dRNA.out.num_reads)
        
            MAKE_CONTAMINATION_REPORT_1(MAP_CONTAMINATION_dRNA.out.id, MAP_CONTAMINATION_dRNA.out.num_reads, MAP_CONTAMINATION_dRNA.out.num_unmapped_reads_before_chm13,
            MAP_CONTAMINATION_dRNA.out.num_unmapped_reads_after_chm13, MAP_CONTAMINATION_dRNA.out.num_contaminant_reads)
        
        }
        
        if ((params.ont_reads_txt != "None") || (params.path != "None")) {
           
            PYCOQC_dRNA(FILTER_BAM.out.id, FILTER_BAM.out.fastq, FILTER_BAM.out.txt, FILTER_BAM.out.bam_unfiltered, FILTER_BAM.out.bai_unfiltered,
                         quality_score, mapq, FILTER_BAM.out.flagstat_unfiltered, FILTER_BAM.out.flagstat_filtered, FILTER_BAM.out.num_pass_reads)

            if (params.trim_dRNA == true) {
            
                MAKE_QC_REPORT_TRIM(PYCOQC_dRNA.out.num_reads_report, quality_score)

            } else {
                
                MAKE_QC_REPORT_NO_TRIM(PYCOQC_dRNA.out.num_reads_report, quality_score)

            } 

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
            RSEQC_TIN(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, CONVERT_GTF_TO_BED12.out.bed12)
            RSEQC_READ_DISTRIBUTION(FILTER_BAM.out.id, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, CONVERT_GTF_TO_BED12.out.bed12)
        }
       
        BAMBU_PREP(FILTER_BAM.out.id, mapq, FILTER_BAM.out.bam_filtered, FILTER_BAM.out.bai_filtered, ref, annotation, MAKE_FAI.out, track_reads)

}
