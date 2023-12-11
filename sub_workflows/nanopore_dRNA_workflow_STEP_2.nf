// Import Modules
include {MAKE_FAI} from '../modules/make_fai'
include {MAKE_INDEX_dRNA ; MAKE_INDEX_dRNA_CONTAMINATION_CHM13} from '../modules/make_index'
include {MAKE_INDEX_dRNA as MAKE_INDEX_CONTAMINANTS} from '../modules/make_index'
include {CHM13_GTF; CHM13_GTF_ERCC} from '../modules/chm13_gff3_to_gtf'
include {PYCOQC_dRNA} from '../modules/pycoqc'
include {MINIMAP2_dRNA; FILTER_BAM} from '../modules/minimap2'
include {RSEQC} from '../modules/rseqc'
include {BAMBU_PREP} from '../modules/bambu'
include {MAP_CONTAMINATION_dRNA} from '../modules/contamination'
include {MAKE_CONTAMINATION_REPORT_1 ; MAKE_CONTAMINATION_REPORT_2} from '../modules/make_contamination_report.nf'


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

    main:
        MAKE_FAI(ref)
        MAKE_INDEX_dRNA(ref)

        MINIMAP2_dRNA(ont_reads_fq,  MAKE_INDEX_dRNA.out, ont_reads_txt)
        FILTER_BAM(MINIMAP2_dRNA.out.id, mapq, MINIMAP2_dRNA.out.bam, MINIMAP2_dRNA.out.bai)
        

        if (params.contamination_ref != "None") {
      
            MAKE_INDEX_CONTAMINANTS(contamination_ref)

            MAKE_INDEX_dRNA_CONTAMINATION_CHM13()
                  
            BAM_AND_INDEX = MINIMAP2_dRNA.out.bam.combine(MAKE_INDEX_CONTAMINANTS.out).combine(MAKE_INDEX_dRNA_CONTAMINATION_CHM13.out)
                              
            MAP_CONTAMINATION_dRNA(MINIMAP2_dRNA.out.id, BAM_AND_INDEX, MINIMAP2_dRNA.out.bai, MINIMAP2_dRNA.out.num_reads)
        
            MAKE_CONTAMINATION_REPORT_1(MAP_CONTAMINATION_dRNA.out.id, MAP_CONTAMINATION_dRNA.out.num_reads, MAP_CONTAMINATION_dRNA.out.num_unmapped_reads_before_chm13,
            MAP_CONTAMINATION_dRNA.out.num_unmapped_reads_after_chm13, MAP_CONTAMINATION_dRNA.out.num_contaminant_reads)

            MAKE_CONTAMINATION_REPORT_2(MAKE_CONTAMINATION_REPORT_1.out.collect())
        }
        
        if (params.ont_reads_txt != "None") {
            PYCOQC_dRNA(MINIMAP2_dRNA.out.id, MINIMAP2_dRNA.out.fastq, MINIMAP2_dRNA.out.txt, MINIMAP2_dRNA.out.bam, MINIMAP2_dRNA.out.bai)
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

        else
        {
            RSEQC(FILTER_BAM.out.id, FILTER_BAM.out.bam, FILTER_BAM.out.bai, housekeeping)
        }
        
        
        BAMBU_PREP(FILTER_BAM.out.id, mapq, FILTER_BAM.out.bam, FILTER_BAM.out.bai, ref, annotation, MAKE_FAI.out, track_reads)

}
