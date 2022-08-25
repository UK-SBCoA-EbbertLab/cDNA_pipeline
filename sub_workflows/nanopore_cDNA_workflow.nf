// Import Modules
include {MAKE_FAI} from '../modules/make_fai'
include {MAKE_INDEX_cDNA} from '../modules/make_index'
include {CHM13_GTF; CHM13_GTF_ERCC} from '../modules/chm13_gff3_to_gtf'
include {PYCHOPPER} from '../modules/pychopper'
include {PYCOQC} from '../modules/pycoqc'
include {MINIMAP2_cDNA; MINIMAP2_QC} from '../modules/minimap2'
include {RSEQC} from '../modules/rseqc'
include {BAMBU_PREP; BAMBU_DISCOVERY} from '../modules/bambu'
include {STRINGTIE_ONT_cDNA} from '../modules/stringtie'
include {GFFCOMPARE} from '../modules/gffcompare'
include {MULTIQC_GRCh38 ; MULTIQC_CHM13} from '../modules/multiqc'


// TODO: MAKE QUALITY CONTROL WORK (RSEQC + PYCOQC + MULTIQC) ONCE WE HAVE FINAL DATA

workflow NANOPORE_cDNA {

    take:
        ref
        annotation
        housekeeping
        ont_reads_txt
        ont_reads_fq
        ercc
        cdna_kit

    main:
        MAKE_FAI(ref)
        MAKE_INDEX_cDNA(ref)
        PYCHOPPER(ont_reads_fq, ont_reads_txt, cdna_kit)
        MINIMAP2_cDNA(PYCHOPPER.out.id, PYCHOPPER.out.fastq,  MAKE_INDEX_cDNA.out)
        MINIMAP2_QC(PYCHOPPER.out.id, PYCHOPPER.out.original_fastq, PYCHOPPER.out.txt, MAKE_INDEX_cDNA.out)
        PYCOQC(MINIMAP2_QC.out.id, MINIMAP2_QC.out.txt, MINIMAP2_QC.out.bam, MINIMAP2_QC.out.bai)

        if (params.is_chm13 == true)
        {
            MULTIQC_CHM13(MINIMAP2_QC.out.multiQC.collect(), PYCOQC.out.multiQC.collect(), PYCHOPPER.out.multiQC.collect())

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
            RSEQC(MINIMAP2_cDNA.out.bam_mapped.collect(), MINIMAP2_cDNA.out.bai_mapped.collect(), housekeeping)
            MULTIQC_GRCh38(MINIMAP2_QC.out.multiQC.collect(), PYCOQC.out.multiQC.collect(), PYCHOPPER.out.multiQC.collect(), RSEQC.out.multiQC.collect())
        }
        
        BAMBU_PREP(MINIMAP2_cDNA.out.bam_mapped, MINIMAP2_cDNA.out.bai_mapped, ref, annotation, MAKE_FAI.out)
        BAMBU_DISCOVERY(BAMBU_PREP.out.collect(), ref, annotation, MAKE_FAI.out)
        STRINGTIE_ONT_cDNA(MINIMAP2_cDNA.out.id, MINIMAP2_cDNA.out.bam_mapped, MINIMAP2_cDNA.out.bai_mapped, BAMBU_DISCOVERY.out.gtf)
        GFFCOMPARE(BAMBU_DISCOVERY.out.gtf, annotation)
}
