// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2


log.info """
 RNA-SEQ ONT CORRECTION AND DISCOVERY PIPELINE
 ===============================================
 nanopore fastq files                   : ${params.ont_reads_fq}
 nanopore sequencing summary files      : ${params.ont_reads_txt}
 reference genome                       : ${params.ref}
 reference annotation                   : ${params.annotation}
 """


// Import Workflows
include {NANOPORE_cDNA} from '../sub_workflows/nanopore_cDNA_workflow'
include {BASECALLER} from '../sub_workflows/basecall_workflow'

// Define initial files and channels
ont_reads_fq = Channel.fromPath(params.ont_reads_fq).map { file -> tuple(file.baseName, file) }
ont_reads_txt = Channel.fromPath(file(params.ont_reads_txt))
ref = file(params.ref)
annotation = file(params.annotation)
fast5_dir = Channel.fromPath(params.fast5_dir)
basecall_config = Channel.from(params.basecall_config)
basecall_id = Channel.from(params.basecall_id)
cdna_kit = Channel.value(params.cdna_kit)

if (params.ercc != "None") {
    ercc = Channel.fromPath(params.ercc)
    }
else {
    ercc = Channel.from(params.ercc)
    }

// Make sure ONT sequencing summary and fastq files are in the same order
ont_reads_txt = ont_reads_txt.toSortedList( { a, b -> a.baseName <=> b.baseName } ).flatten()
ont_reads_fq = ont_reads_fq.toSortedList( { a, b -> a[0] <=> b[0] } ).flatten().buffer(size:2)

workflow {


    if (params.ont_reads_fq != "None"){
        NANOPORE_cDNA(ref, annotation, ont_reads_txt, ont_reads_fq, ercc, cdna_kit)
    }

    if (params.fast5_dir != "None"){
        BASECALLER(fast5_dir, basecall_config, basecall_id)
    }
}
