#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

bam <- args[1]
fa_file <- args[2]
gtf_file <- args[3]
track_reads_input <- args[4] == "true"

bambuAnnotations <- prepareAnnotations(gtf_file)

se <- bambu(reads = bam, annotations = bambuAnnotations, genome = fa_file, rcOutDir = "./bambu_prep/", trackReads=track_reads_input,
                        lowMemory=FALSE, quant=FALSE, discovery=FALSE, yieldSize = 5000000, verbose=TRUE)
