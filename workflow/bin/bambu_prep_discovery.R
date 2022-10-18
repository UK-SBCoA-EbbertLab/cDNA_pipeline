#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

bam <- args[1]
fa_file <- args[2]
gtf_file <- args[3]
NDR_input <- as.double(args[4])

bambuAnnotations <- prepareAnnotations(gtf_file)

se <- bambu(reads = bam, annotations = bambuAnnotations, genome = fa_file, rcOutDir = "./bambu_prep_discovery/", trackReads=TRUE,
                        NDR=1, ncore=8, lowMemory=TRUE, quant=FALSE, discovery=TRUE, yieldSize = 1000000 )

