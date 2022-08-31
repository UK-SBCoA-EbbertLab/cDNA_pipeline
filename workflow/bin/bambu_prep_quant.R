#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

bam <- args[1]
fa_file <- args[2]
gtf_file <- args[3]

bambuAnnotations <- prepareAnnotations(gtf_file)

se <- bambu(reads = bam, annotations = bambuAnnotations, genome = fa_file, rcOutDir = "./bambu_prep_quant/",
                  ncore=8, lowMemory=TRUE, quant=TRUE, discovery=FALSE, yieldSize = 1000000)

