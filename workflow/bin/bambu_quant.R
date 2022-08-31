#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

rc_files <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]

bambuAnnotations <- prepareAnnotations(gtf_file)

se_novel <- bambu(rcFile=rc_files, annotations=bambuAnnotations, genome=fa_file,
                  ncore=8, lowMemory=TRUE, discovery=FALSE, quant=TRUE)

writeBambuOutput(se_novel, path = "./bambu_quant/")
