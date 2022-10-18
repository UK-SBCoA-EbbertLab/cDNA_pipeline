#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

rc_files <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]
NDR_input <- as.double(args[4])

bambuAnnotations <- prepareAnnotations(gtf_file)

se_novel <- bambu(rcFile=rc_files, annotations=bambuAnnotations, genome=fa_file,
                  lowMemory=TRUE, NDR=NDR_input, discovery=TRUE, quant=TRUE, trackReads=TRUE)

writeBambuOutput(se_novel, path = "./bambu_discovery/")

saveRDS(se_novel, file="./bambu_discovery/final_discovery.RDS")
