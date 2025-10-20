#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

rc_file <- args[1]
fa_file <- args[2]
gtf_file <- args[3]
track_reads_input <- args[4] == "true"

name_without_extension <- sub("\\.rds$", "", rc_file)
output_dir = paste("./", name_without_extension, "/", sep="")

bambuAnnotations <- prepareAnnotations(gtf_file)

se_quant <- bambu(reads=rc_file, annotations=bambuAnnotations, genome=fa_file,
                  ncore=12, lowMemory=FALSE, discovery=FALSE, quant=TRUE, verbose=TRUE, trackReads=track_reads_input)

writeBambuOutput(se_quant, path=output_dir)

saveRDS(se_quant, file=paste("./", name_without_extension, "_quant.rds", sep=""))


file_names <- c(paste(output_dir, "CPM_transcript.txt", sep=""), 
                paste(output_dir, "counts_gene.txt", sep=""), 
                paste(output_dir, "counts_transcript.txt", sep=""),
                paste(output_dir, "fullLengthCounts_transcript.txt", sep=""), 
                paste(output_dir, "uniqueCounts_transcript.txt", sep=""))


for (old_path in file_names) {

    base_name <- basename(old_path)
    
    new_path = paste(paste("./", name_without_extension, sep=""), base_name, sep="_")

    file.rename(old_path, new_path)

}

