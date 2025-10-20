#!/usr/bin/Rscript

library("bambu")

args <- commandArgs(trailingOnly = TRUE)

rc_files <- unlist(strsplit(args[1], ","))
fa_file <- args[2]
gtf_file <- args[3]
NDR_input <- args[4]
new_gene_and_isoform_prefix <- args[5]
directory <- args[6]

bambuAnnotations <- prepareAnnotations(gtf_file)


# Initialize a character vector to store the names of the newly created files
new_rc_files <- character(length(rc_files))

# Loop over each file name in the vector
for (i in seq_along(rc_files)) {

    rc_file <- rc_files[i]

    # Load the object from the RDS file
    rc <- readRDS(rc_file)

    # Remove the metadata elements
    metadata(rc)$readNames <- NULL
    metadata(rc)$readId <- NULL

    # Construct a new file name by appending "_notrackreads_version" before the .rds extension.
    new_rc_file <- sub("\\.rds$", "_notrackreads_version.rds", rc_file)

    # Save the modified object to the new file
    saveRDS(rc, file = new_rc_file)

    # Store the new file name in the new_files vector
    new_rc_files[i] <- new_rc_file
}

if (NDR_input == "auto") {
    extended_annotations <- bambu(reads=new_rc_files, annotations=bambuAnnotations, genome=fa_file,
                  lowMemory=TRUE, ncore=10, discovery=TRUE, quant=FALSE, trackReads=FALSE, verbose=TRUE, opt.discovery = list(prefix = new_gene_and_isoform_prefix))
} else {

    NDR_input <- as.double(NDR_input)

    extended_annotations <- bambu(reads=new_rc_files, annotations=bambuAnnotations, genome=fa_file,
        lowMemory=TRUE, ncore=10, NDR=NDR_input, discovery=TRUE, quant=FALSE, trackReads=FALSE, verbose=TRUE, opt.discovery=list(prefix = new_gene_and_isoform_prefix))
}


writeToGTF(extended_annotations, file=paste("./", directory, "/extended_annotations.gtf", sep=""))
saveRDS(extended_annotations, file=paste("./", directory, "/bambu_discovery.rds", sep=""))
