#!/usr/bin/env Rscript

# Load necessary libraries
library(edgeR)
library(DESeq2)
library(tools)

# Define a function for CPM normalization
normalize_cpm <- function(counts) {
	cpm_counts <- cpm(counts, log = FALSE)
	return(cpm_counts)
}

# Define a function for TMM normalization
normalize_tmm <- function(counts) {
	dge <- DGEList(counts = counts)
	dge <- calcNormFactors(dge, method = "TMM")
	tmm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
	return(tmm_counts)
}

# Define a function for Median of Ratios normalization
normalize_deseq2 <- function(counts) {
	counts <- round(counts)
	dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(row.names = colnames(counts)), design = ~1)
	dds <- estimateSizeFactors(dds)
	deseq2_counts <- counts(dds, normalized=TRUE)
	return(deseq2_counts)
}

# Define a function to calculate log(expression + 1)
log_transform <- function(counts) {
	log_counts <- log2(counts + 1)
	return(log_counts)
}

convert_counts_matrix <- function(counts_matrix) {

	# Check if the rownames contain the string "_"
	if (all(grepl("_", rownames(counts_matrix)))) {
		# Split the rownames into GENEID and TRANSCRIPTID
		gene_transcript <- do.call(rbind, strsplit(rownames(counts_matrix), "_"))
		
		# Create a new data frame with GENEID and TRANSCRIPTID as the first two columns
		new_counts_matrix <- cbind(data.frame(transcript_id = gene_transcript[, 2], gene_id = gene_transcript[, 1]),
					   counts_matrix)

		# Remove rownames
		rownames(new_counts_matrix) <- NULL
	} else {
		# Extract GENEID
		geneid <- rownames(counts_matrix)

		# Create a new data frame with GENEID as the first column
		new_counts_matrix <- cbind(data.frame(gene_id = geneid), counts_matrix)

    		# Remove rownames
		rownames(new_counts_matrix) <- NULL
	}

	return(new_counts_matrix)
}

# Create a function to save the counts matrix to a new directory
save_counts_matrix <- function(counts_matrix, file_path) {
	# Extract the directory from the file path
	dir_path <- dirname(file_path)
  
	# Create the directory if it does not exist
	if (!dir.exists(dir_path)) {
		dir.create(dir_path, recursive = TRUE)
	}
  
	# Write the counts matrix to the specified file path
	write.table(counts_matrix, file_path, sep = "\t", row.names = FALSE, quote = FALSE)
}




## -------------------- main function --------------------------- ##


args = commandArgs(trailingOnly = TRUE)

directory = args[1]
output_directory = paste(args[2], "/normalization/", sep="")

files <- list.files(path = directory, pattern = "\\.tsv$", full.names=TRUE)

is_transcript=FALSE

for (file in files) {
	print(paste("PROCCESSING FILE: ", file, sep=""))

	basename = tools::file_path_sans_ext(basename(file))

	counts <- read.table(file, header=TRUE, sep="\t")

	if ("TXNAME" %in% colnames(counts) {
		colnames(counts)[colnames(counts) == "TXNAME"] <- "transcript_id"
	}
        if ("GENEID" %in% colnames(counts) {
		colnames(counts)[colnames(counts) == "GENEID"] <- "gene_id"
	}	    
	if ("transcript_id" %in% colnames(counts)) {
		rownames(counts) <- paste(counts$gene_id, counts$transcript_id, sep="_")
		counts <- counts[ , !(names(counts) %in% c("gene_id", "transcript_id"))]
		print("IS TRANSCRIPT")
	} else {
		rownames(counts) <- counts$gene_id
		counts <- counts[ , !(names(counts) %in% c("gene_id"))]
	}

	# Apply normalization methods
	cpm <- normalize_cpm(counts)
	tmm <- normalize_tmm(counts)
	deseq2 <- normalize_deseq2(counts)

	# Apply log transformation
	log_counts <- log_transform(counts)
	log_cpm <- log_transform(cpm)
	log_tmm <- log_transform(tmm)
	log_deseq2 <- log_transform(deseq2)

	# Fix column names and rownames before outputting
	counts <- convert_counts_matrix(counts)
	cpm <- convert_counts_matrix(cpm)
	tmm <- convert_counts_matrix(tmm)
	deseq2 <- convert_counts_matrix(deseq2)
	log_counts <- convert_counts_matrix(log_counts)
	log_cpm <- convert_counts_matrix(log_cpm)
	log_tmm <- convert_counts_matrix(log_tmm)
	log_deseq2 <- convert_counts_matrix(log_deseq2)

	 # Write output to TSV files
	save_counts_matrix(counts, file=paste(output_directory, basename, "_counts.tsv", sep=""))
	save_counts_matrix(cpm, file=paste(output_directory, basename, "_cpm.tsv", sep=""))
	save_counts_matrix(tmm, file=paste(output_directory, basename, "_tmm.tsv", sep=""))
	save_counts_matrix(deseq2, file=paste(output_directory, basename, "_median_of_ratios.tsv", sep=""))

	save_counts_matrix(log_counts, file=paste(output_directory, basename, "_counts_log2.tsv", sep=""))
	save_counts_matrix(log_cpm, file=paste(output_directory, basename, "_cpm_log2.tsv", sep=""))
	save_counts_matrix(log_tmm, file=paste(output_directory, basename, "_tmm_log2.tsv", sep=""))
	save_counts_matrix(log_deseq2, file=paste(output_directory, basename, "_median_of_ratios_log2.tsv", sep=""))

}
