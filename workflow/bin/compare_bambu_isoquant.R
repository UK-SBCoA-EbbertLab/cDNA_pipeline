#!/usr/bin/Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

bambu_path <- args[1]
isoquant_path <- args[2]
gffcompare_path <- args[3]

bambu_counts <- read_tsv(bambu_path)
isoquant_counts <-read_tsv(isoquant_path)


# Compare novel isoforms ----------------------------
# columns at end -> bambu_gene_id	bambu_transcript_id	X4	isoquant_gene_id	isoquant_transcript_id
gffcompare <- read_tsv(gffcompare_path, col_names = FALSE) %>%
	select(c(X3,X4,X5)) %>%
	filter(if_any(.cols = c(X3,X5), .fns = ~ str_detect(., "Bambu|transcript"))) %>%
	separate(X3, c("bambu_gene_id", "bambu_transcript_id"), sep="\\|") %>%
	separate(X5, c("isoquant_gene_id", "isoquant_transcript_id"), sep="\\|", extra="drop") %>%
	mutate(isoquant_gene_id = str_remove(isoquant_gene_id, "q1:"))


# columns at end -> bambu_gene_id	bambu_transcript_id	X4	isoquant_gene_id	isoquant_transcript_id	shared_id
novel_both <- gffcompare %>%
	filter((grepl("Bambu", bambu_transcript_id)) & (grepl("transcript", isoquant_transcript_id)) & (X4 == "=")) %>%
	mutate(shared_id = paste0("NewShared_", row_number()))

# columns at end -> bambu_gene_id	bambu_transcript_id	X4	isoquant_gene_id	isoquant_transcript_id
novel_bambu <- gffcompare %>%
	filter((grepl("Bambu", bambu_transcript_id)) & (X4 != "=")) %>%

# columns at end -> bambu_gene_id	bambu_transcript_id	X4	isoquant_gene_id	isoquant_transcript_id
novel_isoquant <- gffcompare %>%
	filter((grepl("transcript", isoquant_transcript_id)) & (X4 != "="))

# ----------------------------------------------------

# Grab list of new isoform IDs -----------------------

bambu_ids <- bambu_counts %>% filter(grepl("Bambu", TXNAME)) %>% pull(TXNAME)
isoquant_ids <- isoquant_counts %>% filter(grepl("transcript", transcript_id)) %>% pull(transcript_id)

# ----------------------------------------------------

# Convert to long format -----------------------------
# columns at end -> transcript_id	gene_id	sample	bambu_count
bambu_counts_long <- bambu_counts %>%
	pivot_longer(-c("transcript_id", "gene_id"), names_to = "sample", values_to = "bambu_count")

# columns at end -> transcript_id	gene_id	sample	isoquant_count
isoquant_counts_long <- isoquant_counts %>%
	pivot_longer(-c("transcript_id", "gene_id"), names_to = "sample", values_to = "isoquant_count")

# ----------------------------------------------------

# Compare total number of known isoforms expressed ---------
# columns at end -> sample	bambu
bambu_isoforms_by_sample <- bambu_counts_long %>%
	group_by(sample) %>%
	summarise(bambu = n())

# columns at end -> sample	isoquant
isoquant_isoforms_by_sample <- isoquant_counts_long %>%
	group_by(sample) %>%
	summarise(isoquant = n())

# columns_at_end -> sample	tool	n_isoforms
isoforms_by_sample <- inner_join(bambu_isoforms_by_sample,
				isoquant_isoforms_by_sample,
				by = c("sample")) %>%
	pivot_longer(c("bambu", "isoquant"), names_to = "tool", values_to = "n_isoforms")

ggplot(isoforms_by_sample, aes(sample, n_isoforms, fill = tool) +
	geom_col(position="dodge")

ggsave("known_isoforms_expressed_by_sample.pdf")

# ----------------------------------------------------

# Compare number of isoforms expressed per gene ------
# columns at end -> gene_id	sample	n_isoform_bambu
bambu_isoforms_per_gene <- bambu_counts_long %>%
	group_by(gene_id, sample) %>%
	summarise(n_isoform_bambu = n())

# columns at end -> gene_id	sample	n_isoform_isoquant
isoquant_isoforms_per_gene <- isoquant_counts_long %>%
	group_by(gene_id, sample) %>%
	summarise(n_isoform_isoquant = n())

# columns at end -> gene_id	sample	n_isoform_bambu	n_isoform_isoquant
isoforms_per_gene <- full_join(bambu_isoforms_per_gene,
			       isoquant_isoforms_per_gene,
			       by = c("gene_id", "sample"))

ggplot(isoforms_per_gene, aes(n_isoform_bambu, n_isoform_isoquant, fill=sample)) +
	geom_jitter(height = 0.05, width = 0.05, seed = 12, alpha = 0.3)

ggsave("n_isoforms_per_gene_by_sample.pdf")

# ----------------------------------------------------

# Compare isoform expression -------------------------
# columns at end ->	transcript_id	gene_id	sample	bambu_count	isoquant_count
isoform_expression_comparison <- inner_join(bambu_counts_long,
					   isoquant_counts_long,
					   by = c("transcript_id", "gene_id", "sample"))

ggplot(isoform_expression_comparison, aes(bambu_count, isoquant_count, fill=sample)) +
	geom_point(alpha = 0.3)

ggsave("expression_of_known_isoforms.pdf")

# ----------------------------------------------------


# Compare shared new isoforms expressed --------------
# columns at end -> transcript_id	gene_id	sample	bambu_count	shared_id
shared_bambu <- bambu_counts_long %>%
	inner_join(novel_both %>%
		rename(transcript_id = bambu_transcript_id) %>%
		select(transcript_id, shared_id), by = c("transcript_id")) 

# columns at end -> transcript_id	gene_id	sample	isoquant_count	shared_id
shared_isoquant <- isoquant_counts_long %>%
	inner_join(novel_both %>%
		rename(transcript_id = isoquant_transcript_id) %>%
		select(transcript_id, shared_id), by = c("transcript_id")) 

# columns at end -> gene_id	sample	shared_id	bambu_count	isoquant_count
shared_counts <- shared_bambu %>% 
	select(-transcript_id) %>%
	inner_join(shared_isoquant %>% select(-transcript_id), by = c("shared_id", "gene_id", "sample"))

ggplot(shared_counts, aes(bambu_counts, isoquant_counts, fill = sample)) +
	geom_point(alpha=0.3)

ggsave("expression_of_new_shared_isoforms.pdf")

# ----------------------------------------------------

# Compare expression of new isoforms (overlapping and not)
# columns at end -> transcript_id       gene_id sample  counts	tool
new_bambu <- bambu_counts_long %>%
	filter(transcript_id %in% bambu_ids) %>%
	rename(counts = bambu_counts) %>%
	mutate(tool = "bambu")

# columns at end -> transcript_id       gene_id sample  counts	tool
new_isoquant <- isoquant_counts_long %>%
	filter(transcript_id %in% isoquant_ids) %>%
	rename(counts = isoquant_counts) %>%
	mutate(tool = "isoquant")

new_combined <- bind_rows(new_bambu, new_isoquant)

ggplot(new_combined, aes(x=counts, fill=tool, color=tool)) +
	

