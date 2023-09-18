## This is a documentation on how we built the master_contamination_reference.fa file to use for checking if the unmapped reads 
## in the RNAseq data are coming from Bacteria, Fungi, Viruses, or Nanopore PCS/SQK-RNA kit primers/adapters.


### 1. We took the blast database files from Decontaminer for BACTERIA, FUNGI, and VIRUSES: https://drive.google.com/drive/u/2/folders/1JSpAi6sghxyCrB5XqDOInF-P12WyxcE-
### Link to Decontaminer paper, including explanation on how databases were generated: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2684-x

### 2. We converted each database into fasta files using the command: `blastdbcmd -entry all -db <database label> -out <outfile>`

### 3. We concatenated the 3 fasta files generated from each blast database using the command: `cat virus.fa bacteria.fa fungi.fa >> master_contamination_reference.fa`

### 4. We added nanopore RNA/cDNA kit adapters and primers to the master_contamination_reference.fa file.
### These sequences can be found at: https://community.nanoporetech.com/technical_documents/chemistry-technical-document/v/chtd_500_v1_revaj_07jul2016/rna-sequencing-direct-and-via-cdna

### 5. We checked for any duplicate IDs on the master_contamination_reference.fa file using the command: `cat master_contamination_reference.fa | grep '^>' | sort | uniq -d`

### 6. Duplicate entries were manually removed from the master_contamination_reference.fa, so that there would be only one entry for each sequence ID.

### 7. Spaces were substituted with underscores on sequence names using the following command: `tr ' ' _ < master_contaminant_reference.fasta > test.fasta && mv test.fasta master_contaminant_reference.fasta`

### 8. The final decontamination reference was deposited on Zenodo, under the following link: https://zenodo.org/record/8350277
### DOI: 10.5281/zenodo.8350277

### All data were retrieved on September 15th 2023,
