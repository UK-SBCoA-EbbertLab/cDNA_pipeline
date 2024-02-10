# Pipeline output structure and file description

This is a list of all possible outputs, it does not mean that your particular execution of the pipeline will have that output directory. For example, `mapping_cDNA` will only exists for 
cDNA when the parameter `--is_dRNA` is set to `false` and `mapping_dRNA` will only exist when the `--is_dRNA` parameter is set to `true`.


## Pipeline Output for STEP 1 (Basecalling):

    1. fast5_to_pod5 - One directory per sample that had any fast5 files converted into pod5 files for more efficient basecalling with Dorado.

    2. basecalling_output - Dorado basecalling output. One fastq file per sample and one sequencing summary file per sample. Reads for the
                            same run will  be separated into different fastq files based on barcode where appropriate.


## Pipeline Output for STEP 2 (QC, Alignment, and Bambu Pre-processing):
    
    1. bam_filtering - Samtols ".flagstat" reports and ".idxstat" reports after filtering the ".bam" alignment files by the MAPQ threshold
                        and only keeping primary alignments

    2. bambu_prep - Bambu pre-processed RDS file for each individual sample

    3. CHM13_gtf - Converted CHM13 annotation from GFF3 to GTF, will have ERCC annotations added to them if --ercc pipeline parameter
                    is set to "True"

    4. concatenated_fastq_and_sequencing_summary_files - Unzipped and concatenated fastq files for each file as well as matching sequencing summary files.
                                                            Only output if using the --path pipeline parameter

    5. contamination_report - One directory for each sample processed. Each sample has one ".bam" and ".bai" file of contaminant alignments.
                                Also a ".fastq" with reads that did not map to target alignment genome and a ".tsv" with the number of
                                reads aligned to each contaminant. Contaminants with more reads will be at the top of the file (sorted)

    6. dRNA_adapter_trimming_stats - One ".txt" file per sample containing the porechop output for the adapter trimming

    7. fai - Target genome fasta file index

    8. intermediate_qc_reports - Intermediate quality control reports for each sample separated into 4 directories: 
                                  "read_length", "number_of_reads", "quality_score_thresholds", "contamination".

    9. mapping_cDNA - Minimap2 alignment sorted ".bam" and ".bai" files including all alignments 
                        (no MAPQ filter, includes supplementary and secondary alignment)

    10. mapping_dRNA - Minimap2 alignment sorted ".bam" and ".bai" files including all alignments 
                        (no MAPQ filter, includes supplementary and secondary alignment)

    11. multiQC_input:
        a. RseQC - RseQC 3' bias gene body coverage ".txt" and ".r" files for each sample.
        b. minimap2 - Samtools ".flagstat" reports and ".idxstat" reports (no MAPQ filter, includes supplementary and secondary alignment).
        c. pychopper - Pychopper trimming statistics for each sample.
        d. pycoQC - Individual HTML and JSON files with quality control metrics for each sample.



## Pipeline Output for STEP 3:

    1. bambu_discovery - Output from bambu discovery step, includes gene level counts, transcript level counts,
                         transcript level unique counts, transcript level full-length counts, GTF extended annotations
                         containining known and new transcript coordinates, and RDS file with bambu object for all
                         files combined.

    2. bambu_quant - Output from bambu discovery step, includes gene level counts, transcript level counts,
                         transcript level unique counts, transcript level full-length counts, GTF annotations
                         containining known , and RDS file with bambu object for all files combined.

    3. contamination - Percent_Contaminant_Reads_mqc.tsv with the percent of contaminant reads for each sample.
    
    4. num_reads_report - Three reports, one with number of reads for each sample, other with reads length, and another with
                                MAPQ and PHRED quality scores used to filter the files.


    5. multiQC_output - MultiQC output files, most importantly the ".html" report showing summary statistics for all file.
