# Pipeline output structure and file description

This is a list of all possible outputs, it does not mean that your particular execution of the pipeline will have that output directory. For example, `mapping_cDNA` will only exists for 
cDNA when the parameter `--is_dRNA` is set to `false` and `mapping_dRNA` will only exist when the `--is_dRNA` parameter is set to `true`.


## Pipeline Output for STEP 1 (Basecalling):

    1. fast5_to_pod5

    2. basecalling_output


## Pipeline Output for STEP 2 (QC, Alignment, and Bambu Pre-processing):
    
    1. bam_filtering - <Samtols ".flagstat" reports and ".idxstat" reports after filtering the ".bam" alignment files by the MAPQ threshold
                        and only keeping primary alignments>

    2. bambu_prep - <Bambu pre-processed RDS file for each individual sample>

    3. CHM13_gtf - <Converted CHM13 annotation from GFF3 to GTF, will have ERCC annotations added to them if --ercc pipeline parameter
                    is set to "True">

    4. concatenated_fastq_and_sequencing_summary_files - <Unzipped and concatenated fastq files for each file as well as matching sequencing summary files.
                                                            Only output if using the --path pipeline parameter>

    5. contamination_report - <One directory for each sample processed. Each sample has one ".bam" and ".bai" file of contaminant alignments.
                                Also a ".fastq" with reads that did not map to target alignment genome and a ".tsv" with the number of
                                reads aligned to each contaminant. Contaminants with more reads will be at the top of the file (sorted)>

    6. dRNA_adapter_trimming_stats - <One ".txt" file per sample containing the porechop output for the adapter trimming>

    7. fai - <Target genome fasta file index>

    8. mapping_cDNA - <Minimap2 alignment sorted ".bam" and ".bai" files including all alignments 
                        (no MAPQ filter, includes supplementary and secondary alignment)>

    9. mapping_dRNA - <Minimap2 alignment sorted ".bam" and ".bai" files including all alignments 
                        (no MAPQ filter, includes supplementary and secondary alignment)>

    10. multiQC_input   



## Pipeline Output for STEP 3:

    1. bambu_discovery

    2. bambu_quant

    3. multiQC_output
