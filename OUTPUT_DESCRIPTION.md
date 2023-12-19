# Pipeline output structure and file description

This is a list of all possible outputs, it does not mean that your particular execution of the pipeline will have that output directory. For example, `mapping_cDNA` will only exists for 
cDNA when the parameter `--is_dRNA` is set to false and `mapping_dRNA` will only exist when the `--is_dRNA` parameter is set to true.


## Pipeline Output for STEP 1 (Basecalling):

    1. fast5_to_pod5

    2. basecalling_output


## Pipeline Output for STEP 2 (QC, Alignment, and Bambu Pre-processing):
    
    1. bam_filtering

    2. bambu_prep

    3. concatenated_fastq_and_sequencing_summary_files

    4. contamination_report

    5. dRNA_adapter_trimming_stats

    6. fai

    7. mapping_cDNA

    8. mapping_dRNA

    9. multiQC_input   



## Pipeline Output for STEP 2 from BAM (Pre-processing):


    1. bam_filtering

    2. bambu_prep


## Pipeline Output for STEP 3:

    1. bambu_discovery

    2. bambu_quant

    3. multiQC_output
