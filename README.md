# cDNA_pipeline

## Getting Started

### 1) Have a functioning version of Nextflow in your Path.

- Information on how to install NextFlow can be found [here](https://www.nextflow.io/docs/latest/getstarted.html).
          
### 2) Have a functioning version of Singularity on your Path.

- Information on how to install Singularity cna be found [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
          
          
### 3) Clone this github repo using the command below

          git clone https://github.com/UK-SBCoA-EbbertLab/cDNA_pipeline


### 4) Go into the `./workflows/nextflow.config` file and make any necessary changes:

- Alter slurm job manager (or other job manager) parameters to suit your local environment. I don't recommend changing the memory/cpu/time allocated 
for the job manager.
        

          
### 5) Make sure you have all the sequencing files and reference genomes/assemblies files and annotation files you will need to run the pipeline.
          
- ".fastq" -- Nanopore cDNA sequencing files or ".bam" alignment files.

- "sequencing_summary.txt" -- These files are not necessary for execution, but if not available the PycoQC quality control step will be skipped.

- refecence/assembly ".fa" file.

- annotation ".gtf" file is preffered. Only use ".gff3" if using CHM13. Pipeline has an option to handle this, see `Pipeline parameters for STEP 2`.
          


## Pipeline parameters for STEP 1 (Basecalling)

          --fast5_dir   <path to directory containing fast5 files. example: /sequencing_run/basecaling/fast5/">
          
          --basecall_config     <configuration name for basecalling setting. Example for PCS111 R9.4.1 PromethION flow cell: `dna_r9.4.1_450bps_hac_prom`>
          
          --basecall_id         <Sample ID for the output fastq file that is being basecalled>



## Pipeline parameters for STEP 2

          --ont_reads_fastq   <path to fastq sequencing data, can submit multiple at once using pattern "/path/*.fastq". If you don't specify this parameter
                               the pipeline will not run. An exception is if you specify a ".bam" file input, see Pipeline parameters for STEP 2 from BAM below.
                               Default: "None">
          
          --ont_reads_txt     <path to sequencing summary files, can submit multiple at once. Make sure they follow same naming pattern as fastq files.
                               if a fastq file is named "sample_1.fastq" then its sequencing summary file should be "sample_1.txt" or "sample_1_sequencing_summary.txt".
                               If you don't specify this parameter pipeline will run, but skip PycoQC quality control step. Default: "None">
          
          --ref               <path to reference/assembly ".fa" file. if using ERCC make sure to concatenate it to the end of the file. Default: "None">
  
          --annotation        <path to reference annotation ".gtf" file for GRCh38 or ".gff3" for CHM13. If using GRCh38 and ERCC concatenate the two ".gtf" files 
                               prior to running the pipeline. If using ERCC with CHM13 make sure to set --is_chm13 to "True" and set --ercc to the path of your
                               ERCC ".gtf" file. Default: "None">
  
          --out_dir           <name of output directory for pipeline submission. Will appear under "results/<out_dir>" on the directory the pipeline was submitted from
                               Default: "output_directory/">
  
          --ercc              <path to ERCC annotation file. Only needed if using CHM13 reference + GFF3 annotation file and adding ERCC.
                               Otherwise do not specify this parameter. Default: "None">
  
          --cdna_kit          <option for pychopper trimming using adapters from the specific cDNA library version, current options are "PCS109", "PCS110", "PCS111".
                               Default: "PCS111">
  
          --is_chm13          <logical, set to "True" if using CHM13 and "False" if not. Fixes CHM13 annotation for compatibility with Bambu and converts to ".gtf".
                               Default: "False">
          
          --housekeeping      <path to bed file containing coordinates for housekeeping genes used by RSeQC quality control tool. Currently only supported if 
                               --is_chm13 is set to "False". You can find the bed file with housekeeping gene coordinates for GRCh38 here: 
                               https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38.HouseKeepingGenes.bed.gz/download
                               Default: "None">
                               
           --track_reads      <Logical, set to "True" if you want Bambu to keep track of read assignments to transcripts in the output ".RDS" file from Bambu.
                               Set to "False" if you don't need to keep track of read assignments (smaller files). Default: "False">
            
            --mapq            <Integer, set it to the number you want to be used to filter ".bam" file by mapq. --mapq 10 filters out reads with MAPQ < 10. 
                               set it to 0 if don't want to filter out any reads. Default: 0>
  
  
## Pipeline parameters for STEP 3
  


## Examples of the submissions


### CHM13 without ERCCs

          nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
              --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
              --ref "../../references/chm13v2.0.fa" \
              --annotation "../../references/CHM13.v2.0.gff3" \
              --ercc "None" \
              --out_dir "./CHM13_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "True"  -resume

### CHM13 with ERCCs

          nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
              --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
              --ref "../../references/chm13v2.0_ERCC.fa" \
              --annotation "../../references/CHM13.v2.0.gff3" \
              --ercc "../../references/ERCC92.gtf" \
              --out_dir "./CHM13_ERCC_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "True"  -resume
    
#### Notice that for CHM13 you need to concatenate CHM13 and the ERCC reference prior to submitting the pipeline, but the annotations are entered separately and concatenated by the program itself after converting CHM13 annotation to ".gtf" format.

### GRCh38 without ERCCs

          nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
              --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
              --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
              --annotation "../../references/Homo_sapiens.GRCh38.106.gtf" \
              --out_dir "./GRCh38_test/" \
              --ercc "None" \
              --cdna_kit "PCS111" \
              --is_chm13 "False"  -resume


### GRCh38 with ERCCs

          nextflow ../main.nf --ont_reads_fq "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.fastq" \
              --ont_reads_txt "/mnt/gpfs3_amd/condo/mteb223/bag222/data/cdna_comparison_project/2019_ont_data/test_data/*.txt" \
              --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
              --annotation "../../references/Homo_sapiens.GRCh38.106_ERCC.gtf" \
              --out_dir "./GRCh38_ERCC_test/" \
              --ercc "None" \
              --cdna_kit "PCS111" \
              --is_chm13 "False"  -resume
    
#### Notice that for GRCh38 the `--ercc` is always set to "None" as there the user can easily concatenate both the GRCh38 reference and the annotation to the ERCC reference and annotation.          
