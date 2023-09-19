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

          --step                        <"1". Performs step 1>
          
          --fast5_dir                   <path to directory containing fast5 files. example: /sequencing_run/basecaling/fast5/">
          
          --basecall_config             <configuration name for basecalling setting. Example for PCS111 R9.4.1 PromethION flow cell: `dna_r9.4.1_450bps_hac_prom`>
          
          --basecall_id                 <Sample ID for the output fastq file that is being basecalled>



## Pipeline parameters for STEP 2 (Pre-processing)

          --step                        <"2". Performs step 2>

          --ont_reads_fastq             <path to fastq sequencing data, can submit multiple at once using pattern "/path/*.fastq".  If you don't specify this
                                        parameter the pipeline will not run. An exception is if you specify a  ".bam" file input, see Pipeline parameters for
                                        STEP 2 from BAM below. Default: "None">
          
          --ont_reads_txt               <path to sequencing summary files, can submit multiple at once. Make sure they follow same naming pattern as fastq
                                        files. If a fastq file is named "sample_1.fastq" then its sequencing summary file should be "sample_1.txt"
                                        or "sample_1_sequencing_summary.txt". If you don't specify this parameter pipeline will run, but skip PycoQC
                                        quality control step. Default: "None">
          
          --ref                         <path to reference/assembly ".fa" file. if using ERCC make sure to concatenate it to the end of the file.
                                        Default: "None">
  
          --annotation                  <path to reference annotation ".gtf" file for GRCh38 or ".gff3" for CHM13. If using GRCh38 and ERCC concatenate the
                                        two ".gtf" files prior to running the pipeline. If using ERCC with CHM13 make sure to set --is_chm13 to "True" and
                                        set --ercc to the path of your ERCC ".gtf" file. Default: "None">
  
          --out_dir                     <name of output directory for pipeline submission. Will appear under "results/<out_dir>" on the directory the pipeline
                                        was submitted from Default: "output_directory/">
  
          --ercc                        <path to ERCC annotation file. Only needed if using CHM13 reference + GFF3 annotation file and adding ERCC. Otherwise
                                        do not specify this parameter. Default: "None">
  
          --cdna_kit                    <option for pychopper trimming using adapters from the specific cDNA library version, current options are "PCS109",
                                        "PCS110", "PCS111". Default: "PCS111">
  
          --is_chm13                    <logical, set to "True" if using CHM13 and "False" if not. Fixes CHM13 annotation for compatibility with Bambu and
                                        converts to ".gtf". Default: "False">
          
          --housekeeping                <path to bed file containing coordinates for housekeeping genes used by RSeQC quality control tool. Currently only
                                        supported if --is_chm13 is set to "False". You can find the bed file with housekeeping gene coordinates for GRCh38
                                        here: https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38.HouseKeepingGenes.bed.gz/download
                                        Default: "None">
                               
           --track_reads                <Logical, set to "True" if you want Bambu to keep track of read assignments to transcripts in the output ".RDS" file
                                        from Bambu. Set to "False" if you don't need to keep track of read assignments (smaller files). Default: "False">
            
            --mapq                      <Integer, set it to the number you want to be used to filter ".bam" file by mapq. --mapq 10 filters out reads with MAPQ
                                        < 10. set it to 0 if don't want to filter out any reads. Default: 0>
                              
            --is_dRNA                   <Logical, set to "True if you want to run the minimap2 step with the parameters for "noisy" Nanopore directRNAseq reads.
                                         Set to "False" if you want to run it with minimap2 step for PCR amplified Nanopore cDNA reads. If you set this parameter
                                         to "True" you should omit the "--cdna_kit" parameter as that is only used for the Pychopper step of the pipeline and that
                                         step is skipped for directRNAseq. Note that this is the only step in the pipeline that is specific for 
                                         DirectRNAseq, Step 2 from BAM and Step 3 are not modified for dRNA vs cDNA. Step 1 will need the 
                                         specific basecalling configuration for dRNA instead of cDNA if you are running dRNA. Default: "False">

            --contamination_ref          <Path to a ".fasta" reference containing any contaminants you wish to align the unmapped reads against. For a 
                                          reference containing the contaminants we use in house got to the "contamination_reference_doc" folder in
                                          this GitHub or to this link on Zenodo: https://zenodo.org/deposit/8350277. If you pass a contamination
                                          reference a report will be generated with multiQC showing what percent of reads mapped to a contaminant.
                                          For further detail go into the "contamination_report/" folder under your pipeline execution results
                                          to find out which were the most abundant contaminants in each of your samples.>
  
 


##  Pipeline parameters for STEP 2 from BAM (Pre-processing)

          --step              <"2". Performs step 2>
          
          --bam               <path to ".bam" alignment files, using this will supersed any ".fastq" input>
          
          --bai               <path to ".bai" alignment index files, using this will supersed any ".fastq" input>
          
          --ref               <path to reference/assembly ".fa" file. if using ERCC make sure to concatenate it to the end of the file.
                              Default: "None">
  
          --annotation        <path to reference annotation ".gtf" file for GRCh38 or ".gff3" for CHM13. If using GRCh38 and ERCC concatenate the two
                              ".gtf" files prior to running the pipeline. If using ERCC with CHM13 make sure to set --is_chm13 to "True" and set
                              --ercc to the path of your ERCC ".gtf" file. Default: "None">
  
          --out_dir           <name of output directory for pipeline submission. Will appear under "results/<out_dir>" on the directory the pipeline
                              was submitted from. Default: "output_directory/">
  
          --ercc              <path to ERCC annotation file. Only needed if using CHM13 reference + GFF3 annotation file and adding ERCC.
                              Otherwise do not specify this parameter. Default: "None">
  
          --cdna_kit          <option for pychopper trimming using adapters from the specific cDNA library version, current options are "PCS109",
                              "PCS110", "PCS111". Default: "PCS111">
  
          --is_chm13          <logical, set to "True" if using CHM13 and "False" if not. Fixes CHM13 annotation for compatibility with Bambu and
                              converts to ".gtf". Default: "False">
          
           --track_reads      <logical, set to "True" if you want Bambu to keep track of read assignments to transcripts in the output ".RDS" file
                              from Bambu. Set to "False" if you don't need to keep track of read assignments (smaller files). Default: "False">
            
            --mapq            <integer, set it to the number you want to be used to filter ".bam" file by mapq. --mapq 10 filters out reads with
                              MAPQ < 10. Set it to 0 if don't want to filter out any reads. Default: 0>

## Pipeline parameters for STEP 3
  
          --bambu_rds         <path to individually pre-processed bambu RDS objects (output from step 1). Default: "None"
          
          --ref               <path to reference/assembly ".fa" file. if using ERCC make sure to concatenate it to the end of the file.
                              Default:"None">
          
          --fai               <path to reference index ".fai" file>
  
          --annotation        <path to reference annotation ".gtf" file for GRCh38 or ".gff3" for CHM13. If using GRCh38 and ERCC concatenate the two
                              ".gtf" files prior to running the pipeline. If using ERCC with CHM13 make sure to set --is_chm13 to "True" and set
                              --ercc to the path of your ERCC ".gtf" file. Default: "None">
  
          --out_dir           <name of output directory for pipeline submission. Will appear under "results/<out_dir>" on the directory the pipeline
                              was submitted from. Default: "output_directory/">
          
          --is_discovery      <Logical, if "True" perform transcript discovery and quantification with Bambu, else if "False" perform only
                              quantification based on the given GTF annotation>
  
          --NDR               <NDR parameter for performing transcript discovery with Bambu. Values can range from 0-1, the closer to 0 the more
                              stringent the transcript discovery is. Meaning you will discover less transcripts but be more certain that they are
                              real. "auto" option lets Bambu determine which NDR threshold to use based on the data. Authors of the tool recommend
                              "auto" for most situations. Default: "auto">
  
          
          --track_reads       <logical, set to "True" if you want Bambu to keep track of read assignments to transcripts in the output ".RDS" file
                              from Bambu. Set to "False" if you don't need to keep track of read assignments (smaller files). Default: "False">
            
          --multiqc_input    <path to directory containing multiqc input data from pipeline step 2. Use path and add **.
                              Example: /path/multiqc_input/** Default: "None">
          
          --multiqc_config   <path to multiqc ".yaml" config file. Default: "None". PS: You can find an example of a multiqc config file that works
                              for this pipeline under /cDNA_pipeline/workflow/bin/multiqc_config.yaml>

## Examples of the submissions

### Example for Step 1: Basecalling 
       nextflow ../main.nf --step 1 \
              --fast5_dir "../ont_data/test_data/fast5_directory/" \
              --basecall_config "dna_r9.4.1_450bps_hac_prom" \
              --basecall_id "test_sample_1" 

### Example for step 2: CHM13 without ERCCs

          nextflow ../main.nf --step 2 \
              --ont_reads_fq "../ont_data/test_data/*.fastq" \
              --ont_reads_txt "../ont_data/test_data/*.txt" \
              --ref "../../references/chm13v2.0.fa" \
              --annotation "../../references/CHM13.v2.0.gff3" \
              --ercc "None" \
              --out_dir "./CHM13_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "True" \
              --track_reads "False" \
              --mapq "0" 

### Example for step 2: CHM13 with ERCCs

          nextflow ../main.nf --step 2 \
              --ont_reads_fq "../ont_data/test_data/*.fastq" \
              --ont_reads_txt "../ont_data/test_data/*.txt" \
              --ref "../../references/chm13v2.0_ERCC.fa" \
              --annotation "../../references/CHM13.v2.0.gff3" \
              --ercc "../../references/ERCC92.gtf" \
              --out_dir "./CHM13_ERCC_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "True" \
              --track_reads "False" \
              --mapq "0" 
    
#### Notice that for CHM13 you need to concatenate CHM13 and the ERCC reference prior to submitting the pipeline, but the annotations are entered separately and concatenated by the program itself after converting CHM13 annotation to ".gtf" format.

### Example for step 2 (cDNA): GRCh38 without ERCCs

          nextflow ../main.nf --step 2 \ 
              --ont_reads_fq "../ont_data/test_data/*.fastq" \
              --ont_reads_txt "../ont_data/test_data/*.txt" \
              --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
              --annotation "../../references/Homo_sapiens.GRCh38.106.gtf" \
              --out_dir "./GRCh38_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "False" \
              --track_reads "False" \
              --mapq "0" \
              --housekeeping "../../references/hg38.HouseKeepingGenes.bed"


### Example for step 2 (cDNA): GRCh38 with ERCCs

          nextflow ../main.nf --step 2 \ 
              --ont_reads_fq "../ont_data/test_data/*.fastq" \
              --ont_reads_txt "../ont_data/test_data/*.txt" \
              --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
              --annotation "../../references/Homo_sapiens.GRCh38.106_ERCC.gtf" \
              --out_dir "./GRCh38_ERCC_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "False" \
              --track_reads "False" \
              --mapq "0" \
              --housekeeping "../../references/hg38.HouseKeepingGenes.bed"
              
### Example for step 2 (dRNA): GRCh38 with ERCCs

          nextflow ../main.nf --ont_reads_fq "/scratch/bag222/data/ont_data/2023-06-06_brain_directRNA_intronic_reads/*.fastq" \
                    --ont_reads_txt "/scratch/bag222/data/ont_data/2023-06-06_brain_directRNA_intronic_reads/*.txt" \
                    --ref "../../../../cDNA_pipeline/references/Homo_sapiens.GRCh38_ERCC.fa" \
                    --annotation "../../../../cDNA_pipeline/references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
                    --housekeeping "../../../../cDNA_pipeline/references/hg38.HouseKeepingGenes.bed" \
                    --out_dir "./GRCh38_dRNA_EERCC_test/" \
                    --is_discovery "True" \
                    --bambu_track_reads "True" \
                    --is_dRNA "True" \
                    --mapq "0" \
                    --step "2" \
                    --is_chm13 "False"
 
 
#### Notice that for GRCh38 the `--ercc` is not needed as the user can easily concatenate both the GRCh38 reference and the annotation to the ERCC reference and annotation prior to running the analysis.     


### Example for step 2 from BAM: GRCh38 with ERCCs

          nextflow ../main.nf --step 2 \ 
              --bam "./results/GRCh38_ERCC_test/mapping_cDNA/*.bam" \
              --bai "./results/GRCh38_ERCC_test/mapping_cDNA/*.bai" \
              --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
              --annotation "../../references/Homo_sapiens.GRCh38.106_ERCC.gtf" \
              --out_dir "./GRCh38_ERCC_test_mapq10_track_reads/" \
              --cdna_kit "PCS111" \
              --is_chm13 "False" \
              --track_reads "True" \
              --mapq "10" 
              
### Example for step 3: GRCh38 with ERCCs

          nextflow ../main.nf --step 3 \
              --bambu_rds "./results/GRCh38_ERCC_test/bambu_prep/*.rds" \
              --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
              --fai "./results/GRCh38_ERCC_test/fai/*.fai \
              --annotation "../../references/Homo_sapiens.GRCh38.106_ERCC.gtf" \
              --is_discovery "True" \
              --track_reads "False" \
              --NDR "auto" \
              --multiqc_input "./results/GRCh38_ERCC_test/multiQC_input/**" \
              --multiqc_config "../../references/multiqc_config.yaml" \
              --out_dir "./GRCh38_ERCC_test/" \
              --is_chm13 "False"       
