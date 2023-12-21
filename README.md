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
          
- ".fastq" -- Nanopore cDNA sequencing files or ".bam" alignment files. Could also be the small ".fastq.gz" files for each sample, see
              the `--path` option below if you are processing files as they are outputted from basecalling in an Oxford Nanopore Sequencer.
              Could also be ".pod5" or ".fast5" files if you are starting from the basecalling step. 

- "sequencing_summary.txt" -- These files are not necessary for execution, but if not available the PycoQC quality control step will be skipped.

- refecence/assembly ".fa" file.

- annotation ".gtf" file is preffered. Only use ".gff3" if using CHM13. Pipeline has an option to handle this, see `Pipeline parameters for STEP 2`.


##  
## For documentation of pipeline output click on the `OUTPUT_DESCRIPTION.md` file in this repository.
## 


## Pipeline parameters for STEP 1 (Basecalling)

Many of the parameters for this step are based on dorado basecaller, see their [documentation](https://github.com/nanoporetech/dorado) to understand it better.

PS: I took care to pass the base modification information to the FASTQ headers so that it is still contained even in the fastq files. I prefer outputting
fastq files because it allows more flexibility in downstream alignment and processing.

          --step                        <"1". Performs step 1>
          
          --basecall_path               <path to base directory containing all fast5 and/or pod5 files you want to basecall.
                                        It will automatically separate samples based on naming conventions and 
                                        directory structure. example: /sequencing_run/">

          --basecall_speed              <"fast", "hac", "sup". Default = "hac">

          --basecall_mods               <Comma separated list of base modifications you want to basecall. See dorado docs 
                                        for more information. Example: "5mCG_5hmCG,5mC_5hmC,6mA". Default: "False">

          --basecall_compute            <"gpu", "cpu". Default: "gpu". Allows users to choose to basecall with 
                                        CPU or GPU. CPU basecalling is super slow and should only be used for small
                                        test datasets when GPUs are not available. Also, use --basecall_speed "fast"
                                        when basecalling with CPU to make it less slow. Default: "gpu">

          --basecall_config             <configuration name for basecalling setting. This is not necessary since dorado 
                                        is able to automatically determine the appropriate configuration. When set to "None"
                                        the basecaller will automatically pick the basecall configuration.
                                        Example: "dna_r9.4.1_e8_hac@v3.3". Default: "None">

          --basecall_trim              <"all", "primers", "adapters", "none". Default: "none". Note that as of December 2023
                                        dorado always trims the adapters during basecalling for directRNA data>

          --qscore_thresh              <Mean quality score threshold for basecalled reads. Default: 9>

          
          --demux                       <"True", "False". Whether you want the data to be demultiplexed
                                        setting it to "True" will perform demultiplexing>

          --trim_barcodes               <"True", "False". Only relevant is --demux is set to "True". 
                                        if set to "True" barcodes will be trimmed during demultiplexing
                                        and will not be present in output "fastq" files>

          --prefix                      <Will add a prefix to the beggining of your filenames, good
                                         when wanting to keep track of batches of data.
                                         Example: "Batch_1". Default value is "None" which does not add any prefixes>

          --out_dir                    <Name of output directory. Output files/directories will be output to
                                        "./results/<out_dir>/". Default: "output_directory">



## Pipeline parameters for STEP 2 (QC, Alignment, and Bambu Pre-processing)

          --step                        <"2". Performs step 2>

          --ont_reads_fastq             <path to fastq sequencing data, can submit multiple at once using pattern "/path/*.fastq".  If you don't specify this
                                        parameter the pipeline will not run. An exception is if you specify a  ".bam" file input, see Pipeline parameters for
                                        STEP 2 from BAM below. Default: "None">
          
          --ont_reads_txt               <path to sequencing summary files, can submit multiple at once. Make sure they follow same naming pattern as fastq
                                        files. If a fastq file is named "sample_1.fastq" then its sequencing summary file should be "sample_1.txt"
                                        or "sample_1_sequencing_summary.txt". If fastq files and sequencing summary files don't follow the same naming pattern
                                        they might be paired erroneously, which will lead to erronous results from quality control steps. If you don't specify 
                                        this parameter pipeline will run, but skip PycoQC quality control step. Default: "None">

          --path                        <path to ONT files as exported by the basecaller from the PromethION, GridION, of MinION devices. It takes a path to one or 
                                        multiple samples. For example if you had two samples in "/dir1/dir2/sample_1/" and "/dir1/dir2/sample_2/" you would pass 
                                        "/dir1/dir2/" to have both samples be processed. Using the "path" parameter will make the pipeline unzip all the smaller 
                                        ".fastq.gz" files for each sample and concatenate them into one ".fastq" file for each sample. This parametes will also
                                        automatically search for the "sequencing_summary.txt" files for each sample, this file must be present for all samples
                                        for this parameter to work properly. If you pass this path parameter, you don't have to pass the "ont_reads_txt" and the
                                        "ont_reads_fastq" parameters, in fact those will get ignored even if you pass them. The goal of this parameter is to
                                        make it easier on the user so that you don't have to manually unzip and concatenate files. It also automatically 
                                        matches the naming pattern for the fastq and the sequencing summary file. This parameter assumes that a directory with
                                        the name of the sample is two directories below the "fastq_pass" directory, for example:
                                        "/dir1/dir2/sample_1/unique_flowcell_id/fastq_pass/". The pipeline will still work if this is not the case and each
                                        sample file will still have a unique ID, but the unique sample ID will be preappended with the whatever name the 
                                        directory two steps below the "fastq_pass" directory is named. For example, if your directory structure is:
                                        "/dir1/dir2/dir3/unique_flowcell_id/fastq_pass/", your files will be named "dir3_unique_flowcell_id.fastq"
                                        and "dir3_unique_flowcell_id.txt" (sequencing summary file). The ".fastq.gz" files must be in a folder
                                        named "fastq_pass" and the sequencing_summary.txt file must be in the same directory as the
                                        "fastq_pass" directory... Not inside the "fastq_pass" directory, but in the same directory as
                                        the "fastq_pass". Default: "None">
          
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
                                        "PCS110", "PCS111", "PCS114", "LSK114". Default: "PCS111">
  
          --is_chm13                    <logical, set to "True" if using CHM13 and "False" if not. Fixes CHM13 annotation for compatibility with Bambu and
                                        converts to ".gtf". Default: "False">
          
          --housekeeping                <path to bed file containing coordinates for housekeeping genes used by RSeQC quality control tool. Currently only
                                        supported if --is_chm13 is set to "False". You can find the bed file with housekeeping gene coordinates for GRCh38
                                        here: https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38.HouseKeepingGenes.bed.gz/download
                                        Default: "None">
                               
            --mapq                      <Integer, set it to the number you want to be used to filter ".bam" file by mapq. --mapq 10 filters out reads with MAPQ
                                        < 10. set it to 0 if don't want to filter out any reads. Default: 0>
                              
            --is_dRNA                   <Logical, set to "True" if you want to run the minimap2 step with the parameters for "noisy" Nanopore directRNAseq reads.
                                         Set to "False" if you want to run it with minimap2 step for PCR amplified Nanopore cDNA reads. If you set this parameter
                                         to "True" you should omit the "--cdna_kit" parameter as that is only used for the Pychopper step of the pipeline and that
                                         step is skipped for directRNAseq. Note that this is the only step in the pipeline that is specific for 
                                         DirectRNAseq, Step 2 from BAM and Step 3 are not modified for dRNA vs cDNA. Step 1 will need the 
                                         specific basecalling configuration for dRNA instead of cDNA if you are running dRNA. Default: "False">

          --trim_dRNA                   <Logical, set to "True" if you want to trim dRNA adapters using porechop_ABI.
                                        It will automatically find the adapters in your dRNA adapters, see this link: 
                                        https://github.com/bonsai-team/Porechop_ABI/issues/19 for caveats.
                                        Default: "False">
                                        

            --contamination_ref          <Path to a ".fasta" reference containing any contaminants you wish to align the unmapped reads against. For a 
                                          reference containing the contaminants we use in house got to the "contamination_reference_doc" folder in
                                          this GitHub or to this link on Zenodo: https://zenodo.org/deposit/8350277. If you pass a contamination
                                          reference a report will be generated with multiQC showing what percent of reads mapped to a contaminant.
                                          For further detail go into the "contamination_report/" folder under your pipeline execution results
                                          to find out which were the most abundant contaminants in each of your samples.>
            
            --quality_score              <Minimum mean base quality from basecalled sequence. For example, when set at 9 all reads with 
                                          mean base quality below 9 will be filtered out from the fastq file. This is executed by Pychopper,
                                          therefore it only affects cDNA data analysis. This parameter will always be ignored by dRNA pipeline
                                          execution, there is no quality fileter applied by the pipeline for dRNA. Default: 9>

          --track_reads                  <logical, set to "True" if you want Bambu to keep track of read assignments to transcripts in the output ".RDS" file
                                          from Bambu. Set to "False" if you don't need to keep track of read assignments (smaller files). Default: "False">

          --prefix                      <Will add a prefix to the beggining of your filenames, good
                                         when wanting to keep track of batches of data.
                                         Example: "Batch_1". Default value is "None" which does not add any prefixes>
  
 


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
  
          --is_chm13          <logical, set to "True" if using CHM13 and "False" if not. Fixes CHM13 annotation for compatibility with Bambu and
                              converts to ".gtf". Default: "False">
          
          --track_reads      <logical, set to "True" if you want Bambu to keep track of read assignments to transcripts in the output ".RDS" file
                              from Bambu. Set to "False" if you don't need to keep track of read assignments (smaller files). Default: "False">
                              
          --prefix            <Will add a prefix to the beggining of your filenames, good
                              when wanting to keep track of batches of data.
                              Example: "Batch_1". Default value is "None" which does not add any prefixes>
            
          --mapq            <integer, set it to the number you want to be used to filter ".bam" file by mapq. --mapq 10 filters out reads with
                              MAPQ < 10. Set it to 0 if don't want to filter out any reads. Default: 0>

## Pipeline parameters for STEP 3
  
          --bambu_rds         <path to individually pre-processed bambu RDS objects (output from step 1). Default: "None"
          
          --ref               <path to reference/assembly ".fa" file. if using ERCC make sure to concatenate it to the end of the file.
                              Default:"None">

          --annotation        <path to reference annotation ".gtf" file for GRCh38 or ".gff3" for CHM13. If using GRCh38 and ERCC concatenate the two
                              ".gtf" files prior to running the pipeline. If using ERCC with CHM13 make sure to set --is_chm13 to "True" and set
                              --ercc to the path of your ERCC ".gtf" file. Default: "None">

         --is_chm13          <logical, set to "True" if using CHM13 and "False" if not. Fixes CHM13 annotation for compatibility with Bambu and
                              converts to ".gtf". Default: "False">

          
          --fai               <path to reference index ".fai" file>

  
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

          --intermediate_qc    <path to directory containing intermediate quality control output from pipeline step 2.
                               Example: "/path/intermediate_qc_reports/" Default: "None">
            
          --multiqc_input    <path to directory containing multiqc input data from pipeline step 2. Use path and add **.
                              Example: /path/multiqc_input/** Default: "None">
          
          --multiqc_config   <path to multiqc ".yaml" config file. Default: "None". PS: You can find an example of a multiqc config file that works
                              for this pipeline under /cDNA_pipeline/workflow/bin/multiqc_config.yaml>

## Examples of the submissions

### Example for Step 1: Basecalling 
          nextflow ../main.nf --step 1 --basecall_path ../../test_dRNA_data/basecalling_2/" \
              --basecall_speed "fast" \
              --basecall_mods "False" \
              --basecall_config "False" \
              --basecall_trim "none" \
              --basecall_compute "cpu" \
              --basecall_demux "True" \
              --outdir "test_basecall_dRNA_1"

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

### Example for step 2 (cDNA): GRCh38 with ERCCs and contamination step

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
              --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
              --contamination_ref "../../references/contamination_ref.fasta"
              

### Example for step 2 (cDNA): GRCh38 with ERCCs using path parameter instead of --ont_reads_fq and --ont_reads_txt

          nextflow ../main.nf --step 2 \ 
              --path "../ont_data/" \
              --ref "../../references/Homo_sapiens.GRCh38_ERCC.fa" \
              --annotation "../../references/Homo_sapiens.GRCh38.106_ERCC.gtf" \
              --out_dir "./GRCh38_ERCC_test/" \
              --cdna_kit "PCS111" \
              --is_chm13 "False" \
              --track_reads "False" \
              --mapq "0" \
              --housekeeping "../../references/hg38.HouseKeepingGenes.bed"
              
### Example for step 2 (dRNA): GRCh38 with ERCCs

          nextflow ../main.nf --ont_reads_fq "../../data/ont_data/2023-06-06_brain_directRNA_intronic_reads/*.fastq" \
                    --ont_reads_txt "../..//data/ont_data/2023-06-06_brain_directRNA_intronic_reads/*.txt" \
                    --ref "/../../references/Homo_sapiens.GRCh38_ERCC.fa" \
                    --annotation "../../references/Homo_sapiens.GRCh38.107_ERCC.gtf" \
                    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
                    --out_dir "./GRCh38_dRNA_EERCC_test/" \
                    --bambu_track_reads "True" \
                    --is_dRNA "True" \
                    --trim_dRNA "True" \
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
              --fai "./results/GRCh38_ERCC_test/fai/*.fai" \
              --annotation "../../references/Homo_sapiens.GRCh38.106_ERCC.gtf" \
              --is_discovery "True" \
              --track_reads "False" \
              --NDR "auto" \
              --multiqc_input "./results/GRCh38_ERCC_test/multiQC_input/**" \
              --multiqc_config "../../references/multiqc_config.yaml" \
              --out_dir "./GRCh38_ERCC_test/" \
              --intermediate_qc "./results/GRCh38_ERCC_test/intermediate_qc_reports/" \
              --is_chm13 "False"       


### Example of using resume function to continue a halted submission

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
              --is_chm13 "False" -resume
