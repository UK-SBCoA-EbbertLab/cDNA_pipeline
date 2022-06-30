# cDNA_pipeline

## System requirements:

### 1) Have a functioning version of Nextflow in your Path.

          - On the MCC it can be found under `/project/mteb223_uksr/sequencing_resources/tools/bin/nextflow`
          
          - Make sure to run `module load ccs/java/jdk1.8.0_202` on load Java. Nextflow needs Java to work. I added this command to my `~/.bash_profile` to make life easier.
          
### 2) Get the singularity image used for this pipeline:

          - You can find it on the MCC under: `/project/mteb223_uksr/singularity_files/2022-06-30_cdna_nanopore_pipe.sif`
          
          - Alternatively you can build the singularity file from scratch using the ".def" file contained here. I do not recommend this as tools could have been updated in 
          and not be compatible with the pipeline anymore.
          
          
### 3) Clone this github repo using the command `git clone https://github.com/UK-SBCoA-EbbertLab/cDNA_pipeline`


### 4) Put the singularity ".sif" file in the `singularity_container` folder. You can use a softlink if preferred.

### 5) Go into the `workflows/nextflow.config` file and make any necessary changes:

        - Alter slurm job manager parameters and singularity file path to suit your needs and PI affiliation. I don't recommend changing the memory/cpu/time allocated 
        for the job manager.
        
        - You do not need to change the paths to the files specified by the `params.xxx` variables. Those can be set at the time of executing the pipeline.

          
### 6) Make sure you have all the files, "sequencing_summary.txt" files, reference genomes/assemblies files and annotation files you will need to run the pipeline.
          
          - ".fastq" 
          - "sequencing_summary.txt" --- You can create dummy files if you don't have these available for now. This option will be improved in the future.
          - refecence/assembly ".fa" file.
          - annotation ".gtf" file is preffered. Only use ".gff3" if using CHM13. Pipeline has an option to handle this.
          






## An example of the submission of the pipeline can be seen under `workflow/test_workflow/submit_cDNA_test_workflow.sh"
