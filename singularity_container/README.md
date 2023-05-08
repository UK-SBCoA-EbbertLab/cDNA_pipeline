# Singularity Containers




bambu.def - definition file for singularity container used to run the Bambu R package on the Nextflow Pipeline.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/bambu:latest`


guppy.def - definition file for singularity container used to run guppy basecaller... Was not used due to data being basecalled on the PromethION.

 `pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/guppy:latest`
 


nanopore.def - definition file for singularity container used to run the software nanopore data analysis.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/nanopore:latest`



quality_control.def - definition file for singularity container used to run the quality control software.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/quality_control:latest`



### For more information about software versions see the %help section of definition (.def) files for the singularity images.
