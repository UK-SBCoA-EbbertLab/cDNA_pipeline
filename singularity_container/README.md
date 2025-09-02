# Singularity Containers




2024-04-24_bambu.def - definition file for singularity container used to run the Bambu R package on the Nextflow Pipeline.

pull command: `singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/bambu:sha256.44e2b6d7282a488b95b132198b7c4ca659c9e8d6a83797493e746aa3a87ecfea`


2024-03-18_dorado.def - definition file for singularity container used to run guppy basecaller... Was not used due to data being basecalled on the PromethION.

pull command: `singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/dorado:sha256.4cec2a7db51075e2480ac8b75a7b12c4e77727aa779ccb07925d77ed31283cbd`
 


2024-03-18_nanopore.def - definition file for singularity container used to run the software nanopore data analysis.

pull command: `singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/nanopore:sha256.df8b78d8644d58861e7ea5c82c218b76845559c0d71fdb45e58d271a349fd045`



2024-03-18_quality_control.def - definition file for singularity container used to run the quality control software.

pull command: `singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/quality_control:sha256.45fd0d3aa770abea5c195a734c139250f73af2763f8ae10e03c4751143844bb4`



### For more information about software versions see the %help section of definition (.def) files for the singularity images.
