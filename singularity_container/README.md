# Singularity Containers




2023-12-08_bambu.def - definition file for singularity container used to run the Bambu R package on the Nextflow Pipeline.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/bambu:sha256.84a7c3d4829374db2099567ef0d2d3faa703350f016b4045204ce371499eafb1`


2023-12-08_dorado.def - definition file for singularity container used to run guppy basecaller... Was not used due to data being basecalled on the PromethION.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/dorado:sha256.cb7a2872ced6cdbdaf953931867d701bf02f8a89835937dd86f1d0a8904a2306`
 


2023-12-11_nanopore.def - definition file for singularity container used to run the software nanopore data analysis.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/nanopore:sha256.2f1ba1d359e5a58f6e25cf2f0c8376911465d8cf39342aa9d30d0ec47e7d8eb8`



2023-12-08_quality_control.def - definition file for singularity container used to run the quality control software.

`pull command: singularity pull --arch amd64 library://ebbertlab/nanopore_cdna/quality_control:sha256.05c71be73b714f0f45026e975563d315334f45a753a4119987c42d15b92e8d46`



### For more information about software versions see the %help section of definition (.def) files for the singularity images.
