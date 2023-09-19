#!/bin/bash

singularity exec ../singularity_container/quality_control.sif multiqc -c ../references/multiqc_config.yaml -n multiQC_report.html .
