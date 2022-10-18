#!/bin/bash


./merged_first_half_loose_discovery_GRCH38-107.sh 
./merged_first_half_moderate_discovery_GRCH38-107.sh 
./merged_first_half_stringent_discovery_GRCH38-107.sh
./merged_first_half_super_stringent_discovery_GRCH38-107.sh

./merged_second_half_loose_discovery_GRCH38-107.sh
./merged_second_half_moderate_discovery_GRCH38-107.sh
./merged_second_half_stringent_discovery_GRCH38-107.sh
./merged_second_half_super_stringent_discovery_GRCH38-107.sh

ssh mcc-dtn.ccs.uky.edu

cd /scratch/bag222/cDNA_pipeline/workflow/submission_scripts/

./rclone_transfer.sh
