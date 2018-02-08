#!/bin/sh
#SBATCH --time=2:00:00

sample_attribute_file="$1"
tissue_list_input_file="$2"
pre_process_output_dir="$3"
covariate_directory_v6="$4"

python pre_process.py $sample_attribute_file $tissue_list_input_file $pre_process_output_dir $covariate_directory_v6