#!/bin/sh
#SBATCH --time=05:00:00 --mem=20GB


tissue_type="$1"
covariate_output_dir="$2"
covariate_directory_v6="$3"
covariate_regression_method="$4"
num_pc="$5"
clusters_filter_output_dir="$6"

python generate_covariate_files.py $tissue_type $covariate_output_dir $covariate_directory_v6 $covariate_regression_method $num_pc $clusters_filter_output_dir