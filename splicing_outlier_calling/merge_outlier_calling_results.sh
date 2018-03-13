#!/bin/sh
#SBATCH --time=1:00:00


tissue_type="$1"
total_nodes="$2"
tissue_specific_outlier_root="$3"
tissue_specific_jxn_file="$4"

python merge_outlier_calling_results.py $tissue_type $total_nodes $tissue_specific_outlier_root $tissue_specific_jxn_file


python merge_outlier_calling_jxn_results.py $tissue_type $total_nodes $tissue_specific_outlier_root