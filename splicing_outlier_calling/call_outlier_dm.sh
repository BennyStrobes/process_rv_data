#!/bin/sh
#SBATCH --time=10:00:00 --mem=15GB

tissue_type="$1"
tissue_specific_jxn_file="$2"
tissue_specific_outlier_file="$3"
outlier_calling_dm_output_dir="$4"
max_dm_junctions="$5"
jxn_filter_method="$6"
node_number="$7"
total_nodes="$8"
covariate_regression_method="$9"
covariate_file="${10}"
rna_seq_samples_file="${11}"
num_pc="${12}"
outlier_calling_samples_file="${13}"


date
echo $node_number
echo $num_pc
python call_outlier_dm.py $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $tissue_list_input_file $covariate_file $rna_seq_samples_file $num_pc $outlier_calling_samples_file $lam
date