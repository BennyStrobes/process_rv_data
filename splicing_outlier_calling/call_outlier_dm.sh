#!/bin/sh
#SBATCH --time=24:00:00

tissue_type="$1"
tissue_specific_jxn_file="$2"
tissue_specific_outlier_file="$3"
outlier_calling_dm_output_dir="$4"
max_dm_junctions="$5"
jxn_filter_method="$6"
node_number="$7"
total_nodes="$8"
covariate_regression_method="$9"
sample_attribute_file="${10}"
tissue_list_input_file="${11}"
covariate_file="${12}"
rna_seq_samples_file="${13}"
num_pc="${14}"
filter_global_outlier_method="${15}"
v6_sample_attribute_file="${16}"
outlier_calling_samples_file="${17}"


date
python call_outlier_dm.py $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $sample_attribute_file $tissue_list_input_file $covariate_file $rna_seq_samples_file $num_pc $filter_global_outlier_method $v6_sample_attribute_file $outlier_calling_samples_file
date