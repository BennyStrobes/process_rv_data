#!/bin/sh
#SBATCH --time=3:00:00

tissue_list_input_file="$1"
clusters_filter_output_dir="$2"
min_reads="$3"
outlier_calling_rare_gtex_output_dir="$4"


while read tissue_type alternative_name; do
echo $tissue_type
tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_"$min_reads"_reads.txt"
tissue_specific_rare_outlier_file=$outlier_calling_rare_gtex_output_dir$tissue_type"_hg19.txt"
sh call_outlier_rare_gtex.sh $tissue_specific_jxn_file $tissue_specific_rare_outlier_file

done<$tissue_list_input_file