#!/bin/sh
#SBATCH --time=16:00:00

input_file="$1"
output_file="$2"

python call_outlier_rare_gtex.py $input_file $output_file