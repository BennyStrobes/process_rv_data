#!/bin/sh
#SBATCH --time=16:00:00

tissue_type="$1"
tissue_specific_jxn_file="$2"
tissue_specific_rare_outlier_file_root="$3"
snaptron_directory="$4"
liftover_directory="$5"
output_dir="$6"





hg38_jxn_rareness_file=$tissue_specific_rare_outlier_file_root"hg38.txt"

#for each jxn in $tissue_specific_jxn_file, learn how many samples have that jxn
python call_outlier_rare_extract.py $tissue_type $tissue_specific_jxn_file $hg38_jxn_rareness_file $snaptron_directory


hg19_jxn_rareness_file=$tissue_specific_rare_outlier_file_root"hg19.txt"
python call_outlier_rare_liftover.py $hg38_jxn_rareness_file $hg19_jxn_rareness_file $liftover_directory $output_dir $tissue_type
