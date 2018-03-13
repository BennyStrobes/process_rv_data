#!/bin/sh
#SBATCH --time=20:00:00 --mem=10GB


vcf_file="$1"
subject_phenotype_file="$2"
rna_seq_individuals="$3"
wgs_samples_dir="$4"

#****************************************
python filter_subjects.py $subject_phenotype_file $rna_seq_individuals $wgs_samples_dir"EA_gtex_ids.txt"


sh vcf2bedfiles_serial_part.sh $vcf_file $wgs_samples_dir
