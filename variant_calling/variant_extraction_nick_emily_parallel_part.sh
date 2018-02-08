#!/bin/sh
#SBATCH --time=28:00:00 



vcf_file="$1"
wgs_samples_dir="$2"
job_number="$3"
total_jobs="$4"


bash vcf2bedfiles_parallel_part.sh $vcf_file $wgs_samples_dir $job_number $total_jobs
