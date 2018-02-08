#!/bin/bash

# Goes from GTEx VCF files to a bed file for each individual with that individual's variant sites and their allele frequency

set -o nounset

# set number of processes to run in parallel when gzipping
nproc=15

## takes as input the two vcfs to operate on
usage="usage: vcf2bedfiles.sh <GTEx SNV/indel vcf> <SV vcf>"
if [ $# -ne 2 ]; then
    echo $usage
    exit
fi

vcf=$1
wgs_samples_dir=$2

indincl=$wgs_samples_dir"EA_gtex_ids.txt"

# get prefix of output files
nEA=`wc -l $indincl | awk '{print $1}'`
fileprefix=`basename $vcf`
fileprefix=${fileprefix%.vcf.gz}
fileprefix=${fileprefix}"_"${nEA}"EAonly"

prefix=${wgs_samples_dir}${fileprefix}

## actually run things!
#######################
date
# first get vctools to generate useful information
echo "Processing VCF files with vcftools..."
#*************************************************************************************************
bash vcf2bedfiles_helper_processVCF.sh $vcf $indincl $prefix
echo "Processing VCF files (SNPs/indels) done."
date

# process SNPs
echo
echo "Processing SNPs..."
bash vcf2bedfiles_helper_processVCFtoolsOutput_serial_part.sh SNPs $prefix
sleep 5 # so they don't both create the outdir at the same time


