#!/bin/sh
#SBATCH --time=10:00:00  --mem=5GB

variant_bed_dir="$1"
input_individuals="$2"
all_variants_gene_mapping_file="$3"
gencode_hg19_gene_annotation_file="$4"
gene_mapping_method="$5"


date
python variant_to_gene_mapping.py $variant_bed_dir $input_individuals $all_variants_gene_mapping_file $gencode_hg19_gene_annotation_file $gene_mapping_method
date