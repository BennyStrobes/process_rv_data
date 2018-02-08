#!/bin/sh
#SBATCH --time=05:00:00


tissue_list_input_file="$1"
clusters_filter_output_dir="$2"
gencode_hg19_gene_annotation_file="$3"


echo "Re-generating clusters across tissues"
python generate_cross_tissue_clusters.py $tissue_list_input_file $clusters_filter_output_dir

echo "Mapping clusters to genes across tissues"

python map_clusters_to_genes_cross_tissues.py $tissue_list_input_file $clusters_filter_output_dir $gencode_hg19_gene_annotation_file