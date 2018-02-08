#!/bin/sh
#SBATCH --time=10:00:00 --mem=10GB

tissue_type="$1"
pre_process_dir="$2"
snaptron_gtex_junction_file="$3"
snaptron_gtex_samples_file="$4"
leafcutter_code_dir="$5"
junctions_output_dir="$6"
clusters_output_dir="$7"
clusters_filter_output_dir="$8"
min_reads="$9"
liftover_directory="${10}"
gencode_hg19_gene_annotation_file="${11}"

date 
cwd=`pwd`"/"
####################################################################################
#Generate Junctions
####################################################################################
echo "generating junctions"
python generate_junctions.py $tissue_type $pre_process_dir $snaptron_gtex_junction_file $snaptron_gtex_samples_file $junctions_output_dir


#####################################################################################
#Generate Clusters
#####################################################################################
echo "generating clusters"
leafcutter_input_file=$junctions_output_dir$tissue_type"_junction_file_locations.txt"
leafcutter_output_prefix=$tissue_type"_leafcutter_hg38"
leafcutter_output_dir=$clusters_output_dir
cd $leafcutter_code_dir
#All of the following parameters were selected to be super conservative (ie they leave the most junctions... Wouldn't wanna loose any rare ones)
#m is minimum reads in a cluster (default 30 reads)
#l is maximum intron length in bp (default 100,000bp)
#p is minimum fraction of reads in a cluster that support a junction (default 0.001)
python leafcutter_cluster.py -j $leafcutter_input_file -m 1 -o $leafcutter_output_prefix -l 500000 -r $leafcutter_output_dir -p 0.00000000000001


leafcutter_output_file=$clusters_output_dir$tissue_type"_leafcutter_hg38_perind_numers.counts.gz"

#####################################################################################
#Filter Clusters
#####################################################################################
echo "filtering clusters"
cd $cwd
python filter_clusters.py $tissue_type $leafcutter_output_file $clusters_filter_output_dir $min_reads $liftover_directory