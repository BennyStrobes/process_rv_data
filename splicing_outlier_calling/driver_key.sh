####Ben Strober
####3/28/16
#########These scripts will:
##############1. Generate tissue specific junction files using RAIL aligned split reads
##############2. Extract rareness of each of our observed junctions using SNAPTRON
##############3. Call outliers using DM-MD model



#############################################################
#Input data
#############################################################

#V6 sample attribute file that contains 'flagged samples'
v6_sample_attribute_file="/work-zfs/abattle4/lab_data/GTEx_v6p/sample_annotations/GTEx_Data_2014-06-13_Annotations_SampleAttributesDS.txt"

#List of gtex tissues. First colum is traditional GTEx tissue name. Second column is GTEx tissue name as listed in $sample_attribute_file
tissue_list_input_file="/work-zfs/abattle4/lab_data/GTEx_v6p/tissue_names/sr_tissue_converter.txt"

#File that contains all jxns in found in GTEx. Note: it is aligned to hg38 using RAIL.
#Each row is a jxn and column 11 contains all samples and their corresponding counts of the given junction
#Downloaded from snaptron data website: "http://snaptron.cs.jhu.edu/data/gtex/junctions.bgz"
snaptron_gtex_junction_file="/work-zfs/abattle4/lab_data/snaptron_data/gtex_hg38_jxns.bgz"

#File that will be used for conversion from SNAPTRON IDs to GTEx sample IDs
#Downloaded from snaptron website: "http://snaptron.cs.jhu.edu/data/gtex/samles.tsv"
snaptron_gtex_samples_file="/work-zfs/abattle4/lab_data/snaptron_data/gtex_samples.tsv"


#Directory that contains Leafcutter code (https://github.com/davidaknowles/leafcutter)
leafcutter_code_dir="/work-zfs/abattle4/bstrober/tools/leafcutter/clustering/"

#Gencode hg19 gene annotation file
gencode_hg19_gene_annotation_file="/work-zfs/abattle4/lab_data/hg19/gencode_gene_annotations/gencode.v19.annotation.gtf.gz"

#Directory that contains necessary liftover information.
##Specifically, it must contain:
#####1. 'liftOver'   --> the executable
#####2. 'hg19ToHg38.over.chain.gz'   --> for converting from hg19 to hg38
#####2. 'hg38ToHg19.over.chain.gz'   --> for converting from hg38 to hg19
liftover_directory="/work-zfs/abattle4/bstrober/tools/liftOver_x86/"

#v6p covariate directory. Also contains samples used in v6p analysis
covariate_directory_v6="/work-zfs/abattle4/lab_data/GTEx_v6p/covariates/"


#############################################################
#Parameters
#############################################################
#We only consider junctions that have at least one sample with greater than or equal to $min_reads
min_reads="10"
#Some genes have over 100 hundred jxns. Need to filter this in some way.
#Restrict to a maximum of $max_dm_junctions junctions for dm outlier calling 
max_dm_junctions="20"
#How to Select the $max_dm_junctions.
#extreme_junctions is the method of taking the $max_dm_junctions/2 with the highest read count. and the $max_dm_junctions/2 with the smallest read count 
#ignore_genes is the method of ignoring any genes that have more than max_dm_junctions
jxn_filter_method="ignore_genes"
#How to deal with covariates when fitting Dirichlet multinomial
#covariate_regression is the method of fitting GLM DM with covariates; then using integer residuals as 'corrected counts'
#stan_pca_regression, remove_top, and v6p_eqtl_covariates
covariate_regression_method="v6p_eqtl_covariates"
##Number of principle components to use
num_pc="10"
##Method to filter "global" splicing outlier individuals
#'only_flagged' removes all samples flagged in v6 sample annotation file
#'filter1' removes all samples flagged in v6 sample annotation file as well as removing global outliers in a specific tissue. Global outliers were removed by:
######a. Removing samples greater than 1 stdev away from the total number of nonzero junctions
######b. Removing samples greater than 1 stdev away from the total number of unique junctions
#'filter2' removes all samples flagged in v6 sample annotation file, all samples we do not have covariates for, and global outliers in a specific tissue. Global outliers were removed by:
######a. Removing samples greater than 1 stdev away from the total number of nonzero junctions
######b. Removing samples greater than 1 stdev away from the total number of unique junctions
filter_global_outlier_method="filter2"

#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################
# Root directories that all other directories will be placed in 
splicing_output_root="/work-zfs/abattle4/bstrober/rare_variant/rare_splice/process_rv_data/splicing_outlier_calling/"
#output_dir for pre_process.sh
pre_process_output_dir=$splicing_output_root"pre_process/"
#output_dir for pre_process.sh
filter_pre_process_output_dir=$splicing_output_root"filter_pre_process/"
#output_dir for generate_junctions.sh (contains junction files / sample)
junctions_output_dir=$splicing_output_root"junctions/"
#Another output_dir for generate_junctions.sh (contains direct output of leafcutter)
clusters_output_dir=$splicing_output_root"clusters/"
#Another output_dir for generate_junctions.sh (contains filtered versions of leafcutter output)
clusters_filter_output_dir=$splicing_output_root"clusters_filter/"
#Output_dir for call_outlier_dm.sh
outlier_calling_dm_output_dir=$splicing_output_root"outlier_calling_dm/"
#Output dir that contains a covariate file for each gtex tissue type
covariate_output_dir=$splicing_output_root"covariates/"
# Output dir that contains samples used in outlier calling analysis
outlier_calling_samples_dir=$splicing_output_root"outlier_calling_samples/"










#################################################################################################################
#Sample filtering
##################################################################################################################
#pre_process.sh produces lists of what gtex samples will be used in which tissue
#Output files made (all we be written in pre_process_output_dir)
##1. *_rnaseq_sample_ids.txt where * is the conventional_tissue_name. These files are simply a list of all GTEx RNASEQ samples used for this tissue.
##2. all_individuals.txt. A list of all of the individuals across all tissues. (going to ignore this for now. Because we are limiting to individuals with RNA seq samples)
##3. tissue_sample_counts.txt. A table that contains information on how many gtex samples are in each tissue. Column 0 is the gtex tissue name and column 1 is the number of samples in that tissue.
#NOTE: By samples, I mean RNA-seq samples. By individuals, I mean actual people. So for gtex 1 individual may have multiple samples. But 1 sample only has 1 individual. 
#Runtime: < 1 min
###################################################################################################################
if false; then
sh pre_process.sh $v6_sample_attribute_file $tissue_list_input_file $pre_process_output_dir $covariate_directory_v6
fi


#################################################################################################################
#Temporary (Remove once we figure out issue of why there are 4 samples used in gtex not in snaptron)
##################################################################################################################
#filter_pre_process.sh produces lists of what gtex samples will be used in which tissue.
#This list is simply a filtering of the above step. It is done seperately for modularity purposes
#The filtering filters samples/individuals to those that SNAPTRON has RNA seq samples for (ie. not all of v6p  --> There are 4 SAMPLES NOT In snaptron)
#Output files made (all we be written in filter_pre_process_output_dir)
##1. *_rnaseq_sample_ids.txt where * is the conventional_tissue_name. These files are simply a list of all GTEx RNASEQ samples used for this tissue.
##2. all_individuals.txt. A list of all of the individuals across all tissues. (going to ignore this for now. Because we are limiting to individuals with RNA seq samples)
##3. tissue_sample_counts.txt. A table that contains information on how many gtex samples are in each tissue. Column 0 is the gtex tissue name and column 1 is the number of samples in that tissue.
#NOTE: By samples, I mean RNA-seq samples. By individuals, I mean actual people. So for gtex 1 individual may have multiple samples. But 1 sample only has 1 individual. 
#Runtime: < 1 min
###################################################################################################################
if false; then
sh filter_pre_process.sh $pre_process_output_dir $snaptron_gtex_samples_file $filter_pre_process_output_dir $tissue_list_input_file
fi


##################################################################################################################
#JXN PREPERATION
##################################################################################################################
#generate_junctions.sh produces the jxn files, filters the jxn files, converts to hg19, and maps to genes
#As such, this script has multiple subparts
###########################################
##Subpart_1: 'generate_junctions.py'
####This part produces a junction file for each gtex rna-seq sample
####Ouput files made (all will be written to #junction_output_dir):
#######1. *.junc where * is a GTEx sample ID. This file is the exact format as the output of step 1 of leafcutter
#######2. *_junction_file_locations.txt where * is a tissue type. Each line is the absolute location of all of the &.junc files in this tissue. Where & is a GTEx sample id.
#spot checked
############################################
##Subpart_2: 'leafcutter_cluster.py'
####This script was taken from the leafcutter github (https://github.com/davidaknowles/leafcutter)
####It is precisely step 2 of leafcutter (ie it will generate clusters of jxns from the jxns files)
#Output files made (all will be written in $clusters_output_dir):
##1. *_leafcutter_hg38_perind_numers.counts.gz  where * is gtex tissue type 
############################################
##Subpart_3: 'filter_clusters.py '
#A
##This script does a number of different cleaning/organization related activities. It works directly on the leafcutter output (*_leafcutter_hg38_perind_numers.counts.gz)
##It first filters junctions if:
#####1. junction contains no individuals with more than min_reads mapped to the junction
#####2. Junction is on a non-autosomal chromosome
##After these filters, remove jxns that belong to a cluster with only 1 junction.
#NOTE: There is some discongruence amongst different methods on what defines a jxn position.
##Snaptron defines a jxn as (a,b). Running leafcutter on that result yields (a,b+1). Derek's leafcutter yields (a-1,b+1).
#This step will convert into Derek's format!
#Also change header to alphabetically ordered individual ids as opposed to sample ids
##output saved as *_hg38_filtered_&_reads.txt where * is tissue type and & is min_reads. Saved in 'clusters_filter_output_dir'
#B
##It then converts from hg38 --> hg19 using Chris Wilk's liftover parameters (#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
##As some jxns are not 'liftoverable', a small fraction of jxns don't map. 
#Therefor, remove jxns that belong to a cluster with only 1 junction.
##output saved as *_hg19_filtered.txt where * is tissue type. Saved in 'clusters_filter_output_dir'
##################################################################################################################
if false; then
while read tissue_type alternative_name; do
	sbatch generate_junctions.sh $tissue_type $filter_pre_process_output_dir $snaptron_gtex_junction_file $snaptron_gtex_samples_file $leafcutter_code_dir $junctions_output_dir $clusters_output_dir $clusters_filter_output_dir $min_reads $liftover_directory $gencode_hg19_gene_annotation_file
done<$tissue_list_input_file
fi
##################################################################################################################
#generate_junctions_cross_tissues.sh performs the final steps of junction file generation that must be done by taking all tissues into account (as we want clusters to be the same across tissues)
#As such, this script has multiple subparts:
##########################
##Subpart_1: generate_cross_tissue_clusters.py
#####This section creates leafcutter clusters ACROSS TISSUES based on our filtered junctions.
#####It will remove junctions to belong to clusters with only one junction (It will do this in a tissue specific and global manner. So really only tissue specific)
#####It will output a file of the form $clusters_filter_output_dir$tissue_type"_hg19_filtered_xt_reclustered.txt"
##Subpart_2: map_clusters_to_genes_cross_tissues.py
#####This section will map a cluster to a given gencode v19 protein coding gene
#####A cluster is mapped to a gene if one of its junctions lies within the gene body
#####We will filter out clusters that are not mapped to any genes
#####It will output a file of the form $clusters_filter_output_dir$tissue_type"_hg19_filtered_xt_reclustered_gene_mapped.txt"
#####It will also output a file of the form $clusters_filter_output_dir"cluster_info.txt"
if false; then
sbatch generate_cross_tissue_clusters.sh $tissue_list_input_file $clusters_filter_output_dir $gencode_hg19_gene_annotation_file
fi





##################################################################################################################
#Covariate Preperation
##################################################################################################################
num_pc="3"
covariate_regression_method="v6p_eqtl_covariates"

covariate_regression_method="junction_pc"
tissue_type="Adipose_Subcutaneous_Analysis"
if false; then
sbatch generate_covariate_files.sh $tissue_type $covariate_output_dir $covariate_directory_v6 $covariate_regression_method $num_pc $clusters_filter_output_dir
fi






#################################################################################################################
#DM OUTLIER CALLING
##################################################################################################################
total_nodes="50"
tissue_type="Adipose_Subcutaneous_Analysis"
covariate_regression_method_data="junction_pc"
tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_xt_reclustered_gene_mapped.txt"
rna_seq_samples_file=$filter_pre_process_output_dir$tissue_type"_rnaseq_sample_ids.txt"
outlier_calling_samples_file=$outlier_calling_samples_dir$tissue_type"_"$filter_global_outlier_method"_used_samples.txt"

covariate_regression_method="none"
num_pc="0"
covariate_file=$covariate_output_dir$tissue_type"_"$covariate_regression_method_data"_"$num_pc"_covariates.txt"
if false; then
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    tissue_specific_outlier_file=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$node_number
    sbatch call_outlier_dm.sh $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $covariate_file $rna_seq_samples_file $num_pc $outlier_calling_samples_file
done
fi

covariate_regression_method="junction_pc_no_regress_reg"
num_pc="1"
covariate_file=$covariate_output_dir$tissue_type"_"$covariate_regression_method_data"_"$num_pc"_covariates.txt"
if false; then
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    tissue_specific_outlier_file=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$node_number
    sbatch call_outlier_dm.sh $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $covariate_file $rna_seq_samples_file $num_pc $outlier_calling_samples_file
done
fi

covariate_regression_method="junction_pc_no_regress_reg"
num_pc="3"
covariate_file=$covariate_output_dir$tissue_type"_"$covariate_regression_method_data"_"$num_pc"_covariates.txt"
if false; then
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    tissue_specific_outlier_file=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$node_number
    sbatch call_outlier_dm.sh $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $covariate_file $rna_seq_samples_file $num_pc $outlier_calling_samples_file
done
fi


covariate_regression_method="junction_pc_only_no_regress_reg"
num_pc="1"
covariate_file=$covariate_output_dir$tissue_type"_"$covariate_regression_method_data"_"$num_pc"_covariates.txt"
if false; then
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    tissue_specific_outlier_file=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$node_number
    sbatch call_outlier_dm.sh $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $covariate_file $rna_seq_samples_file $num_pc $outlier_calling_samples_file
done
fi

covariate_regression_method="junction_pc_only_no_regress_reg"
num_pc="3"
covariate_file=$covariate_output_dir$tissue_type"_"$covariate_regression_method_data"_"$num_pc"_covariates.txt"
if false; then
for node_number in $(seq 0 `expr $total_nodes - "1"`); do
    tissue_specific_outlier_file=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$node_number
    sbatch call_outlier_dm.sh $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $covariate_file $rna_seq_samples_file $num_pc $outlier_calling_samples_file
done
fi

if false; then
total_nodes="50"
tissue_type="Adipose_Subcutaneous_Analysis"
covariate_regression_method="junction_pc_no_regress_reg"

num_pc="1"
tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_xt_reclustered_gene_mapped.txt"
tissue_specific_outlier_root=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"
sbatch merge_outlier_calling_results.sh $tissue_type $total_nodes $tissue_specific_outlier_root $tissue_specific_jxn_file

num_pc="3"
tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_xt_reclustered_gene_mapped.txt"
tissue_specific_outlier_root=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"
sbatch merge_outlier_calling_results.sh $tissue_type $total_nodes $tissue_specific_outlier_root $tissue_specific_jxn_file


total_nodes="50"
tissue_type="Adipose_Subcutaneous_Analysis"
covariate_regression_method="junction_pc_only_no_regress_reg"

num_pc="1"
tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_xt_reclustered_gene_mapped.txt"
tissue_specific_outlier_root=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"
sbatch merge_outlier_calling_results.sh $tissue_type $total_nodes $tissue_specific_outlier_root $tissue_specific_jxn_file

num_pc="3"
tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_xt_reclustered_gene_mapped.txt"
tissue_specific_outlier_root=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"
sbatch merge_outlier_calling_results.sh $tissue_type $total_nodes $tissue_specific_outlier_root $tissue_specific_jxn_file

total_nodes="50"
tissue_type="Adipose_Subcutaneous_Analysis"
covariate_regression_method="none"

num_pc="0"
tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_xt_reclustered_gene_mapped.txt"
tissue_specific_outlier_root=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"
sbatch merge_outlier_calling_results.sh $tissue_type $total_nodes $tissue_specific_outlier_root $tissue_specific_jxn_file
fi


