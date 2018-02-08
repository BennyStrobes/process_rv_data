####Ben Strober
####2/8/18
#########These scripts will:
##############1. Generate tissue specific junction files using RAIL aligned split reads
##############2. Extract rareness of each of our observed junctions using SNAPTRON
##############3. Call outliers using DM-MD model



#############################################################
#Input data
#############################################################
#File containing which GTEx samples are to be used for each tissue. Filter list on val(column[27]) == 'RNASEQ'
sample_attribute_file="/work-zfs/abattle4/lab_data/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SampleAttributesDS.txt"

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

#File containing TPM scores for all RNA-seq samples in V7 for all genes
read_counts_tpm_file="/scratch1/battle-fs1/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt.gz"

#v6 covariate directory
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
#In performing outlier calling, we filter samples. This is not possible to do earlier on b/c it requires junction files. The final list of RNA-seq samples used is:
outlier_calling_samples_dir=$splicing_output_root"outlier_calling_samples/"

#################################################################################################################
##################################################################################################################
#pre_process.sh produces lists of what gtex samples will be used in which tissue
#Output files made (all we be written in pre_process_output_dir)
##1. *_rnaseq_sample_ids.txt where * is the conventional_tissue_name. These files are simply a list of all GTEx RNASEQ samples used for this tissue.
##2. all_rna_samples_file_locations.txt. This file is just a list of all the absolute file locations from part 1. Used for downstream analysis
##3. all_individuals.txt. A list of all of the individuals across all tissues. (going to ignore this for now. Because we are limiting to individuals with RNA seq samples)
##4. tissue_sample_counts.txt. A table that contains information on how many gtex samples are in each tissue. Column 0 is the gtex tissue name and column 1 is the number of samples in that tissue.
#NOTE: By samples, I mean RNA-seq samples. By individuals, I mean actual people. So for gtex 1 individual may have multiple samples. But 1 sample only has 1 individual. 
#Runtime: < 1 min
#Spot checked
###################################################################################################################
if false; then
sh pre_process.sh $sample_attribute_file $tissue_list_input_file $pre_process_output_dir
fi

#A file generated from pre_process.sh.
#Each line in the file corresponds to a given tissue's absolute location of a file that contains all of the gtex sample ids for that tissue
all_rna_samples_file_locations=$pre_process_output_dir"/all_rna_samples_file_locations.txt"

#################################################################################################################
#Temporary (remove once we get all v7 samples)
##################################################################################################################
#filter_pre_process.sh produces lists of what gtex samples will be used in which tissue.
#This list is simply a filtering of the above step. It is done seperately for modularity purposes
#The filtering filters samples/individuals to those that SNAPTRON has RNA seq samples for (ie. not all of v7, though a few more than v6p)
#Output files made (all we be written in filter_pre_process_output_dir)
##1. *_rnaseq_sample_ids.txt where * is the conventional_tissue_name. These files are simply a list of all GTEx RNASEQ samples used for this tissue.
##2. all_rna_samples_file_locations.txt. This file is just a list of all the absolute file locations from part 1. Used for downstream analysis
##3. all_individuals.txt. A list of all of the individuals across all tissues. (going to ignore this for now. Because we are limiting to individuals with RNA seq samples)
##4. tissue_sample_counts.txt. A table that contains information on how many gtex samples are in each tissue. Column 0 is the gtex tissue name and column 1 is the number of samples in that tissue.
#NOTE: By samples, I mean RNA-seq samples. By individuals, I mean actual people. So for gtex 1 individual may have multiple samples. But 1 sample only has 1 individual. 
#Runtime: < 1 min
#Spot checked
###################################################################################################################
if false; then
sh filter_pre_process.sh $pre_process_output_dir $snaptron_gtex_samples_file $filter_pre_process_output_dir
fi

#A file generated from filter_pre_process.sh.
#Each line in the file corresponds to a given tissue's absolute location of a file that contains all of the gtex sample ids for that tissue
all_rna_samples_file_locations=$filter_pre_process_output_dir"/all_rna_samples_file_locations.txt"


##################################################################################################################
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
###########################################
##Subpart_4: 'map_clusters_to_genes.py '
#This script takes that output from subpart_e (*_hg19_filtered_&_reads.txt) and maps clusters to hg19 protein coding genes
#ie, a cluster is assigned to a gene if one of the cluster's jxns falls within the gencode 'gene' position for that gene
#Ie, it is very possible for 1 cluster to be assigned to more than one gene
#Also, possible for a cluster to not be assigned to any gene (Then that cluster and its corresponding junctions) are removed from downstream
##output_saved as *_hg19_filtered_gene_mapped.txt where * is the tissue type. Saved in 'clusters_filter_output_dir'
##################################################################################################################
if false; then
tissue_type="Adipose_Subcutaneous_Analysis"
alternative_name="Adipose - Subcutaneous"
sbatch generate_junctions.sh $tissue_type $filter_pre_process_output_dir $snaptron_gtex_junction_file $snaptron_gtex_samples_file $leafcutter_code_dir $junctions_output_dir $clusters_output_dir $clusters_filter_output_dir $min_reads $liftover_directory $gencode_hg19_gene_annotation_file
fi

if false; then
while read tissue_type alternative_name; do
	sbatch generate_junctions.sh $tissue_type $filter_pre_process_output_dir $snaptron_gtex_junction_file $snaptron_gtex_samples_file $leafcutter_code_dir $junctions_output_dir $clusters_output_dir $clusters_filter_output_dir $min_reads $liftover_directory $gencode_hg19_gene_annotation_file
done<$tissue_list_input_file
fi

##################################################################################################################
##################################################################################################################




total_nodes="5"

if false; then
tissue_type="Adipose_Subcutaneous_Analysis"
node_number="4"
	covariate_file=$covariate_directory_v6$tissue_type".covariates.txt"
	tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_gene_mapped.txt"
	rna_seq_samples_file=$filter_pre_process_output_dir$tissue_type"_rnaseq_sample_ids.txt"
	outlier_calling_samples_file=$outlier_calling_samples_dir$tissue_type"_"$filter_global_outlier_method"_used_samples.txt"
		tissue_specific_outlier_file=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$filter_global_outlier_method"_"$node_number
		sh call_outlier_dm.sh $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $sample_attribute_file $tissue_list_input_file $covariate_file $rna_seq_samples_file $num_pc $filter_global_outlier_method $v6_sample_attribute_file $outlier_calling_samples_file
 fi

if false; then
while read tissue_type; do 

	covariate_file=$covariate_directory_v6$tissue_type".covariates.txt"
	tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg19_filtered_gene_mapped.txt"
	rna_seq_samples_file=$filter_pre_process_output_dir$tissue_type"_rnaseq_sample_ids.txt"
	outlier_calling_samples_file=$outlier_calling_samples_dir$tissue_type"_"$filter_global_outlier_method"_used_samples.txt"

	for node_number in $(seq 0 `expr $total_nodes - "1"`); do
		tissue_specific_outlier_file=$outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$filter_global_outlier_method"_"$node_number
		if false; then
		sbatch call_outlier_dm.sh $tissue_type $tissue_specific_jxn_file $tissue_specific_outlier_file $outlier_calling_dm_output_dir $max_dm_junctions $jxn_filter_method $node_number $total_nodes $covariate_regression_method $sample_attribute_file $tissue_list_input_file $covariate_file $rna_seq_samples_file $num_pc $filter_global_outlier_method $v6_sample_attribute_file $outlier_calling_samples_file
		fi
	done
	echo $tissue_type
	cat $outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$filter_global_outlier_method"_"*"_pvalue.txt" > $outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$filter_global_outlier_method"_pvalue.txt"

	cat $outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$filter_global_outlier_method"_"*"_md.txt" > $outlier_calling_dm_output_dir$tissue_type"_"$jxn_filter_method"_"$covariate_regression_method"_"$num_pc"_"$filter_global_outlier_method"_md.txt"


done<"/scratch1/battle-fs1/bstrober/rare_variants/rare_splice/outlier_calling/downloaded_data/gtex_tissues.txt"
fi



















##################################################################################################################
#OLD CODE!!! CURRENTLY IGNORE
##################################################################################################################

#output_dir for call_outlier_rare.sh
outlier_calling_rare_output_dir="/scratch1/battle-fs1/bstrober/rare_variants/rare_splice/outlier_calling/v2/outlier_calling_rare/"
#output_dir for call_outlier_rare_gtex.sh
outlier_calling_rare_gtex_output_dir="/scratch1/battle-fs1/bstrober/rare_variants/rare_splice/outlier_calling/v2/outlier_calling_rare_gtex/"



##################################################################################################################
#recount rareness calling
##################################################################################################################

if false; then
while read tissue_type alternative_name; do

tissue_specific_jxn_file=$clusters_filter_output_dir$tissue_type"_hg38_filtered_"$min_reads"_reads.txt"
tissue_specific_rare_outlier_file_root=$outlier_calling_rare_output_dir$tissue_type"_"
sbatch call_outlier_rare.sh $tissue_type $tissue_specific_jxn_file $tissue_specific_rare_outlier_file_root $snaptron_directory $liftover_directory $outlier_calling_rare_output_dir

done<$tissue_list_input_file
fi


##################################################################################################################
#gtex rareness calling
##################################################################################################################
if false; then
sbatch call_outlier_rare_gtex_shell.sh $tissue_list_input_file $clusters_filter_output_dir $min_reads $outlier_calling_rare_gtex_output_dir
fi



