####Ben Strober
####3/1/17
#########These scripts will:
##############1. Extract rare variants from vcf files
##############2. Filter those variants to leave only those variants within certain areas around the gene (TBD)

#############################################################
#Input data
#############################################################
##This is the GTEx v7 VCF file
vcf_file="/work-zfs/abattle4/lab_data/GTEx_v7/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz"

#Use this file to filter on only european indiduals.
#Filter is applied by
subject_phenotype_file="/work-zfs/abattle4/lab_data/GTEx_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SubjectPhenotypesDS.txt"

#As Snaptron does not have all of the rnaseq samples for v7, we have to limit to those in snaptron (for now)
#Those individuals found in snaptron can be found here
rna_seq_individuals="/work-zfs/abattle4/bstrober/rare_variant/rare_splice/process_rv_data/splicing_outlier_calling/filter_pre_process/all_individuals.txt"

#Gencode hg19 gene annotation file
gencode_hg19_gene_annotation_file="/work-zfs/abattle4/lab_data/hg19/gencode_gene_annotations/gencode.v19.annotation.gtf.gz"

# Directory containing 1K genomes data (downloaded from ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/)
k_genomes_input_dir="/work-zfs/abattle4/lab_data/1k_genomes/"



##############################################################
#Used directories
###############################################################
variant_calling_root_dir="/work-zfs/abattle4/bstrober/rare_variant/rare_splice/process_rv_data/variant_calling/"

#There is a file for each individual (*_SNPS.bed where * is individual id), that contains all of their variants after the filtering pipeline (by emily and joe)
variant_bed_dir=$variant_calling_root_dir"variant_beds/"

#Contains output of variant_to_gene_mapping.sh
variant_gene_mapping_output_dir=$variant_calling_root_dir"variant_gene_mapping/"

# Contains names of all samples that we have WGS for
wgs_samples_dir=$variant_calling_root_dir"wgs_samples/"

# Contains processed 1KG genomes data
k_genomes_process_dir=$variant_calling_root_dir"1k_genomes_processed/"

##############################################################
#Parameters
###############################################################
#gene_mapping_method refers to the way we map variants to genes.
#method_1: Limit to protein coding. A gene considered to be what gencode defines as 'gene', as well as 1 KB upstream and 1KB downstream
gene_mapping_method="method_1"


#################################################################################################################
##################################################################################################################
#variant_extraction_nick_emily.sh extracts variants from vcf after applying some filters to it (1000 genomes,etc)
#It will produce:
##1.$input_individuals: These are all of the individuals we tried to extract WGS from. It is a cross between filtering on only rna_seq_individuals and european ancenstry individuals
##2. $tested_individuals : These are all of the individuals we actually extracted WGS from (we lost ~60)
##3.$variant_bed_dir* : where * is the individual id. This is a bed file of all the variants produced.
#NOTE: This code was produced by Joe Davis and Emily Tsang
#It relies on set directories in ~/.bash_rc. So in future, if you want to change anything, need to change ~/.bashrc
#See https://github.com/joed3/GTExV6PRareVariation for reference on their code
total_jobs="80"

if false; then
sbatch variant_extraction_nick_emily_serial_part.sh $vcf_file $subject_phenotype_file $rna_seq_individuals $wgs_samples_dir
fi

job_number="0"
if false; then
sbatch variant_extraction_nick_emily_parallel_part.sh $vcf_file $wgs_samples_dir $job_number $total_jobs
fi


if false; then
job_number="0"
sbatch variant_extraction_nick_emily_parallel_part.sh $vcf_file $wgs_samples_dir $job_number $total_jobs
fi

if false; then
for job_number in $(seq 1 `expr $total_jobs - "1"`); do
	echo $job_number
	sbatch variant_extraction_nick_emily_parallel_part.sh $vcf_file $wgs_samples_dir $job_number $total_jobs
done
fi



#################################################
# Process 1K-genomes data
##################################################
if false; then
sbatch extract.1kg.AF.sh $k_genomes_input_dir $k_genomes_process_dir $wgs_samples_dir
fi

#File created all individuals used
#Basically a result of cross filtering the subject_phenotype_file file and rna_seq_individuals file
input_individuals=$wgs_samples_dir"EA_gtex_VCFids.txt"
#################################################################################################################
##################################################################################################################
#variant_to_gene_mapping.sh maps variants (discovered in variant_extraction_nick_emily.sh) to genes (if they lie on the gene)
#It produces the output file "all_variants_gene_mapping_"$gene_mapping_method".txt" (in $variant_gene_mapping_output_dir)
#$gene_mapping_method refers to the way by which a variant is matpped to a given gene. So far, gene_mapping_method's that are currently implemented:
###method_1: Limit to protein coding. A gene considered to be what gencode defines as 'gene', as well as 1 KB upstream and 1KB downstream
all_variants_gene_mapping_file=$variant_gene_mapping_output_dir"all_variants_gene_mapping_"$gene_mapping_method".txt"
sbatch variant_to_gene_mapping.sh $wgs_samples_dir"individuals/" $input_individuals $all_variants_gene_mapping_file $gencode_hg19_gene_annotation_file $gene_mapping_method


