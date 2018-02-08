import numpy as np
import os
import sys
import pdb


##Get list of tuples. Where tuple[0] is GTEx tissue name conventional string and tuple[1] is the string of the tissue name that appears
##in sample_attribute_file
def get_tissue_pairs(tissue_list_input_file):
    f = open(tissue_list_input_file)
    tissue_pairs = []
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        tissue_pairs.append((data[0], data[1]))
    return tissue_pairs


##Extract all GTEx RNAseq samples that belong to this tissue and write them to output_file
def extract_tissue_specific_samples(sample_attribute_file, conventional_tissue_name, attribute_file_tissue_name, indis, output_file):
    f = open(sample_attribute_file)
    t = open(output_file, 'w')
    header_count = 0
    used_tissue_indis = {}  # keep track of which (tissue,indi_ids) have been used to ensure we do not have more than one individual per tissue
    count = 0  # Keep track of how many samples are in each tissue
    sample_ids = {}
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if header_count == 0:  # skip the header
            header_count = header_count + 1
            continue
        tiss = data[13]
        sample_type = data[27]
        sample_id = data[0]  # ID corresponding to RNA seq sample
        info = data[0].split('-')   # Extracting individual id from sample id
        indi_id = info[0] + '-' + info[1]  # individual id
        tissue_indi = tiss + '_' + indi_id
        if sample_type != 'RNASEQ' or tiss != attribute_file_tissue_name:  # We only care about sample_type == RNASEQ and the correct tissue type
            continue
        if tissue_indi in used_tissue_indis:  # keep track of which (tissue,indi_ids) have been used to ensure we do not have more than one individual per tissue
            print('ERRRROOROROROR')
        indis[indi_id] = 1  # keep track of which individuals are used
        used_tissue_indis[tissue_indi] = 1
        count = count + 1
        if sample_id in sample_ids:  # make sure there are no repeats in sample_ids per tissue (or ever!)
            print('erororoor')
        sample_ids[sample_id] = 1
    #Write to output in alphabetical order
    for sample_id in sorted(sample_ids.keys()):
        t.write(sample_id + '\n')
    t.close()
    return indis, count


def write_all_individuals(output_file, indis):
    t = open(output_file, 'w')
    for indi in sorted(indis):
        t.write(indi + '\n')
    t.close()

sample_attribute_file = sys.argv[1]  # File containing which GTEx samples are to be used for each tissue. Filter list on val(column[27]) == 'RNASEQ'
tissue_list_input_file = sys.argv[2]  # List of gtex tissues. First colum is traditional GTEx tissue name. Second column is GTEx tissue name as listed in $sample_attribute_file
pre_process_output_dir = sys.argv[3]  # output_dir location


##################################################################################################################
#This script produces lists of what gtex samples will be used in which tissue
#Output files made (all we be written in pre_process_output_dir)
##1. *_rnaseq_sample_ids.txt where * is the conventional_tissue_name. These files are simply a list of all GTEx RNASEQ samples used for this tissue.
##2. all_rna_samples_file_locations.txt. This file is just a list of all the absolute file locations from part 1. Used for downstream analysis
##3. all_individuals.txt. A list of all of the individuals across all tissues.
##4. tissue_sample_counts.txt. A table that contains information on how many gtex samples are in each tissue. Column 0 is the gtex tissue name and column 1 is the number of samples in that tissue.

#NOTE: By samples, I mean RNA-seq samples. By individuals, I mean actual people. So for gtex 1 individual may have multiple samples. But 1 sample only has 1 individual. 
###################################################################################################################

tissue_pairs = get_tissue_pairs(tissue_list_input_file)
indis = {}  # Keep track of all individuals we are going to test
t1 = open(pre_process_output_dir + 'tissue_sample_counts.txt', 'w')
t2 = open(pre_process_output_dir + 'all_rna_samples_file_locations.txt', 'w')
#Loop through tisssues
for tissue_pair in tissue_pairs:
    conventional_tissue_name = tissue_pair[0]
    attribute_file_tissue_name = tissue_pair[1]  # gtex id format used in sample_attribute_file
    tissue_specific_sample_file = pre_process_output_dir + conventional_tissue_name + '_rnaseq_sample_ids.txt'  # tissue_specific file that contains the sample ids in that tissue
    #main script to extract the samples in this tissue
    indis, tissue_specific_sample_count = extract_tissue_specific_samples(sample_attribute_file, conventional_tissue_name, attribute_file_tissue_name, indis, tissue_specific_sample_file)
    t1.write(conventional_tissue_name + '\t' + str(tissue_specific_sample_count) + '\n')
    t2.write(tissue_specific_sample_file + '\n')
t2.close()
t1.close()
#Create a file that contains all individuals across all of the tissues.
all_indi_output_file = pre_process_output_dir + 'all_individuals.txt'
write_all_individuals(all_indi_output_file, indis)
