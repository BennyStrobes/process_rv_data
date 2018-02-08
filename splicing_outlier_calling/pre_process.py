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

#Extract individuals (first line) of the covariate file in the form of a dictionary
def get_eqtl_indi_ids_from_covariate_file(file_name):
    ids = {}
    f = open(file_name)
    for line in f:
        line = line.rstrip()
        data = line.split()
        for ele in data[1:]:
            ids[ele] = 1
        break
    return ids


##Extract all GTEx RNAseq samples that belong to this tissue and write them to output_file
def extract_tissue_specific_samples(sample_attribute_file, conventional_tissue_name, attribute_file_tissue_name, indis, output_file, covariate_directory_v6):
    f = open(sample_attribute_file)
    t = open(output_file, 'w')
    # Get list of gtex individual ids used in eqtl analysis. Extract this from the covariate files used in the eqtl analysis
    eqtl_indi_ids = get_eqtl_indi_ids_from_covariate_file(covariate_directory_v6 + conventional_tissue_name + '.covariates.txt')
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
        if len(data) < 31:  # Skip Short lines (they are not RNA seq)
            continue
        tiss = data[13]
        sample_id = data[0]  # ID corresponding to RNA seq sample
        info = data[0].split('-')  # Extracting individual id from sample id
        indi_id = info[0] + '-' + info[1]  # individual id
        sample_annotation = data[28]  # 'USE ME' or 'FLAGGED'
        tissue_indi = tiss + '_' + indi_id

        if indi_id not in eqtl_indi_ids:  # Ignore samples not used in eqtl analysis
            continue
        if tiss != attribute_file_tissue_name:  # We only care about samples with the correct tissue type
            continue
        if sample_annotation != 'USE ME':  # Ignore flagged samples
            continue
        if tissue_indi in used_tissue_indis:  # keep track of which (tissue,indi_ids) have been used to ensure we do not have more than one individual per tissue
            print('ERRRROOROROROR')
            pdb.set_trace()

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

    if len(sample_ids) != len(eqtl_indi_ids):  # Double check to ensure that the sample ids we have extracted are the same length as individual ids used in eqtl analysis
        print('EROROROROROR')
        pdb.set_trace()
    return indis, count


#Create a file that contains all individuals across all of the tissues.
def write_all_individuals(output_file, indis):
    t = open(output_file, 'w')
    for indi in sorted(indis):
        t.write(indi + '\n')
    t.close()

sample_attribute_file = sys.argv[1]  # File containing which GTEx samples are to be used for each tissue. Filter list on val(column[28]) == 'USE ME'. Also filter if length(line) < 30 (not full/complete)
tissue_list_input_file = sys.argv[2]  # List of gtex tissues. First colum is traditional GTEx tissue name. Second column is GTEx tissue name as listed in $sample_attribute_file
pre_process_output_dir = sys.argv[3]  # output_dir location
covariate_directory_v6 = sys.argv[4]  # Used to filter tissue specific samples.
##################################################################################################################
#This script produces lists of what gtex samples will be used in which tissue
#Output files made (all we be written in pre_process_output_dir)
##1. *_rnaseq_sample_ids.txt where * is the conventional_tissue_name. These files are simply a list of all GTEx RNASEQ samples used for this tissue.
##2. all_individuals.txt. A list of all of the individuals across all tissues.
##3. tissue_sample_counts.txt. A table that contains information on how many gtex samples are in each tissue. Column 0 is the gtex tissue name and column 1 is the number of samples in that tissue.

#NOTE: By samples, I mean RNA-seq samples. By individuals, I mean actual people. So for gtex 1 individual may have multiple samples. But 1 sample only has 1 individual.
###################################################################################################################

tissue_pairs = get_tissue_pairs(tissue_list_input_file)  # Extract array of tuples. Where each tuple is the two ways to spell the tissue type
indis = {}  # Keep track of all individuals we are going to test
t1 = open(pre_process_output_dir + 'tissue_sample_counts.txt', 'w')
#Loop through tisssues
for tissue_pair in tissue_pairs:
    conventional_tissue_name = tissue_pair[0]
    attribute_file_tissue_name = tissue_pair[1]  # gtex id format used in sample_attribute_file
    tissue_specific_sample_file = pre_process_output_dir + conventional_tissue_name + '_rnaseq_sample_ids.txt'  # tissue_specific file that contains the sample ids in that tissue
    #main script to extract the samples in this tissue
    indis, tissue_specific_sample_count = extract_tissue_specific_samples(sample_attribute_file, conventional_tissue_name, attribute_file_tissue_name, indis, tissue_specific_sample_file, covariate_directory_v6)
    t1.write(conventional_tissue_name + '\t' + str(tissue_specific_sample_count) + '\n')
t1.close()
#Create a file that contains all individuals across all of the tissues.
all_indi_output_file = pre_process_output_dir + 'all_individuals.txt'
write_all_individuals(all_indi_output_file, indis)
