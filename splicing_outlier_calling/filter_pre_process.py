import numpy as np
import os
import sys
import pdb
import gzip


#Extract list dictionary of all gtex samples that are present in snaptron
def get_samples_in_snaptron(file_name):
    count = 0  # for header
    hits = {}
    f = open(file_name)
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if count == 0:  # skip header
            count = count + 1
            continue
        gtex_sample_id = data[7]
        hits[gtex_sample_id] = 1
    return hits


def get_tissues(file_name):
    f = open(file_name)
    arr = []
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        arr.append(data[0])
    return arr

pre_process_output_dir = sys.argv[1]
snaptron_gtex_samples_file = sys.argv[2]
filter_pre_process_output_dir = sys.argv[3]
tissue_list_input_file = sys.argv[4]


tissues = get_tissues(tissue_list_input_file)


snaptron_samples = get_samples_in_snaptron(snaptron_gtex_samples_file)

t_tiss_counts = open(filter_pre_process_output_dir + 'tissue_sample_counts.txt', 'w')  # File handle

indiz = {}  # Keep track of all gtex individuals analyzed

for tissue in tissues:
    input_file_name = pre_process_output_dir + tissue + '_rnaseq_sample_ids.txt'  # Sample id file before removing samples not in snaptron
    output_file_name = filter_pre_process_output_dir + tissue + '_rnaseq_sample_ids.txt'  # Sample id file after removing samples not in snaptron
    f = open(input_file_name)
    t = open(output_file_name, 'w')

    tissue_count = 0  # Keep track of how many samples are in each tissue (after removing samples not in snaptron)
    for line in f:
        sample_id = line.strip()
        if sample_id not in snaptron_samples:  # Only include samples that are present in snaptron
            print('Sample not contained in SNAPTRON ' + tissue + ' ' + sample_id + '.. (this is OK for now)')
            continue
        indi_id = sample_id.split('-')[0] + '-' + sample_id.split('-')[1]  # Convert sample id -> individual id
        indiz[indi_id] = 1
        tissue_count = tissue_count + 1
        t.write(sample_id + '\n')
    t_tiss_counts.write(tissue + '\t' + str(tissue_count) + '\n')
t_tiss_counts.close()


t = open(filter_pre_process_output_dir + 'all_individuals.txt', 'w')
for indi in sorted(indiz):
    t.write(indi + '\n')
t.close()
