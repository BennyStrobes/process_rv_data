import numpy as np
import os
import sys
import pdb
import gzip

def get_all_v7_samples(all_rna_samples_file_locations):
	gtex_sample_ids = {}
	f = open(all_rna_samples_file_locations)
	for i,line in enumerate(f):
		tissue_specific_file_name = line.rstrip()
		g = open(tissue_specific_file_name)
		for line2 in g:
			line2 = line2.rstrip()
			gtex_sample_ids[line2] = 1
	return gtex_sample_ids

def get_samples_in_snaptron(file_name,v7_samples):
	count = 0 #for header
	hits = {}
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if count == 0: #skip header
			count = count + 1
			continue
		rail_sample_id = data[0] 
		gtex_sample_id = data[7]
		if gtex_sample_id not in v7_samples:
			continue
		hits[gtex_sample_id] = 1
	return hits
def print_helper_indi(hits,output_file):
	indi = {}
	t = open(output_file,'w')
	for key in hits.keys():
		info = key.split('-')
		indi[info[0] + '-' + info[1]] =1
	for ind in sorted(indi.keys()):
		t.write(ind + '\n')
	t.close()

def print_helper2(snaptron_samples,all_rna_samples_file_locations,filter_pre_process_output_dir):
	t = open(filter_pre_process_output_dir + 'all_rna_samples_file_locations.txt','w')
	t1 = open(filter_pre_process_output_dir + 'tissue_sample_counts.txt','w')
	f = open(all_rna_samples_file_locations)
	for line in f:
		line = line.rstrip()
		core_string = line.split('/')[-1].split('_rnaseq')[0]
		new_output_file = filter_pre_process_output_dir + core_string + '_rnaseq_sample_ids.txt'
		t.write(new_output_file + '\n')
		t2 = open(new_output_file,'w')
		g = open(line)
		count = 0
		for liner in g:
			liner = liner.rstrip()
			indi = liner.split('-')[0] + '-' + liner.split('-')[1]
			if liner in snaptron_samples:
				count = count + 1
				t2.write(liner + '\n') 
		g.close()
		t2.close()
		t1.write(core_string + '\t' + str(count) + '\n')
	t1.close()
	t.close()
	f.close()


def get_tissues(file_name):
	f = open(file_name)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		arr.append(data[0])
	return arr


def get_individuals_used_in_dm_splicing(tissues,dm_outlier_samples_directory):
	dicti = {}
	dicti['all_tissues'] = {}
	for tissue in tissues:
		input_file = dm_outlier_samples_directory + tissue + '_filter2_used_samples.txt'
		dicti[tissue] = {}
		f = open(input_file)
		for line in f:
			indi = line.rstrip()
			dicti[tissue][indi] = 1
			dicti['all_tissues'][indi] = 1
		f.close()
	return dicti

pre_process_output_dir= sys.argv[1]
snaptron_gtex_samples_file= sys.argv[2]
filter_pre_process_output_dir = sys.argv[3]
tissue_list_input_file= sys.argv[4]


tissues = get_tissues(tissue_list_input_file)



all_rna_samples_file_locations = pre_process_output_dir + 'all_rna_samples_file_locations.txt'
v7_samples = get_all_v7_samples(all_rna_samples_file_locations)
snaptron_samples = get_samples_in_snaptron(snaptron_gtex_samples_file,v7_samples)
print_helper_indi(snaptron_samples,filter_pre_process_output_dir + 'all_individuals.txt')
print_helper2(snaptron_samples,all_rna_samples_file_locations,filter_pre_process_output_dir)
