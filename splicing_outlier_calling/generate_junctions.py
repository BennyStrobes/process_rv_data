import numpy as np
import os
import sys
import pdb
import gzip
'''

def extract_all_gtex_sample_ids(all_rna_samples_file_locations,tissue_num):
	gtex_sample_ids = {}
	gtex_sample_ids_tiss_spec = {}
	f = open(all_rna_samples_file_locations)
	for i,line in enumerate(f):
		tissue_specific_file_name = line.rstrip()
		g = open(tissue_specific_file_name)
		for line2 in g:
			line2 = line2.rstrip()
			gtex_sample_ids[line2] = 1
			if i == tissue_num:
				gtex_sample_ids_tiss_spec[line2] =1
	return gtex_sample_ids,gtex_sample_ids_tiss_spec


def extract_sample_ids(all_rna_samples_file_locations,snaptron_gtex_samples_file,tissue_num):
	rail_sample_ids = {}
	gtex_to_rail_sample_id_converter = {}
	gtex_sample_ids,gtex_sample_ids_tiss_spec = extract_all_gtex_sample_ids(all_rna_samples_file_locations,tissue_num)
	f = open(snaptron_gtex_samples_file)
	count = 0 #for header
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if count == 0: #skip header
			count = count + 1
			continue
		rail_sample_id = data[0]
		gtex_sample_id = data[7]
		if gtex_sample_id not in gtex_sample_ids_tiss_spec:
			continue
		rail_sample_ids[rail_sample_id] = []
		gtex_to_rail_sample_id_converter[gtex_sample_id] = rail_sample_id
	return rail_sample_ids,gtex_to_rail_sample_id_converter,gtex_sample_ids

def print_helper_rna_seq_individuals(gtex_to_rail_sample_id_converter,output_file):
	t =open(output_file,'w')
	individuals = {}
	for key in gtex_to_rail_sample_id_converter.keys():
		info = key.split('-')
		individuals[info[0] + '-' + info[1]] =1
	for individual in sorted(individuals.keys()):
		t.write(individual + '\n')
	t.close()

def extract_junctions(snaptron_gtex_junction_file,rail_sample_ids):
	f = gzip.open(snaptron_gtex_junction_file)
	count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		samples_with_jxn = data[11].split(',')[1:]
		for sample in samples_with_jxn:
			info = sample.split(':')
			sample_id = info[0]
			sample_counts = info[1]
			if sample_counts != '0' and sample_id in rail_sample_ids:
				jxn_file_string = data[1] + '\t' + data[2] + '\t' + data[3] + '\t.\t' + sample_counts + '\t' + data[5] + '\n'
				rail_sample_ids[sample_id].append(jxn_file_string)
		count = count + 1
		if count %100000 == 0:
			print(count)
	return rail_sample_ids
		
def check_no_missing_samples(rail_sample_ids):
	missing = 0
	for rail_sample_id in sorted(rail_sample_ids.keys()):
		arr = rail_sample_ids[rail_sample_id]
		if len(arr) == 0:
			missing = missing + 1
	if missing > 0:
		print('ERROR: There are ' + str(missing) + ' rna seq samples with no junctions. This seems wrong')


def print_helper_samples_files(all_rna_samples_file_locations,rail_sample_ids,gtex_to_rail_sample_id_converter,junctions_output_dir):
	t = open(junctions_output_dir + 'all_used_rna_samples_file_locations.txt','w')
	f = open(all_rna_samples_file_locations)
	for file_name in f:
		file_name = file_name.rstrip()
		g = open(file_name)
		new_rnaseq_sample_file = junctions_output_dir + file_name.split('_rnaseq')[0].split('/')[-1] + '_used_rnaseq_sample_ids.txt'
		t.write(new_rnaseq_sample_file + '\n')
		t1 = open(new_rnaseq_sample_file,'w')
		for line in g:
			line = line.rstrip()
			gtex_sample_id = line
			if gtex_sample_id in gtex_to_rail_sample_id_converter:
				t1.write(gtex_sample_id + '\n')
		t1.close()
	t.close()

def print_helper_junctions_files(all_rna_samples_file_locations,rail_sample_ids,gtex_to_rail_sample_id_converter,junctions_output_dir):
	f = open(all_rna_samples_file_locations)
	tissue_level_file_locations_file = junctions_output_dir + 'tissue_jxn_file_locations.txt'
	t = open(tissue_level_file_locations_file,'w')
	for file_name in f:
		file_name = file_name.rstrip()
		g = open(file_name)
		specific_tissue_junction_file_locations = junctions_output_dir + file_name.split('_rnaseq')[0].split('/')[-1] + '_jxn_file_locations.txt'
		t.write(specific_tissue_junction_file_locations + '\n')
		t1 = open(specific_tissue_junction_file_locations,'w')
		for line in g:
			line = line.rstrip()
			gtex_sample_id = line
			if gtex_sample_id in gtex_to_rail_sample_id_converter:
				sample_specific_jxn_file = junctions_output_dir + gtex_sample_id + '.junc'
				t1.write(sample_specific_jxn_file + '\n')
				t2 = open(sample_specific_jxn_file,'w')
				rail_sample_id = gtex_to_rail_sample_id_converter[gtex_sample_id]
				string_array = rail_sample_ids[rail_sample_id]
				for string in string_array:
					t2.write(string)
				t2.close()
		t1.close()
	t.close()


#Main driver script for this function
def run_version_0(tissue_list_input_file,all_rna_samples_file_locations,snaptron_gtex_junction_file,snaptron_gtex_samples_file,junctions_output_dir,tissue_num):
	#extract all rail ids that correspond to gtex sample ids in any of the gtex tissues
	#Also extract a converter that changes from gtex sample ids to rail ids
	rail_sample_ids,gtex_to_rail_sample_id_converter,gtex_sample_ids = extract_sample_ids(all_rna_samples_file_locations,snaptron_gtex_samples_file,tissue_num)
	print_helper_rna_seq_individuals(gtex_sample_ids,junctions_output_dir + 'rna_seq_individuals.txt')
	rail_sample_ids = extract_junctions(snaptron_gtex_junction_file,rail_sample_ids)
	check_no_missing_samples(rail_sample_ids)
	#print_helper_samples_files(all_rna_samples_file_locations,rail_sample_ids,gtex_to_rail_sample_id_converter,junctions_output_dir)
	#print_helper_junctions_files(all_rna_samples_file_locations,rail_sample_ids,gtex_to_rail_sample_id_converter,junctions_output_dir)
'''

#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################
#########################################################################################################################################################################

#Get gtex sample ids
def get_gtex_sample_ids(tissue_specific_gtex_samples_file):
	f = open(tissue_specific_gtex_samples_file)
	sample_ids = {}
	for line in f:
		line = line.rstrip()
		sample_ids[line] = 1
	return sample_ids

#create dictionary that converts from gtex sample id to Rail (snaptron) sample id
def get_gtex_to_rail_sample_id_dictionary(gtex_sample_ids,snaptron_gtex_samples_file):
	f = open(snaptron_gtex_samples_file)
	count = 0 #for header
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if count == 0: #skip header
			count = count + 1
			continue
		rail_sample_id = data[0]
		gtex_sample_id = data[7]
		if gtex_sample_id not in gtex_sample_ids: #Ignore those lines that are not gtex samples for this specific tissue
			continue
		dicti[gtex_sample_id] = rail_sample_id
	if len(dicti) != len(gtex_sample_ids): #They must be equal or else snaptron_gtex_samples_file is missing some lines
		print('erorororo')
		pdb.set_trace()
	return dicti


#First extract dictionary of gtex sample ids from tissue_specific_gtex_samples_file
#Then create a dictionary that converts from gtex sample id to Rail sample id (Used by snaptron)
def extract_sample_ids(tissue_specific_gtex_samples_file,snaptron_gtex_samples_file):
	gtex_sample_ids = get_gtex_sample_ids(tissue_specific_gtex_samples_file) #extract gtex_sample ids
	#create dictionary that converts from gtex sample id to Rail (snaptron) sample id
	gtex_to_rail_dictionary = get_gtex_to_rail_sample_id_dictionary(gtex_sample_ids,snaptron_gtex_samples_file)
	return gtex_sample_ids,gtex_to_rail_dictionary

#In order to stream the junction file and not save all of the data (too much memory), I created this data structure (file_handler)
#file_handler is a dictionary with keys that are the rail sample ids for this tissue. And the values of the keys correspond to file handles for each samle
def create_file_handle_object(gtex_sample_ids,gtex_to_rail_dictionary,junctions_output_dir,tissue_type):
	file_handler = {} #Core dictionary
	#Create output file that contains the absolute location of all junction sample files for this tissue
	t = open(junctions_output_dir + tissue_type + '_junction_file_locations.txt','w')
	#Loop through gtex_sample_ids for this tissue
	for gtex_sample_id in sorted(gtex_sample_ids.keys()):
		rail_sample_id = gtex_to_rail_dictionary[gtex_sample_id] #gtex sample id -> rail sample id
		output_file = junctions_output_dir + gtex_sample_id + '.junc' #absolute location of a sample's junction file
		file_handler[rail_sample_id] = open(output_file,'w') #add file handle to dictionary for this sample
		t.write(output_file + '\n')
	t.close()
	return file_handler

#Stream snaptron junction file and write junctions to files using file_handler
#Each sample has its own file. Files are in the leafcutter step 1 output format (https://github.com/davidaknowles/leafcutter)
def reformat_junction_data(snaptron_gtex_junction_file,file_handler):
	f = gzip.open(snaptron_gtex_junction_file)
	#loop through junctions
	for line in f:
		line = line.rstrip()
		data = line.split()
		#Each sample is seperated by comma
		samples_with_jxn = data[11].split(',')[1:]
		for sample in samples_with_jxn:
			#sample_id and the read counts for thaat sample are seperated by ':'
			info = sample.split(':')
			sample_id = info[0]
			sample_counts = info[1]
			#There is a bug in Chris's code so we simply have to ignore samples with 0 reads mapped (not an issue besides for this)
			#Also, only consider RAIL ids in the file_handler object
			if sample_counts != '0' and sample_id in file_handler:
				#format at leafcutter step 1 output
				jxn_file_formatted_string = data[1] + '\t' + data[2] + '\t' + data[3] + '\t.\t' + sample_counts + '\t' + data[5] + '\n'
				#write to file
				file_handler[sample_id].write(jxn_file_formatted_string)
	f.close()
	#close handles
	for sample_id in file_handler.keys():
		file_handler[sample_id].close()


#main driver script
def run_version_1(tissue_type,tissue_specific_gtex_samples_file,snaptron_gtex_junction_file,snaptron_gtex_samples_file,junctions_output_dir):
	#First extract dictionary of gtex sample ids from tissue_specific_gtex_samples_file
	#Then create a dictionary that converts from gtex sample id to Rail sample id (Used by snaptron)
	gtex_sample_ids,gtex_to_rail_dictionary = extract_sample_ids(tissue_specific_gtex_samples_file,snaptron_gtex_samples_file)
	#In order to stream the junction file and not save all of the data (too much memory), I created this data structure (file_handler)
	#file_handler is a dictionary with keys that are the rail sample ids for this tissue. And the values of the keys correspond to file handles for each samle
	file_handler = create_file_handle_object(gtex_sample_ids,gtex_to_rail_dictionary,junctions_output_dir,tissue_type)
	#Stream snaptron junction file and write junctions to files using file_handler
	#Each sample has its own file. Files are in the leafcutter step 1 output format (https://github.com/davidaknowles/leafcutter)
	reformat_junction_data(snaptron_gtex_junction_file,file_handler)




##Subpart_1 of 'generate_junctions.sh': 'generate_junctions.py'
####This part produces a junction file for each gtex rna-seq sample
####Ouput files made (all will be written to #junction_output_dir):
#######1. *.junc where * is a GTEx sample ID. This file is the exact format as the output of step 1 of leafcutter
#######2. *_junction_file_locations.txt where * is a tissue type. Each line is the absolute location of all of the &.junc files in this tissue. Where & is a GTEx sample id.


tissue_type = sys.argv[1]
pre_process_dir = sys.argv[2] #input directory
#File from downloaded from snaptron that contains all GTEx junctions
snaptron_gtex_junction_file = sys.argv[3] 
#File downloaded from snaptron that contains all GTEx samples in SNAPTRON
snaptron_gtex_samples_file = sys.argv[4]
#Ouptput directory
junctions_output_dir = sys.argv[5]

#Input file where each line is a GTEx rna sample id for this tissue.
tissue_specific_gtex_samples_file = pre_process_dir + tissue_type + '_rnaseq_sample_ids.txt'

#main driver script
run_version_1(tissue_type,tissue_specific_gtex_samples_file,snaptron_gtex_junction_file,snaptron_gtex_samples_file,junctions_output_dir)







