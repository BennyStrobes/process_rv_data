import sys
import pdb
import os
import numpy as np 
import time
import gzip
from scipy import stats
import scipy.stats as ss

#Extract dictionary (junction_counts) of all junctions that we wish to know how 'rare' they are. Keys will be junction ids and values will be counts of number of samples that express that junction
#Converts from leafcutter coordinates to snaptron coordinates
def extract_junctions(tissue_specific_jxn_file):
	f = open(tissue_specific_jxn_file)
	junction_counts = {}
	count =0
	for line in f:
		line = line.rstrip()
		data = line.split()
		#ignore header
		if count == 0:
			count = count + 1
			continue
		jxn_id = data[0]
		info = jxn_id.split(':')
		#convert from leafcutter coordinates to snaptron coordinates
		new_jxn_id = info[0] + ':' + str(int(info[1]) + 1) + ':' + str(int(info[2]) -1)
		#Add jxn to dictionary
		junction_counts[new_jxn_id] = 0
	return junction_counts

#Update number of samples that express this junction (junction_counts) based on this study
def update_junction_counts(snaptron_jxn_file,junction_counts,study):
	f = gzip.open(snaptron_jxn_file)
	#loop through junctions
	for line in f:
		line = line.rstrip()
		data = line.split()
		jxn_id = data[1] + ':' + data[2] + ':' + data[3]
		#Each sample is seperated by comma
		samples_with_jxn = data[11].split(',')[1:]
		#Gtex study
		if study == 'gtex':
			#There is a bug in Chris's code so we simply have to ignore samples with 0 reads mapped (not an issue besides for this)
			#Also, only consider RAIL ids in the file_handler object
			nonzero_samples = 0
			for sample in samples_with_jxn:
				#sample_id and the read counts for thaat sample are seperated by ':'
				info = sample.split(':')
				sample_id = info[0]
				sample_counts = info[1]
				if int(sample_counts) > 0:
					nonzero_samples = nonzero_samples + 1
		#everything else besides for gtex
		else:
			nonzero_samples = len(samples_with_jxn)
		if jxn_id in junction_counts:
			junction_counts[jxn_id] = junction_counts[jxn_id] + nonzero_samples
	f.close()
	return junction_counts


#main driver script for analysis
def start_rareness_analysis(tissue_type,tissue_specific_jxn_file,hg38_jxn_rareness_file,snaptron_directory,studies_to_use):
	#Extract dictionary (junction_counts) of all junctions that we wish to know how 'rare' they are. Keys will be junction ids and values will be counts of number of samples that express that junction
	#converts from leafcutter coordinates to snaptron coordinates
	junction_counts = extract_junctions(tissue_specific_jxn_file)

	#Loop through all specified studies (tcga,sra,gtex,etc)
	for study in studies_to_use:
		#file used as input to snaptron
		snaptron_jxn_file =snaptron_directory + study + '_hg38_jxns.bgz'
		#Update number of samples that express this junction (junction_counts) based on this study
		junction_counts = update_junction_counts(snaptron_jxn_file,junction_counts,study)

		#error checking
		for jxn in junction_counts.keys():
			if junction_counts[jxn] == 0:
				print('errororo')
				print(jxn)

	#print output
	t = open(hg38_jxn_rareness_file,'w')
	for jxn in junction_counts.keys():
		t.write(jxn + '\t' + str(junction_counts[jxn]) + '\n')
	t.close()


tissue_type = sys.argv[1]
#Input file
tissue_specific_jxn_file = sys.argv[2]
#output_file
hg38_jxn_rareness_file = sys.argv[3]
#Directory containing files that contain information on which samples are present for each junction
snaptron_directory = sys.argv[4]

#main driver script
start_rareness_analysis(tissue_type,tissue_specific_jxn_file,hg38_jxn_rareness_file,snaptron_directory,['gtex','srav2'])