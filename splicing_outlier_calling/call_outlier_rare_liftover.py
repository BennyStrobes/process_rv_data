import sys
import pdb
import os
import numpy as np 
import fit_mvpolya as dm 
import time
import gzip
from scipy import stats
import scipy.stats as ss



#First make tempory bed file of our existing hg38 cluster positions. Will be used as input to liftover
def make_input_bed_file(cluster_file_hg38,temporary_hg38_bed_file):
	t = open(temporary_hg38_bed_file,'w')
	f = open(cluster_file_hg38)
	count = 0 #header
	for line in f:
		line = line.rstrip()
		data = line.split()

		jxn_id = data[0]
		chromer = jxn_id.split(':')[0]
		start_pos = jxn_id.split(':')[1]
		end_pos = jxn_id.split(':')[2]
		t.write(chromer + '\t' + start_pos + '\t' + end_pos + '\n')
	t.close()

#Run liftOver with parameters specified by Chris Wilks (https://github.com/ChristopherWilks/snaptron/blob/master/scripts/liftover_snaptron_coords.pl)
def run_liftover(input_file,output_file,missing_file,liftover_directory):
	stringer = liftover_directory + 'liftOver -minMatch=1.00 -ends=2 ' + input_file + ' ' + liftover_directory + 'hg38ToHg19.over.chain.gz ' + output_file + ' ' + missing_file
	os.system(stringer)

#Some jxns were not able to be mapped with liftover (lacked confidence). So first extract those unmapped jxns
def get_unmapped_jxns(temporary_missing_file):
	f = open(temporary_missing_file)
	unmapped_jxns = {}
	for line in f:
		if line.startswith('#'):
			continue
		line = line.rstrip()
		data = line.split()
		jxn_id = data[0] + ':' + data[1] + ':' + data[2]
		unmapped_jxns[jxn_id] = 1
	return unmapped_jxns

#If raw_leafcutter_file has n+1 lines (therefor n jxns), compute binary vector of length n called pass_filter
#If pass_filter == 1 then that jxn has a mapping in hg19:
def determine_jxns_that_pass_unmappable_filter(cluster_file_hg38,unmapped_jxns):
	pass_filter = []
	count = 0 #for header
	total_jxns = 0 #keep track of number of jxns in the file
	unmapped_jxn_count = 0 #keep track of unmapped jxns
	cluster_counts = {} #keep track of number of times a cluster is called
	valid_clusters = {} #only keep clusters that have more than one jxn (after unmapping filter)
	f = open(cluster_file_hg38)
	for line in f:
		line = line.rstrip()
		data = line.split()
		total_jxns = total_jxns + 1
		jxn_id_info = data[0].split(':')
		cluster_id = jxn_id_info[-1]
		jxn_identifier = jxn_id_info[0] + ':' + jxn_id_info[1] + ':' + jxn_id_info[2]
		if jxn_identifier in unmapped_jxns: #jxn was unmapped
			pass_filter.append(0)
			unmapped_jxn_count = unmapped_jxn_count + 1
		else: #jxn mapped
			pass_filter.append(1)
			if cluster_id not in cluster_counts:
				cluster_counts[cluster_id] =1
			else:
				cluster_counts[cluster_id] = cluster_counts[cluster_id] + 1
	if len(pass_filter) != total_jxns:  #simple check
		print('Length error')
	if unmapped_jxn_count != len(unmapped_jxns): #simple check
		print('lerngther errror')
	for cluster_id in cluster_counts.keys():
		if cluster_counts[cluster_id] > 1: #more than 1 jxn are in this cluster after filter
			valid_clusters[cluster_id] = 1
	return pass_filter,valid_clusters

#Extract converted hg19 coordinates
def get_hg19_coordinates(temporary_hg19_bed_file):

	hg19_coordinates=[]
	f = open(temporary_hg19_bed_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		cluster_id = data[0] + ':' + data[1] + ':' + data[2]
		hg19_coordinates.append(cluster_id)
	return hg19_coordinates

#print these new coordinates while filtering on those unmappable ones/ ones with clusters with only one jxn
def print_helper_liftover(cluster_file_hg38,cluster_file_hg19,pass_filter,valid_clusters,hg19_coordinates):
	f = open(cluster_file_hg38)
	t = open(cluster_file_hg19,'w')
	count = 0 #header
	jxn_count = 0 
	mapped_jxn_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if pass_filter[jxn_count] == 0:  #jxn was unmapped in hg19
			jxn_count = jxn_count + 1
			continue
		jxn_count = jxn_count + 1
		hg19_coordinate = hg19_coordinates[mapped_jxn_count]
		mapped_jxn_count = mapped_jxn_count + 1
		cluster_id = data[0].split(':')[-1]
		#if cluster_id not in valid_clusters: #cluster did not have more than 1 jxn
	#		continue
		###Passed all filters
		info = hg19_coordinate.split(':')
		start = int(info[1]) -1
		end = int(info[2]) +1
		new_hg19_coordinate = info[0] + ':' + str(start) + ':' + str(end)
		t.write(new_hg19_coordinate +'\t' + '\t'.join(data[1:]) + '\n')
	t.close()


#Convert data from hg38 to hg19
def liftover_from_hg38_to_hg19(cluster_file_hg38,cluster_file_hg19,liftover_directory,output_dir,tissue_type):
	#First make tempory bed file of our existing hg38 cluster positions. Will be used as input to liftover
	temporary_hg38_bed_file = output_dir + tissue_type + 'temp_hg38.bed'
	make_input_bed_file(cluster_file_hg38,temporary_hg38_bed_file)
	temporary_hg19_bed_file = output_dir + tissue_type + 'temp_hg19.bed' #temporary liftover output file
	temporary_missing_file = output_dir + tissue_type + 'temp_liftover_missing.bed' #temporary liftover missing values file
	run_liftover(temporary_hg38_bed_file,temporary_hg19_bed_file,temporary_missing_file,liftover_directory)
	#Some jxns were not able to be mapped with liftover (lacked confidence). So first extract those unmapped jxns
	unmapped_jxns = get_unmapped_jxns(temporary_missing_file)
	#If raw_leafcutter_file has n+1 lines (therefor n jxns), compute binary vector of length n called pass_filter
	#If pass_filter == 1 then that jxn has a mapping in hg19:
	#Also keep track of which clusters are valid (ie have more than one jxn after unmappability filtering)
	pass_filter,valid_clusters = determine_jxns_that_pass_unmappable_filter(cluster_file_hg38,unmapped_jxns)
	#Extract converted hg19 coordinates
	hg19_coordinates = get_hg19_coordinates(temporary_hg19_bed_file)
	#print these new coordinates while filtering on those unmappable ones/ ones with clusters with only one jxn
	print_helper_liftover(cluster_file_hg38,cluster_file_hg19,pass_filter,valid_clusters,hg19_coordinates)
	##########################
	#REMOVE TEMP FILES 
	##########################
	os.system('rm ' + temporary_hg19_bed_file)
	os.system('rm ' + temporary_missing_file)
	os.system('rm ' + temporary_hg38_bed_file)


#input_file
hg38_jxn_rareness_file = sys.argv[1]
#output_file
hg19_jxn_rareness_file = sys.argv[2]
#directory containing all necessary files to perform liftover
liftover_directory = sys.argv[3]

output_dir = sys.argv[4]
tissue_type = sys.argv[5]
liftover_from_hg38_to_hg19(hg38_jxn_rareness_file,hg19_jxn_rareness_file,liftover_directory,output_dir,tissue_type)