import numpy as np
import os
import sys
import pdb
import gzip


#produce list of autosomal chromosomes
def get_valid_chromosomes():
	valid_chromosomes = {}
	for num in range(1,23):
		stringer = 'chr' + str(num)
		valid_chromosomes[stringer] = 1
	return valid_chromosomes

#If raw_leafcutter_file has n+1 lines (therefor n jxns), compute binary vector of length n called pass_read_filter
#If pass_read_filter == 1 then that jxn has:
#####1. at least 1 individual with greater than or equal to min_reads
#####2. Is on an autosomal
def determine_jxns_that_pass_read_filter(raw_leafcutter_file,min_reads):
	pass_read_filter = []
	f = gzip.open(raw_leafcutter_file)
	#produce list of autosomal chromosomes
	valid_chromosomes = get_valid_chromosomes()
	count = 0 #for header
	total_jxns = 0 #keep track of number of jxns in the file
	for line in f:
		line = line.rstrip()
		data = line.split()
		if count == 0:  #skip header
			count = count + 1
			continue
		read_counts_across_samples = np.asarray(data[1:]).astype(float)
		#determine which samples have more reads mapped than min_reads
		expressing_samples = np.where(read_counts_across_samples >= min_reads)[0]
		jxn_chrom_num = data[0].split(':')[0]
		if len(expressing_samples) > 0 and jxn_chrom_num in valid_chromosomes: #more than zero samples have more than min_reads and on autosomal chromosome
			pass_read_filter.append(1)
		else: #Jxn not on autosomal chromosome or no samples have more than min_reads
			pass_read_filter.append(0)
		total_jxns =total_jxns + 1
	if len(pass_read_filter) != total_jxns: ##Just a double check
		print('Error in number of jxns in this file')
		pdb.set_trace()
	return pass_read_filter

#Now determine all of the clusters that have at least 2 valid jxns (after the above filter).
#We will call these 'valid_clusters'. 
def determine_clusters_with_more_than_one_junction(raw_leafcutter_file,pass_read_filter):
	#cluster_counts will keep track of how many junctions are assigned to each cluster. ie key is the cluster, value is the number of jxns (int)
	cluster_counts = {}
	#A list (dictionary) of all valid_clusters 
	valid_clusters = {}
	f = gzip.open(raw_leafcutter_file)
	count = 0 #for header
	jxn_count = 0 #keep track of number jxn we are on in the file
	for line in f:
		if count == 0:  #skip header
			count = count + 1
			continue
		if pass_read_filter[jxn_count] == 0: #This jxn doesn't pass our filters. Ignore
			jxn_count = jxn_count + 1
			continue
		jxn_count = jxn_count + 1
		line = line.rstrip()
		data = line.split()
		cluster_id = data[0].split(':')[-1] #leafcutter cluster id
		#Keep track of how many jxns each cluter has
		if cluster_id not in cluster_counts: 
			cluster_counts[cluster_id] = 1
		else:
			cluster_counts[cluster_id] = cluster_counts[cluster_id] + 1
	f.close()
	#Extract only clusters with more than one jxn
	for cluster in cluster_counts.keys():
		if cluster_counts[cluster] > 1:
			valid_clusters[cluster] =1
	return valid_clusters

#Convert sample_ids to individual ids
#Put individual ids in alphabetical order
#Then create array of length(ids) where each element corresponds to the sample's/individual's position in the alphabetical array
def get_alphabetical_ordering(sample_ids):
	individual_ids = []
	#loop through sample ids
	for sample_id in sample_ids:
		#extract individual id from sample id (its just the first two strings)
		individual_id = sample_id.split('-')[0] + '-' + sample_id.split('-')[1]
		individual_ids.append(individual_id)
	#alphabetically sort individaul ids
	individual_ids = np.asarray(individual_ids)
	sorted_individual_ids = np.asarray(sorted(individual_ids))
	alphabetical_ordering = []
	#loop through sorted individual ids
	for individual_id in sorted_individual_ids:
		new_position = np.where(individual_ids == individual_id)[0][0]
		alphabetical_ordering.append(new_position)
	if np.array_equal(individual_ids[alphabetical_ordering],sorted_individual_ids) == False:
		print('Error in creating alphabetical ordering') 
	return alphabetical_ordering,sorted_individual_ids

#Subtract 1 from the 5' splice site position in order for us to equal Derek's leafcutter output 
def reformat_jxn_id(old_id):
	info = old_id.split(':')
	if len(info) != 4:
		print('error in reformatting junction id')
	old_5_pos = int(info[1])
	new_5_pos = str(old_5_pos -1)
	new_jxn_id = info[0] + ':' + new_5_pos + ':' + info[2] + ':' + info[3]
	return new_jxn_id

#Now we want to print this filtered list of jxns/ clusters
#We are also (for organization) going to:
####1. Reorder the columns and put them in alphabetical order by individual (not sample) id
####2.Subtract 1 from the 5' splice site position in order for us to equal Derek's leafcutter output 
############(NOTE: this means our coordinates are different from SNAPTRON). If snaptron is coordinate (start,end). We are now (start-1,end+1)
def print_filtered_junctions(raw_leafcutter_file,pass_read_filter,valid_clusters,output_file):
	f = gzip.open(raw_leafcutter_file)
	count = 0 #for header
	jxn_count = 0 #keep track of number jxn we are on in the file
	t = open(output_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split()
		if count == 0:  #header
			count = count + 1
			alphabetical_ordering,sorted_individual_ids = get_alphabetical_ordering(data)
			t.write('Jxn_id\t' + '\t'.join(sorted_individual_ids) + '\n')
			continue
		if pass_read_filter[jxn_count] == 0: #This jxn doesn't pass our filters. Ignore
			jxn_count = jxn_count + 1
			continue
		jxn_count = jxn_count + 1
		cluster_id = data[0].split(':')[-1]
		if cluster_id not in valid_clusters: #This jxn doesn't pass our filters. Ignore
			continue
		#Subtract 1 from the 5' splice site position in order for us to equal Derek's leafcutter output 
		new_jxn_id = reformat_jxn_id(data[0])
		count_array = np.asarray(data[1:])
		#Re-arrange array of counts so it is in the same order as the alphabetically ordered header
		alphabetical_order_count_array = count_array[alphabetical_ordering]
		#print line
		t.write(new_jxn_id + '\t' + '\t'.join(alphabetical_order_count_array) + '\n')
	t.close()

def filter_low_expressing_junctions(tissue_type,raw_leafcutter_file,output_file,min_reads):
	#If raw_leafcutter_file has n+1 lines (therefor n jxns), compute binary vector of length n called pass_read_filter
	#If pass_read_filter == 1 then that jxn has:
	#####1. at least 1 individual with greater than or equal to min_reads
	#####2. Is on an autosomal!!
	pass_read_filter = determine_jxns_that_pass_read_filter(raw_leafcutter_file,min_reads)
	#Now determine all of the clusters that have at least 2 valid jxns (after the above filter).
	#We will call these 'valid_clusters'. 
	valid_clusters = determine_clusters_with_more_than_one_junction(raw_leafcutter_file,pass_read_filter)
	#Now we want to print this filtered list of jxns/ clusters
	#We are also (for organization) going to:
	####1. Reorder the columns and put them in alphabetical order by individual (not sample) id
	####2.Subtract 1 from the 5' splice site position in order for us to equal Derek's leafcutter output 
	############(NOTE: this means our coordinates are different from SNAPTRON). If snaptron is coordinate (start,end). We are now (start-1,end+1)
	print_filtered_junctions(raw_leafcutter_file,pass_read_filter,valid_clusters,output_file)

#First make tempory bed file of our existing hg38 cluster positions. Will be used as input to liftover
def make_input_bed_file(cluster_file_hg38,temporary_hg38_bed_file):
	t = open(temporary_hg38_bed_file,'w')
	f = open(cluster_file_hg38)
	count = 0 #header
	for line in f:
		line = line.rstrip()
		data = line.split()
		if count == 0: #skip header
			count = count + 1
			continue
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
		if count == 0:  #skip header
			count = count + 1
			continue
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
		if count == 0:
			count = count + 1
			t.write(line + '\n')
			continue
		if pass_filter[jxn_count] == 0:  #jxn was unmapped in hg19
			jxn_count = jxn_count + 1
			continue
		jxn_count = jxn_count + 1
		hg19_coordinate = hg19_coordinates[mapped_jxn_count]
		mapped_jxn_count = mapped_jxn_count + 1
		cluster_id = data[0].split(':')[-1]
		if cluster_id not in valid_clusters: #cluster did not have more than 1 jxn
			continue
		###Passed all filters
		t.write(hg19_coordinate + ':' + cluster_id + '\t' + '\t'.join(data[1:]) + '\n')
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





tissue_type = sys.argv[1]
raw_leafcutter_file = sys.argv[2] #unprocessed leafcutter output file (step 2)
clusters_filter_output_dir = sys.argv[3]
#For a jxn to be kept, it must have at least one sample with greater than or equal to min_reads
min_reads = int(sys.argv[4])
#Directory that contains necessary liftover information.
##Specifically, it must contain:
#####1. 'liftOver'   --> the executable
#####2. 'hg19ToHg38.over.chain.gz'   --> for converting from hg19 to hg38
#####2. 'hg38ToHg19.over.chain.gz'   --> for converting from hg38 to hg19
liftover_directory = sys.argv[5]

#######################################
#Filter out jxns that are:
###1. On non-autosomal chromosomes
###2. Have no individuals with greater than min_reads mapped
#We are also (for organization) going to:
####1. Reorder the columns and put them in alphabetical order by individual (not sample) id
####2.Subtract 1 from the 5' splice site position in order for us to equal Derek's leafcutter output 
############(NOTE: this means our coordinates are different from SNAPTRON). If snaptron is coordinate (start,end). We are now (start-1,end+1)
#Always apply filter after to ensure constraint that a cluster has more than one junction. If cluster has only 1 jxn, then remove
output_file_filter_clusters = clusters_filter_output_dir + tissue_type + '_hg38_filtered.txt'
filter_low_expressing_junctions(tissue_type,raw_leafcutter_file,output_file_filter_clusters,min_reads)
#Convert from hg38--> hg19 using liftover.
#Note: some reads are unmappable. This will leave some clusters with only one jxn. Filter again on clusters to make sure they all have more than 1 jxn.
output_file_filter_clusters_lift_over = clusters_filter_output_dir + tissue_type + '_hg19_filtered.txt'
liftover_from_hg38_to_hg19(output_file_filter_clusters,output_file_filter_clusters_lift_over,liftover_directory,clusters_filter_output_dir,tissue_type)




