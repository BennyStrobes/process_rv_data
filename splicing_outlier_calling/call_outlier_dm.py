import sys
import pdb
import os
import numpy as np 
import time
import gzip
import dm_glm
from scipy import stats
import scipy.stats as ss



#simple helper method to add key to dictionary
def add_jxn_to_gene_jxn_data_structure(gene_jxn_data_structure,gene,jxn_read_counts):
    if gene not in gene_jxn_data_structure:
        gene_jxn_data_structure[gene] = []
    gene_jxn_data_structure[gene].append(jxn_read_counts)
    return gene_jxn_data_structure

#Extract mapping from gene to NXK jxn matrix. 
#Also extract ordered array of samples corresponding to max possible samples
def extract_raw_gene_jxn_data_strucutre(tissue_specific_jxn_file):
    f = open(tissue_specific_jxn_file)
    count =0
    gene_jxn_data_structure = {}
    #loop through jxns
    for line in f:
        line = line.rstrip()
        data = line.split()

        #header
        if count == 0:
            count = count + 1
            #unfiltered samples
            temp_samples = data[1:]
            samples = []
            indices = []
            #filter temp_samples based on used_samples
            for i,val in enumerate(temp_samples):
                samples.append(val)
                indices.append(i)
            if np.array_equal(np.asarray(temp_samples)[indices],np.asarray(samples)) == False:
                print('erororor')
            
                pdb.set_trace()
            continue
            


        #################
        #Extract jxn info from current line
        key = data[0]
        info = key.split(':')
        chrom_num = info[0] 
        start = info[1] #5' ss
        end = info[2] #3' ss
        cluster_id = info[3] #name of cluster
        genes = info[4].split(',') #array of genes the jxn is mapped to
        jxn_read_counts = np.asarray(data[1:]).astype(float)[indices]
        #Loop through all genes for this jxn

        cluster = info[3]
        gene_jxn_data_structure = add_jxn_to_gene_jxn_data_structure(gene_jxn_data_structure,cluster,jxn_read_counts)
    f.close()
    #convert list of arrays to matrix.
    for gene in gene_jxn_data_structure.keys():
        jxn_matrix = np.asmatrix(gene_jxn_data_structure[gene])
        gene_jxn_data_structure[gene] = {}
        gene_jxn_data_structure[gene]['samples'] = np.asarray(samples)
        gene_jxn_data_structure[gene]['jxn_matrix'] = np.transpose(jxn_matrix)
    return gene_jxn_data_structure,np.asarray(samples)

#extreme_junctions: Take max_dm_junctions/2 that have the highest read count and take the max_dm_junctions/2 that have the smallest read counts.
def max_number_of_jxns_filter_extreme_junctions(gene_jxn_data_structure,max_dm_junctions):
    new_jxn_structure = {}
    for gene in gene_jxn_data_structure.keys():
        samples = gene_jxn_data_structure[gene]['samples']
        jxn_mat = gene_jxn_data_structure[gene]['jxn_matrix']
        N,K = jxn_mat.shape
        if K > max_dm_junctions: #More jxns than max... Need to filter
            num_top_jxns = int(max_dm_junctions/2)
            num_bottom_jxns = int(max_dm_junctions/2)
            jxn_counts = np.squeeze(np.asarray(np.sum(jxn_mat,axis=0)))
            #sort jxns based on counts across individuals
            #sorted_jxns[0] is the jxn with the fewest number of reads mapped
            #sorted_jxns[-1] is the jxn with the most number of reads mapped
            sorted_jxns = jxn_counts.argsort()
            #new indices
            new_jxn_indices = sorted(np.hstack((sorted_jxns[-num_top_jxns:],sorted_jxns[0:num_bottom_jxns])))
            new_jxn_mat = jxn_mat[:,new_jxn_indices]
            new_jxn_structure[gene] = {}
            new_jxn_structure[gene]['samples'] = samples
            new_jxn_structure[gene]['jxn_matrix'] = new_jxn_mat
        else: #Less or equal to
            new_jxn_structure[gene] = {}
            new_jxn_structure[gene]['samples'] = samples
            new_jxn_structure[gene]['jxn_matrix'] = jxn_mat
    return new_jxn_structure

#ignore_genes: Don't include genes if they contain more jxns than max_dm_junctions.
def max_number_of_jxns_filter_ignore_genes(gene_jxn_data_structure,max_dm_junctions):
    new_jxn_structure = {}
    for gene in gene_jxn_data_structure.keys():
        samples = gene_jxn_data_structure[gene]['samples']
        jxn_mat = gene_jxn_data_structure[gene]['jxn_matrix']
        N,K = jxn_mat.shape
        if K > max_dm_junctions: #More jxns than max... Need to filter
            pass
        else: #Less or equal to
            new_jxn_structure[gene] = {}
            new_jxn_structure[gene]['samples'] = samples
            new_jxn_structure[gene]['jxn_matrix'] = jxn_mat
    return new_jxn_structure

#At the gene level, require an individuals to have at least $min_reads_per_individual.
#If there are less than $min_individuals_per_gene, discard gene
def min_reads_per_individual_filter(gene_jxn_data_structure,min_reads_per_individual,min_individuals_per_gene):
    new_jxn_structure = {}
    for gene in gene_jxn_data_structure.keys():
        samples = gene_jxn_data_structure[gene]['samples']
        jxn_mat = gene_jxn_data_structure[gene]['jxn_matrix']
        N,K = jxn_mat.shape
        new_samples = []
        new_jxn_mat = []
        for n in range(N):
            #vector spanning jxns of number of read counts for this sample
            sample_counts = np.squeeze(np.asarray(jxn_mat[n,:]))
            sample_id = samples[n]
            total_sample_read_counts = np.sum(sample_counts)
            if total_sample_read_counts > min_reads_per_individual:
                new_samples.append(sample_id)
                new_jxn_mat.append(sample_counts)
        #convert from list of arrays to matrix
        new_jxn_mat = np.asmatrix(new_jxn_mat)
        new_samples = np.asarray(new_samples)
        N_new,K_new = new_jxn_mat.shape
        #WIf there are less than $min_individuals_per_gene, don't include gene
        if N_new >= min_individuals_per_gene:
            if N_new != len(new_samples):
                print('EROROROOROROR')
            new_jxn_structure[gene] = {}
            new_jxn_structure[gene]['samples'] = new_samples
            new_jxn_structure[gene]['jxn_matrix'] = new_jxn_mat
    return new_jxn_structure

####1. Limit K to be at max: max_dm_junctions. Implement this filtering by jxn_filter_method
def max_number_of_jxns_filter(gene_jxn_data_structure,max_dm_junctions,jxn_filter_method):
    if jxn_filter_method == 'extreme_junctions':
        gene_jxn_data_structure = max_number_of_jxns_filter_extreme_junctions(gene_jxn_data_structure,max_dm_junctions)
    elif jxn_filter_method == 'ignore_genes':
        gene_jxn_data_structure = max_number_of_jxns_filter_ignore_genes(gene_jxn_data_structure,max_dm_junctions)
    return gene_jxn_data_structure

def order_data(sample_to_pos,data_vec,arr):
    new_arr = []
    for data in arr:
        pos = sample_to_pos[data[0]]
        new_info = data_vec[pos]
        data = data + (new_info,)
        new_arr.append(data)
    return new_arr
def create_ordered_sample_arr(ordered_samples_file):
    f = open(ordered_samples_file)
    arr = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        arr.append((data[0],data[1]))
    return arr

def order_data_update(arr,flagged):
    new_arr = []
    for data in arr:
        binary_flagged = flagged[data[0]]
        new_info = str(binary_flagged)
        data = data + (new_info,)
        new_arr.append(data)
    return new_arr



#get samples pre-filtering
def extract_original_samples(jxn_file):
    f = open(jxn_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        return np.asarray(data[1:])

#Remove samples flagged in the v6 sample attribute file
#Return dictionary of individual ids that is subset of original samples
def filter_samples_on_v6_attribute_file(original_samples,v6_sample_attribute_file,sample_attribute_file_tissue_type,outlier_calling_samples_file,node_number):
    #convert original_samples array into dictionary
    dicti = {}
    for ele in original_samples:
        dicti[ele] = 0
    #Loop through v6 sample atttribute file
    f = open(v6_sample_attribute_file)
    flagged = {}
    count = 0
    for line in f:
        line= line.rstrip()
        data = line.split('\t')
        if count == 0: #skip header
            count = count + 1
            continue
        #convert v6_sample id to indiviudal id
        indi_id = data[0].split('-')[0] + '-' + data[0].split('-')[1]
        #if individual id in original sample and sample id in specified tissue type
        if indi_id in dicti and data[13] == sample_attribute_file_tissue_type and len(data) == 75:
            #not flagged
            if data[28] == 'USE ME':
                flagged[indi_id] =0 
            #glagged
            else:
                flagged[indi_id] = 1
    #simple error checking
    if len(flagged) != len(dicti):
        print('assumption error')
        pdb.set_trace()
    #create dictionary that only contains unflagged samples
    samples_to_use = {}
    for ele in original_samples:
        if flagged[ele] == 0:
            samples_to_use[ele] =1
    return samples_to_use


#Remove samples with greater than 1 STDEV from population mean in terms of number of unique junctions and number of non-zero junctions
def filter_global_outlier_samples(samples_to_use,jxn_file):
    f = open(jxn_file)
    count = 0
    #loop through junction file
    for line in f:
        line = line.rstrip()
        data =line.split()
        #if header == true
        if count == 0:
            count = count + 1
            #extract total individuals
            indiz = np.asarray(data[1:])
            indi_to_ele = {}
            #create dictionary from individaul id to array index (ie unique_counts,nonzero_counts)
            for i,indi in enumerate(indiz):
                indi_to_ele[indi] = i
            #initialize count arrays
            #keep track of number of unique junctions (no other individuals are expressed at this junctions) this individual has
            unique_counts = np.zeros(len(indiz))
            #Keep track of number of nonzero junctions a particular individual has
            nonzero_counts = np.zeros(len(indiz))
            #keep track of total reads mapped to this individual (junction library size)
            total_counts = np.zeros(len(indiz))
            continue
        #Junction identifier
        key = data[0]
        info = key.split(':')
        chrom_num = info[0] 
        start = info[1] #5' ss
        end = info[2] #3' ss
        cluster_id = info[3] #name of cluster
        genes = info[4].split(',') #array of genes the jxn is mapped to
        #The raw read counts for this junction
        jxn_read_counts = np.asarray(data[1:]).astype(float) 
        nonzero_positions = np.where(jxn_read_counts > 0)[0]
        #keep track of how many nonzero junctions each individual has
        for nonzero_position in nonzero_positions:
            nonzero_counts[nonzero_position] = nonzero_counts[nonzero_position] + 1
        #A unique junction has occured
        if len(nonzero_positions) == 1:
            unique_counts[nonzero_positions[0]] = 1 + unique_counts[nonzero_positions[0]]
    f.close()
    #Reduce array dimension from len(individuals) to len(filtered individuals) based on samples_to_use dictionary for purpose of computing accurate mean and standard deviation
    unique_counts2 = []
    nonzero_counts2 = []
    new_samples_to_use = {}
    for indi in indiz:
        if indi not in samples_to_use:
            continue
        unique_counts2.append(unique_counts[indi_to_ele[indi]])
        nonzero_counts2.append(nonzero_counts[indi_to_ele[indi]])
    #We are going to filter using individuals greater than 1 stdevation away from the mean of both unique counts and nonzero counts
    mean_u = np.mean(unique_counts2) #mean of unique counts
    stdev_u = np.std(unique_counts2) #1 stdev of unique counts
    mean_nonz = np.mean(nonzero_counts2) #mean of nonzero counts
    stdev_nonz = np.std(nonzero_counts2) #1 stdev of nonzero counts
    for indi in indiz:
        #only consider individuals in samples_to_use
        if indi not in samples_to_use:
            continue
        #filter
        unique = unique_counts[indi_to_ele[indi]]
        nonzero = nonzero_counts[indi_to_ele[indi]]
        if unique < mean_u + 1.96*stdev_u and unique > mean_u - 1.96*stdev_u and nonzero < mean_nonz + 1.96*stdev_nonz and nonzero > mean_nonz -1.96*stdev_nonz:
            new_samples_to_use[indi] = 1
    return new_samples_to_use


#Remove any samples not in tissue specific v6p covariate file
def filter_samples_on_v6p_covariate_file(samples_to_use,v6p_covariate_file):
    f = open(v6p_covariate_file)
    count = 0
    new_samples_to_use = {}
    for line in f:
        if count == 0: #header
            count = count + 1
            line = line.rstrip()
            data = line.split()
            #extract list of individuals in the v6p covariate file
            indi_cov = data[1:]
            #intersect indi_cov with our previous samples_to_use
            for indi in indi_cov:
                if indi in samples_to_use:
                    new_samples_to_use[indi] = 1
    return new_samples_to_use




#We wish to remove samples that are "global outliers" (ie outliers in many different genes)
#The way we filter these samples (we are actually using individual ids) is determined by the parameter $filter_global_outlier_method
def filter_samples(tissue_specific_jxn_file,filter_global_outlier_method,v6_sample_attribute_file,outlier_calling_samples_file,sample_attribute_file_tissue_type,v6p_covariate_file,node_number):
    #Extract samples pre-filter
    original_samples = extract_original_samples(tissue_specific_jxn_file)


    #Which $filter_global_outlier_method are we to use:...
    #'only_flagged' removes samples flagged by the v6_samples_attribute_file
    if filter_global_outlier_method == 'only_flagged':
        samples_to_use = filter_samples_on_v6_attribute_file(original_samples,v6_sample_attribute_file,sample_attribute_file_tissue_type,outlier_calling_samples_file,node_number)
    #'filter1' removes samples flagged by the v6_samples_attribute_file, but also removes samples with greater than 1 STDEV from population mean in terms of number of unique junctions and number of non-zero junctions
    elif filter_global_outlier_method == 'filter1':
        #remove samples flagged by the v6_samples_attribute_file
        samples_to_use = filter_samples_on_v6_attribute_file(original_samples,v6_sample_attribute_file,sample_attribute_file_tissue_type,outlier_calling_samples_file,node_number)
        #Remove samples with greater than 1 STDEV from population mean in terms of number of unique junctions and number of non-zero junctions
        samples_to_use = filter_global_outlier_samples(samples_to_use,tissue_specific_jxn_file)
    elif filter_global_outlier_method == 'filter2':
        #remove samples flagged by the v6_samples_attribute_file
        samples_to_use = filter_samples_on_v6_attribute_file(original_samples,v6_sample_attribute_file,sample_attribute_file_tissue_type,outlier_calling_samples_file,node_number)
        #Remove samples not in v6p tissue specific covariate file
        samples_to_use = filter_samples_on_v6p_covariate_file(samples_to_use,v6p_covariate_file)
        #Remove samples with greater than 1 STDEV from population mean in terms of number of unique junctions and number of non-zero junctions
        samples_to_use = filter_global_outlier_samples(samples_to_use,tissue_specific_jxn_file)
    #we now wish to print the list of individuals we are using for this analysis
    if node_number == 0: #only print on one of the parrallelized nodes per ex
        t = open(outlier_calling_samples_file,'w')
        for key in sorted(samples_to_use.keys()):
            t.write(key + '\n')
        t.close()
    return samples_to_use

#Extract the covariate matrix from the tissue specific covariate file. Extract for all samples we are currently looking at in this tissue
def extract_covariates_from_v6p_covariate_file(samples,v6p_covariate_file,num_pc):
    f = open(v6p_covariate_file)
    count = 0
    covs = [] #output data structure
    #loop through covariate file
    for line in f:
        line = line.rstrip()
        data = line.split()
        #If header
        if count == 0:
            count = count + 1
            #create mapping (indices) from v6p samples array (indiz) to full samples array (samples)
            indices = []
            indiz = np.asarray(data[1:])
            for sample in samples:
                for i,indi in enumerate(indiz):
                    if indi == sample:
                        indices.append(i)
            #check to make sure that mapping is correct
            if np.array_equal(indiz[indices],samples) == False:
                print('erroorororrorr')
                pdb.set_trace()
            continue
        #Only take top $num_pc principle components
        if data[0].startswith('InferredCov'):
            pc_number = int(data[0].split('dCov')[1])
            if pc_number > num_pc:
                continue
        covariate = np.asarray(data[1:]).astype(float)
        covs.append(covariate)
    covs = np.asmatrix(covs)
    return np.transpose(covs)

#mean center and center columns of the covariate matrix (ie the covariates)
def standardize_matrix(matrix):
    rows,cols = matrix.shape
    new_matrix = np.zeros((rows,cols))
    for i in range(rows):
        for j in range(cols):
            new_matrix[i,j] = (matrix[i,j] - np.mean(matrix[:,j]))/np.std(matrix[:,j])
    return new_matrix

#Add covariate matrix to gene_jxn_data_structure when the covariate regression method == 'v6p_eqtl_covariates'
def add_covariates_to_gene_junction_data_structure_v6p_eqtl_covariates_version(gene_jxn_data_structure,samples,v6p_covariate_file,num_pc):
    #Extract the covariate matrix from the tissue specific covariate file. Extract for all samples we are currently looking at in this tissue
    covariate_matrix = extract_covariates_from_v6p_covariate_file(samples,v6p_covariate_file,num_pc)
    #standardize matrix
    covariate_matrix = standardize_matrix(covariate_matrix)
    #Loop through relevent genes
    for gene in gene_jxn_data_structure.keys():
        #Extract array of samples defined for this gene
        gene_samples = gene_jxn_data_structure[gene]['samples']
        #create mapping (indices) from full samples array (samples) to gene specific samples array (gene_samples)
        indices = []
        for i,sample in enumerate(samples):
            if sample in gene_samples:
                indices.append(i)
        #check to make sure I performed mapping correctly
        if np.array_equal(samples[indices],gene_samples) == False:
            print('eororoorororrororo')
            pdb.set_trace()
        #Add to gene_jxn_data_structure
        gene_jxn_data_structure[gene]['covariate_matrix'] = covariate_matrix[indices,:]
    return gene_jxn_data_structure

def add_covariates_to_gene_jxn_data_structure_general(gene_jxn_data_structure, samples, covariate_file):
    covariate_matrix = np.transpose(np.loadtxt(covariate_file, dtype=str)[1:,1:].astype(float))
    #standardize matrix
    covariate_matrix = standardize_matrix(covariate_matrix)
    #Loop through relevent genes
    for gene in gene_jxn_data_structure.keys():
        #Extract array of samples defined for this gene
        gene_samples = gene_jxn_data_structure[gene]['samples']
        #create mapping (indices) from full samples array (samples) to gene specific samples array (gene_samples)
        indices = []
        for i,sample in enumerate(samples):
            if sample in gene_samples:
                indices.append(i)
        #check to make sure I performed mapping correctly
        if np.array_equal(samples[indices],gene_samples) == False:
            print('eororoorororrororo')
            pdb.set_trace()
        #Add to gene_jxn_data_structure
        gene_jxn_data_structure[gene]['covariate_matrix'] = covariate_matrix[indices,:]
    return gene_jxn_data_structure

#Add covariate matrix to gene_jxn_data_structure. Key is called 'covariate_matrix'
def add_covariates_to_gene_junction_data_structure(gene_jxn_data_structure,samples,covariate_regression_method,covariate_file):
    #Select which version of covariates we wish to use
    if covariate_regression_method != 'none':
        #v6p_eqtl_covariates uses all covariates (except only first num_pc) of the gtex v6p eqtl covariates
        #gene_jxn_data_structure = add_covariates_to_gene_junction_data_structure_v6p_eqtl_covariates_version(gene_jxn_data_structure,samples,v6p_covariate_file,num_pc)
        gene_jxn_data_structure = add_covariates_to_gene_jxn_data_structure_general(gene_jxn_data_structure, samples, covariate_file)
    return gene_jxn_data_structure

#Convert jxn file into a gene based data structure.
#Object is a dictionary with keys that are genes. Values are a dictionary with keys 'samples' and 'jxn_matrix'. 
#Values of jxn_matrix are NxK matrices where N is the number of samples and K is the number of junctions
##Filtering:
####1. Limit K to be at max: max_dm_junctions. Implement this filtering by jxn_filter_method
####2. At the gene level, require an individuals to have at least $min_reads_per_individual.
########Following up on filter 2, if there are less than $min_individuals_per_gene, discard gene

#@Return: gene_jxn_data_structure
#@RETURN: samples --> ordered array of individual ids corresponding to all samples that passed sample specific filters. 
def create_gene_based_data_structure(tissue_specific_jxn_file, max_dm_junctions, jxn_filter_method, min_reads_per_individual, min_individuals_per_gene, outlier_calling_samples_file, covariate_file, covariate_regression_method, node_number):   

    #Get raw data structure
    #Also get samples, this is the maximum possible samples after filtering
    gene_jxn_data_structure, samples = extract_raw_gene_jxn_data_strucutre(tissue_specific_jxn_file)

    #apply max_junction_filter
    gene_jxn_data_structure = max_number_of_jxns_filter(gene_jxn_data_structure, max_dm_junctions, jxn_filter_method)
    #Apply min_reads_per_individual filter
    gene_jxn_data_structure = min_reads_per_individual_filter(gene_jxn_data_structure, min_reads_per_individual, min_individuals_per_gene)

    #Add covariate matrix to gene_jxn_data_structure. Key is called 'covariate_matrix'
    gene_jxn_data_structure = add_covariates_to_gene_junction_data_structure(gene_jxn_data_structure, samples, covariate_regression_method, covariate_file)


    return gene_jxn_data_structure,samples


def create_expression_matrix(data_struct,samples_file,tpm_file,all_samples):
    indi= all_samples
    genes = {}
    for gene in data_struct.keys():
        genes[gene] = 1
    samples = {}
    #create dictionary of len(all_samples) that contains sample ids (as opposed to individual ids present in all_samples)
    f = open(samples_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        indi_id = line.split('-')[0] + '-' + line.split('-')[1]
        if indi_id in all_samples:
            samples[line] = 1
    #error checking
    if len(samples) != len(indi):
        print('ERROR')
    f = gzip.open(tpm_file)
    count = 0
    tpm_counts = {}
    #loop through tpm file
    for line in f:
        line = line.rstrip()
        data = line.split()
        #header
        if count == 0:
            count = count + 1
            indices = []
            for ind in indi:
                #indices converts from tpm order to order of all_samples (indi)
                for i,ele in enumerate(data):
                    if ele in samples and ele.split('-')[0] + '-' + ele.split('-')[1] == ind:
                        indices.append(i)
            if len(indices) != len(indi):
                print('erroororo')
                pdb.set_trace()
            continue
        #gene_id = data[1]
        if data[1] in genes:
            data = np.asarray(data)
            if data[1] not in tpm_counts:
                tpm_counts[data[1]] = data[indices].astype(float)
            else:
                tpm_counts[data[1]] = tpm_counts[data[1]] + data[indices].astype(float)
    tpm_mat = np.zeros((len(tpm_counts),len(indices)))
    for i,gene in enumerate(tpm_counts.keys()):
        tpm_mat[i,:] = tpm_counts[gene]
    return np.transpose(tpm_mat)

def quantile_normalize(count_matrix):
    count_matrix = np.transpose(count_matrix)
    K,N = count_matrix.shape
    new_count_mat = np.zeros((K,N))
    for k in range(K):
        expr_vec = np.squeeze(np.asarray(count_matrix[k,:]))
        ranks = ss.rankdata(expr_vec)/(len(expr_vec) + 1)
        new_exprr = stats.norm.ppf(ranks)
        new_count_mat[k,:] = new_exprr
    return np.transpose(new_count_mat)



#Extract mapping from gene to NXK jxn matrix. No filters here
def extract_raw_counts(tissue_specific_jxn_file):
    f = open(tissue_specific_jxn_file)
    count =0
    jxn_mat = []
    #loop through jxns
    for line in f:
        line = line.rstrip()
        data = line.split()
        if count == 0:
            count = count + 1
            samples = data[1:]
            continue
        #################
        #Extract jxn info from current line
        key = data[0]
        info = key.split(':')
        chrom_num = info[0] 
        start = info[1] #5' ss
        end = info[2] #3' ss
        cluster_id = info[3] #name of cluster
        genes = info[4].split(',') #array of genes the jxn is mapped to
        jxn_read_counts = np.asarray(data[1:]).astype(float) 
        jxn_mat.append(jxn_read_counts)
    jxn_mat = np.asmatrix(jxn_mat)
        #Loop through all genes for this jxn
    return jxn_mat,np.asarray(samples)
#Compute rank of outlier sample. Rank refers to the rank of total read counts of outlier sample divided by length of samples
def get_max_sample_ranking(read_counts,read_count_max):
    ranking = len(np.where(read_counts < read_count_max)[0])
    total = float(len(read_counts))
    return ranking/total


def compute_dm_covariance_matrix(n,alpha):
    alpha_0 = np.sum(alpha)
    p = alpha/alpha_0
    cov = n*((n+alpha_0)/(1+alpha_0))*(np.diag(p) - np.dot(np.transpose(np.asmatrix(p)),np.asmatrix(p)))
    return cov


#Compute mahalanobis distance for one sample
def mahalanobis_distance(x,alpha):
    n = np.sum(x)
    alpha_0 = np.sum(alpha)
    cov = compute_dm_covariance_matrix(n,alpha)
    mu = n*alpha/alpha_0
    diff_mat = np.asmatrix(x-mu)
    distance = np.dot(np.dot(diff_mat,np.linalg.pinv(cov)),np.transpose(diff_mat))[0,0]
    if distance < 0:
        print('ERROR: Mahalanobis distance is less than zero')
        pdb.set_trace()
    return np.sqrt(distance)

#Compute mahalanobis distance for one sample
def mahalanobis_distance_computer(x,alpha):
    n = np.sum(x)
    alpha_0 = np.sum(alpha)
    cov = compute_dm_covariance_matrix(n,alpha)
    K = cov.shape[0]
    cov = cov + np.eye(K)*.1 # padd the covariance matrix
    mu = n*alpha/alpha_0
    diff_mat = np.asmatrix(x-mu)
    p_inv_cov = np.linalg.pinv(cov)
    distance = np.dot(np.dot(diff_mat,p_inv_cov),np.transpose(diff_mat))[0,0]
    error_term = False
    if distance < 0:
        print('ERROR: Mahalanobis distance is less than zero')
        print('cov:')
        print(cov)
        print('x:')
        print(x)
        print('alpha:')
        print(alpha)
        error_term = True
    return np.sqrt(distance), error_term

#Compute mahalanobis distance for each sample. Return array of len == num_samples
def get_mahala_disty_shell(alpha,X):
    N,k = X.shape
    distyz = []
    for n in range(N):
        distyz.append(mahalanobis_distance(X[n,:],alpha))
    return np.asarray(distyz)

#Take num_draws for the fitted DM. For each draw, compute mahalanobis distance
def sample(alpha,num_samples,N_count):
    x = []
    sample_alphas = np.random.dirichlet(alpha,size=num_samples)
    for i in range(num_samples):
        x.append(np.random.multinomial(N_count,sample_alphas[i,:]))
    x = np.asmatrix(x)
    distances = get_mahala_disty_shell(alpha,x)
    return np.squeeze(np.asarray(distances))

def adaptive_sample(alpha,N_count,max_observed_distance):
    num_samples = 1000
    sample_distances = sample(alpha, num_samples, N_count)
    if max(sample_distances) < max_observed_distance:
        num_samples = 10000
        sample_distances = sample(alpha,num_samples,N_count)
        if max(sample_distances) < max_observed_distance:
            num_samples = 100000
            sample_distances = sample(alpha,num_samples,N_count)
        #if max(sample_distances) < max_observed_distance:
            #num_samples = 1000000
            #sample_distances = sample(alpha,num_samples,N_count)
    return sample_distances

def sample_glm(alpha, num_samples, N_count, mu, inv_cov):
    distyz = []
    sample_alphas = np.random.dirichlet(alpha,size=num_samples)

    x = []
    for i in range(num_samples):
        temp_sample_x = np.random.multinomial(N_count,sample_alphas[i,:])
        x.append(temp_sample_x)
    x = np.asmatrix(x)
    if num_samples <= 1000:  # If number of samples is 1000 or less, do not split up the matrix
        diff_mat = x-mu
        distances = np.diag(np.dot(np.dot(diff_mat, inv_cov),np.transpose(diff_mat)))
    else:  # if number of samples is > 1000, split up the matrix into 10 splits (for memory reasons)
        x_split = np.split(x,10)  # Create list of split matrices
        distances = []
        for split_iter in range(len(x_split)):  # For each split matrix, compute the mahalanobis distances
            x_iter = x_split[split_iter]  # The split matrix
            diff_mat = x_iter - mu
            distances.append(np.diag(np.dot(np.dot(diff_mat, inv_cov),np.transpose(diff_mat))))
        distances = np.hstack(distances)
    return np.sqrt(distances)


def adaptive_sample_glm(alpha, md, N_count):
    # N_count = 20000
    alpha_0 = np.sum(alpha)
    cov = compute_dm_covariance_matrix(N_count, alpha)
    K = cov.shape[0]
    cov = cov + np.eye(K)*.1 # padd the covariance matrix
    mu = N_count*alpha/alpha_0
    inv_cov = np.linalg.pinv(cov)

    num_samples = 1000
    sample_distances = sample_glm(alpha, num_samples, N_count, mu, inv_cov)
    if max(sample_distances) < md:
        num_samples = 10000
        sample_distances = sample_glm(alpha, num_samples, N_count, mu, inv_cov)
        #if max(sample_distances) < md:
        #    num_samples = 100000
        #    sample_distances = sample_glm(alpha, num_samples, N_count, mu, inv_cov)
    pval = len(np.where(md <= sample_distances)[0])/float(len(sample_distances))
    return pval

#Take num_draws for the fitted DM. For each draw, compute mahalanobis distance
def sample_rand(alpha,num_samples,counts):
    x = []
    rand_index = np.random.randint(0,len(counts),num_samples)
    sample_alphas = np.random.dirichlet(alpha,size=num_samples)
    for i in range(num_samples):
        N_count = counts[rand_index[i]]
        x.append(np.random.multinomial(N_count,sample_alphas[i,:]))
    x = np.asmatrix(x)
    distances = get_mahala_disty_shell(alpha,x)
    return np.squeeze(np.asarray(distances))


#Run outlier analysis using Dirichlet multinomial for 1 gene
def outlier_analysis_dm(X,samples):
    #Fit dirichlet multinomial to X (learning parameter alpha)
    alpha = dm_glm.dirichlet_multinomial_fit(X)
    N,K = X.shape
    if K != len(alpha):
        print('erorororo')
    #Compute mahalanobis distance for each sample. Return array of len == num_samples
    distances = get_mahala_disty_shell(alpha,X)
    #Compute total number of read counts (summed across jxns) for each sample
    read_counts = np.squeeze(np.asarray(np.sum(X,axis=1)))
    #Most extreme point
    max_observed_distance = max(distances)
    outlier_sample_num = np.where(distances==max_observed_distance)[0][0]
    #Compute rank of outlier sample. Rank refers to the rank of total read counts of outlier sample divided by length of samples
    #outlier_sample_ranking = get_max_sample_ranking(read_counts,read_counts[outlier_sample_num])
    #outlier_sample_id = samples[outlier_sample_num]
    #########################################################
    #estimate emperical distribution
    #########################################################
    #Take num_draws for the fitted DM. For each draw, compute mahalanobis distance
    
    #num_draws=10000
    #sample_distances = sample(alpha,num_draws,read_counts[outlier_sample_num])
    #sample_distances = sample_rand(alpha,num_draws,read_counts)'
    sample_distances = adaptive_sample(alpha,20000,max_observed_distance)

    #Compute nominal pvalues for each sample
    pvalz = [] #array of length number of samples
    for distance in distances:
        pval = len(np.where(distance <= sample_distances)[0])/float(len(sample_distances))
        if pval > 1:
            print('pvalue eroroororor')
            pdb.set_trace()
        pvalz.append(pval)






    return distances,pvalz,alpha

def outlier_analysis_dm_glm(X, alpha_matrix):
    N1,K1 = X.shape
    N2,K2 = alpha_matrix.shape
    N = N1
    if N1 != N2 or K1 != K2:  # make sure shape of X and alpha are the same
        print('Shape mismatch erroror!!')
        pdb.set_trace()
    read_counts = np.squeeze(np.asarray(np.sum(X,axis=1)))  # Compute n (sum(x_i)) for each sample
    distances = []  # Keep track of raw mahalanobis distance for each sample
    pvalues = []  # Keep track of pvalues for each sample
    error = False
    for n in range(N):  # Loop through the samples
        alpha_n = np.squeeze(np.asarray(alpha_matrix[n,:]))  # alpha for this sample
        num_reads = int(np.sum(X[n,:]*500))
        md_n, error_term = mahalanobis_distance_computer(X[n,:]*500, alpha_n)

        if error_term == True:
            print('mahab error')
            num = 5.0/0
        distances.append(md_n)
        pvalue_n = adaptive_sample_glm(alpha_n, md_n, num_reads)
        pvalues.append(pvalue_n)
    return distances, pvalues 



###################################################################################################################################################
#Regress out effects of covariates using Direchlet multinomial GLM 
#Essentially fit GLM. Predict mean of each sample. Take residuals of samples
#Code curtosy of leafcutter (https://github.com/davidaknowles/leafcutter)
def pystan_covariate_regression(X,samples,sample_attribute_file,tissue_type):
    #create matrix of covariates. This fxn will need to be changed
    cov = extract_covariates(samples,sample_attribute_file,tissue_type)
    ##error check
    N1,K1 = X.shape
    N2,K2 = cov.shape
    if N1 != N2:
        print('eororororoor')
    return dm_glm.regress_covariates_with_dm_glm(X,np.hstack((np.ones((N1,1)),cov)))

    




def filter_mat(mat):
    N,K = mat.shape
    good_colz = []
    for k in range(K):
        if np.sum(np.isnan(mat[:,k])) == 0:
            good_colz.append(k)
    mat = mat[:,good_colz]
    stds = np.squeeze(np.asarray(np.std(mat,axis=0)))
    good_colz = []
    for i,std in enumerate(stds):
        if std != 0.0:
            good_colz.append(i)
    mat = mat[:,good_colz]
    N,K = mat.shape
    good_vecs = []
    good_columns = []
    for k in range(K):
        independent = 1
        current_vec = np.squeeze(np.asarray(mat[:,k]))
        for vec in good_vecs:
            corr = abs(np.corrcoef(vec,current_vec)[0,1])
            if corr > .7:
                independent = 0
        if independent == 1:
            good_columns.append(k)
            good_vecs.append(current_vec)
    mat = mat[:,good_columns]
    mat = standardize_mat(mat)
    return mat
def standardize_mat(mat):
    N,K = mat.shape
    for k in range(K):
        vec = mat[:,k]
        mean = np.mean(vec)
        std = np.std(vec)
        mat[:,k] = (vec-mean)/std
    return mat

def extract_covariates(samples,sample_attribute_file,tissue_type):
    sample_dict = {}
    for sample in samples:
        sample_dict[sample] = 1
    count = 0
    f = open(sample_attribute_file)
    for line in f:
        line =line.rstrip()
        data = line.split('\t')
        if count == 0:
            count = count + 1
            continue
        indi_id = data[0].split('-')[0] + '-' + data[0].split('-')[1]
        line_tissue = data[13]
        if line_tissue == tissue_type and indi_id in sample_dict and data[27] == 'RNASEQ':
            data = np.asarray(data)
            if len(data) != 74:
                print('erororor')
            feature_vec = np.hstack((data[[8,16,17]],data[29:74]))
            for i,ele in enumerate(feature_vec):
                if ele == '':
                    feature_vec[i] = 'NaN'
            sample_dict[indi_id] = np.squeeze(np.asarray(feature_vec.astype(float)))
    arr = []
    for indi in samples:
        arr.append(sample_dict[indi])
    mat = np.asmatrix(arr)
    mat = filter_mat(mat)
    return mat








#######################################################################################################################################################


#Print outlier calling dm results to output file
def outlier_calling_print_helper(arr,samples,all_samples,t,gene):
    #simple error checking
    if len(arr) != len(samples):
        print('Print helper error')
        pdb.set_trace()
    counter = 0
    #print row id
    t.write(gene)
    #loop through all samples
    for sample in all_samples:
        #If sample in gene specific samples
        if sample in samples:
            #print value corresponding to that sample
            ele = str(arr[counter])
            counter = counter + 1
        #If sample was filter out (not in gene specific sample). Print NaN
        else:
            ele = 'NaN'
        t.write('\t' + ele)
    t.write('\n')
    if counter != len(samples):
        print('Print helper error 2')
        pdb.set_trace()
    t.flush()
    return t

#call outliers for each gene using a fitted dirichlet multinomial
#Also write to output
def call_outliers_with_dirichlet_multinomial(tissue_specific_outlier_file_root, gene_jxn_data_structure, all_samples, start_number, end_number, covariate_regression_method):

    #Initialize output files
    t_MD = open(tissue_specific_outlier_file_root + '_md.txt','w')#file handle for matrix of mahalnobis distances
    t_pvalue = open(tissue_specific_outlier_file_root + '_emperical_pvalue.txt','w') #file handle for matrix of pvalues 
    #Write headers for output files
    t_MD.write('CLUSTER_ID\t' + '\t'.join(all_samples) + '\n')
    t_pvalue.write('CLUSTER_ID\t' + '\t'.join(all_samples) + '\n')



    start_time = time.time()

    for counter,gene in enumerate(sorted(gene_jxn_data_structure.keys())):
        #Parrallelization Stuff
        if counter < start_number or counter > end_number:
            continue
        print('$$$$$$$$$$$$$$$$$$$$$$$$$$')
        print(counter)
        print('$$$$$$$$$$$$$$$$$$$$$$$$$$')
        ####################################################################
        #time related stuff
        ####################################################################
        curr_time = time.time()
        diff = curr_time - start_time
        start_time = curr_time
        time_min = diff/60.0
        print(time_min)
        print('##################')
        print('START ' + gene)
        ####################################################################
        #Actual analysis
        ####################################################################
        #Extract jxn matrix for gene
        X = gene_jxn_data_structure[gene]['jxn_matrix']
        #Extract genes samples
        samples = gene_jxn_data_structure[gene]['samples']



        working = True
        iteration = 1
        while working:
            try:

                if covariate_regression_method == 'none':
                    #Run outlier analysis
                    N1, J = X.shape
                    # Covariance matrix is just intercept
                    cov_mat = np.ones((N1, 1))

                elif covariate_regression_method == 'junction_pc_no_regress_reg':
                    cov_mat_no_intercept = gene_jxn_data_structure[gene]['covariate_matrix']
                    N1, K = cov_mat_no_intercept.shape
                    N2, J = X.shape
                    if N1 != N2:
                        print('ASSUMPTION ERROROR')
                    cov_mat = np.hstack((np.ones((N1, 1)), cov_mat_no_intercept))
                elif covariate_regression_method == 'junction_pc_only_no_regress_reg':
                    cov_mat_no_intercept = gene_jxn_data_structure[gene]['covariate_matrix']
                    N1, K = cov_mat_no_intercept.shape
                    N2, J = X.shape
                    if N1 != N2:
                        print('ASSUMPTION ERROROR')
                    cov_mat_no_intercept = cov_mat_no_intercept[:,(-cov_mat_no_intercept.shape[1]):-4]
                    cov_mat = np.hstack((np.ones((N1, 1)), cov_mat_no_intercept))

                alpha_matrix, beta_mat, conc_vec = dm_glm.dirichlet_multinomial_glm_fit(X, cov_mat)
                distances, pvalz = outlier_analysis_dm_glm(X, alpha_matrix)
                working = False
                converged = True
            except:
                print('iteration ' + str(iteration))
                iteration = iteration + 1
                if iteration > 20:
                    working = False
                    converged = False
                    print('convergence error for ' + str(counter))

        if converged == True:
            #print Mahalanobis distance results to output file
            t_MD = outlier_calling_print_helper(distances, samples, all_samples, t_MD, gene)
            #print emperical pvalu results to output file
            t_pvalue = outlier_calling_print_helper(pvalz, samples, all_samples, t_pvalue, gene)




def parallelization_start_and_end(total,node_number,total_nodes):
    genes_per_node = (total/total_nodes) +1
    gene_start = node_number*genes_per_node
    gene_end = (node_number + 1)*genes_per_node -1
    return gene_start,gene_end


def call_outliers(tissue_type, tissue_specific_jxn_file, tissue_specific_outlier_file, outlier_calling_dm_output_dir, max_dm_junctions, jxn_filter_method, node_number, total_nodes, covariate_regression_method, covariate_file, rna_seq_samples_file, num_pc, outlier_calling_samples_file):
    min_reads_per_individual = 5
    min_individuals_per_gene = 50
    #Convert jxn file into a gene based data structure.
    #Object is a dictionary with keys that are genes. Values are NxK matrices where N is the number of samples and K is the number of junctions
    ##Filtering:
    ####1. Limit K to be at max: max_dm_junctions. Implement this filtering by jxn_filter_method
    ####2. At the gene level, require an individuals to have at least $min_reads_per_individual.
    ########Following up on filter 2, if there are less than $min_individuals_per_gene, discard gene
    gene_jxn_data_structure, all_samples = create_gene_based_data_structure(tissue_specific_jxn_file, max_dm_junctions, jxn_filter_method, min_reads_per_individual, min_individuals_per_gene, outlier_calling_samples_file, covariate_file, covariate_regression_method, node_number)

    #For parallelization purposes 
    start_number,end_number = parallelization_start_and_end(len(gene_jxn_data_structure),node_number,total_nodes)
    #call outliers for each gene using a fitted dirichlet multinomial
    #Also write to output
    call_outliers_with_dirichlet_multinomial(tissue_specific_outlier_file,gene_jxn_data_structure,all_samples,start_number,end_number, covariate_regression_method)



tissue_type = sys.argv[1]
tissue_specific_jxn_file = sys.argv[2] #input file. Made in generate_junctions.sh
tissue_specific_outlier_file = sys.argv[3] #main output file
outlier_calling_dm_output_dir = sys.argv[4] #output_directory
max_dm_junctions = float(sys.argv[5]) #Maximum number of jxns we will allow a gene to have
#Method to select max_dm_junctions.Current methods include:
###1. extreme_junctions: Take max_dm_junctions/2 that have the highest read count and take the max_dm_junctions/2 that have the smallest read counts.
###2. ignore_genes: Disregard genes that have more than $max_dm_junctions
jxn_filter_method = sys.argv[6]
##############################
#For parallelization purposes
node_number = int(sys.argv[7])
total_nodes = int(sys.argv[8])
##############################
#How to deal with covariates
covariate_regression_method= sys.argv[9]


# covariate file
covariate_file = sys.argv[10]

#all rna-seq sample ids used in this tissue
rna_seq_samples_file = sys.argv[11]
#Number of principle components to use. Only applicable if covariate_regression_method includes pcs
num_pc = int(sys.argv[12])


#output file to write actual samples used in this analysis
outlier_calling_samples_file = sys.argv[13]



call_outliers(tissue_type, tissue_specific_jxn_file, tissue_specific_outlier_file, outlier_calling_dm_output_dir, max_dm_junctions, jxn_filter_method, node_number, total_nodes, covariate_regression_method, covariate_file, rna_seq_samples_file, num_pc, outlier_calling_samples_file)
