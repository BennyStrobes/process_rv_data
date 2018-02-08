import numpy as np
import os
import sys
import pdb
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA

# Return array of individuals, corresponding to the order of individuals in the top line of the junction file
def get_ordered_individuals_from_junction_file(jxn_file):
    indi = []
    f = open(jxn_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        for ele in data[1:]:
            indi.append(ele)
        break  # Only care about header (where individuals are)
    return indi


# Extract matrix of covariates from v6p covariate file. Make sure individuals are in same order as 'ordered_individuals' array
def extract_non_pc_elements_of_v6p_cov_matrix(covariate_file_v6p, ordered_individuals):
    f = open(covariate_file_v6p)
    head_count = 0  # for header
    index = []  # Initialize index array that contains order of columns of the covariate_file_v6p so they match the order of ordered_individuals
    cov_output = []  # Initialize output array
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # header
            # In header, we will fill in index array so the ordered of columns int he covariate file match ordered_individuals
            head_indi = np.asarray(data[1:])  # array of individuals in header line
            head_count = head_count + 1
            for indi in ordered_individuals:  # Loop through ordered indiduals
                pos = np.where(indi == head_indi)[0][0]
                index.append(pos)

            if np.array_equal(head_indi[index], ordered_individuals) is False:  # Check to make sure re-ordering worked correctly
                print('erroror')
            continue
        covariate_name = data[0]  # Name of the covariate for this line
        covariate_array = np.asarray(data[1:])[index]  # Values of this covariate across all samples (ordered)
        ordered_line = np.hstack((covariate_name,covariate_array))
        if covariate_name == 'gender' or covariate_name == 'C1' or covariate_name == 'C2' or covariate_name == 'C3':  # The only covariates from that file we care about
            cov_output.append(ordered_line)
    cov_matrix = np.asmatrix(cov_output)  # Convert array of arrays into matrix form
    return cov_matrix


# Join all of the pieces together into a full covariate matrix
def merge_matrices(ordered_individuals, loadings_standardized, v6p_cov_matrix):
    # FORM HEADER
    header = np.hstack((['Individual'], ordered_individuals))
    row, col = loadings_standardized.shape

    # FORM JXN PC MATRIX
    new_arr = []
    for i in range(row):
        new_arr.append('Jxn_PC_' + str(i + 1))
    new_arr = np.asarray(new_arr)
    jxn_pc_mat = np.hstack((np.transpose(np.asmatrix(new_arr)), loadings_standardized))

    # FORM FINAL, matrix
    covariate_matrix = np.vstack((header, jxn_pc_mat, v6p_cov_matrix))
    return covariate_matrix

# Main driver function when covariate_regression == 'junction_pc'
# Run PCA on jxn matrix. Use top $num_pc principal component loadings as covariates. Also use non-PC covariates in covariate_file_v6p (race, gender, platform, etc)
def extract_junction_principal_components(tissue_type, covariate_output_file, jxn_file, covariate_file_v6p, num_pc, variance_explained_output_file):
    jxn_mat = np.transpose(np.loadtxt(jxn_file, dtype=str, delimiter='\t')[1:, 1:].astype(float))  # Extract junction matrix
    rows, cols = jxn_mat.shape
    jxn_mat_standardized = StandardScaler().fit_transform(jxn_mat)  # standardize each column of the junction matrix (PCA assumes gaussian distributed random variables)

    # Run PCA
    pca_object = sklearnPCA(n_components=20)  # Initialize sklearn pca object (20 was selected b/c it is large, and wanted to view large distribution of variance explained. Will filter to num_pc in following line)
    scores = pca_object.fit_transform(jxn_mat_standardized)[:, :num_pc]  # Run PCA
    U,S,V = np.linalg.svd(jxn_mat_standardized,full_matrices=False)
    loadings = np.transpose(U[:, :num_pc])
    variance_explained = pca_object.explained_variance_  # Extract array of amount of variance explained by each pc
    np.savetxt(variance_explained_output_file, variance_explained, delimiter='\n')

    # Print loadings to new output file
    ordered_individuals = get_ordered_individuals_from_junction_file(jxn_file)  # Return array of individuals, corresponding to the order of individuals in the top line of the junction file
    v6p_cov_matrix = extract_non_pc_elements_of_v6p_cov_matrix(covariate_file_v6p, ordered_individuals)  # Extract matrix of covariates from v6p covariate file. Make sure individuals are in same order as 'ordered_individuals' array
    new_covariate_matrix = merge_matrices(ordered_individuals, loadings, v6p_cov_matrix)  # Join all of the pieces together into a full covariate matrix

    np.savetxt(covariate_output_file, new_covariate_matrix, delimiter='\t', fmt="%s")


# Use non-PC covariates in covariate_file_v6p (race, gender, platform, etc), as well as top num_pc PCs from covariate file
def extract_v6p_eqtl_covariates(tissue_type, covariate_output_file, covariate_file_v6p, num_pc, jxn_file):
    ordered_individuals = get_ordered_individuals_from_junction_file(jxn_file)  # Return array of individuals, corresponding to the order of individuals in the top line of the junction file
    v6p_cov_matrix = extract_elements_of_v6p_cov_matrix(covariate_file_v6p, ordered_individuals, num_pc)  # Extract matrix of covariates from v6p covariate file. Make sure individuals are in same order as 'ordered_individuals' array

    # Merge matrices
    header = np.hstack((['Individual'], ordered_individuals))
    new_covariate_matrix = np.vstack((header, v6p_cov_matrix))
    np.savetxt(covariate_output_file, new_covariate_matrix, delimiter='\t', fmt="%s")


# Extract matrix of covariates from v6p covariate file. Make sure individuals are in same order as 'ordered_individuals' array
def extract_elements_of_v6p_cov_matrix(covariate_file_v6p, ordered_individuals, num_pc):
    f = open(covariate_file_v6p)
    head_count = 0  # for header
    index = []  # Initialize index array that contains order of columns of the covariate_file_v6p so they match the order of ordered_individuals
    cov_output = []  # Initialize output array
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # header
            # In header, we will fill in index array so the ordered of columns int he covariate file match ordered_individuals
            head_indi = np.asarray(data[1:])  # array of individuals in header line
            head_count = head_count + 1
            for indi in ordered_individuals:  # Loop through ordered indiduals
                pos = np.where(indi == head_indi)[0][0]
                index.append(pos)

            if np.array_equal(head_indi[index], ordered_individuals) is False:  # Check to make sure re-ordering worked correctly
                print('erroror')
            continue
        covariate_name = data[0]  # Name of the covariate for this line
        covariate_array = np.asarray(data[1:])[index]  # Values of this covariate across all samples (ordered)
        ordered_line = np.hstack((covariate_name,covariate_array))
        if covariate_name == 'gender' or covariate_name == 'C1' or covariate_name == 'C2' or covariate_name == 'C3':  # The only covariates from that file we care about
            cov_output.append(ordered_line)
        elif covariate_name.startswith('Inf'):
            if int(covariate_name.split('Cov')[1]) <= num_pc:
                cov_output.append(ordered_line)
    cov_matrix = np.asmatrix(cov_output)  # Convert array of arrays into matrix form
    return cov_matrix




# Command line args
tissue_type = sys.argv[1]
covariate_output_dir = sys.argv[2]  # output directory
covariate_directory_v6 = sys.argv[3]  # Input directory containing covariate files used in v6p eqtl analysis (one file for each tissue)
covariate_regression_method = sys.argv[4]  # What type of covariates do we want to include?
num_pc = int(sys.argv[5])  # Number of principal components to use
clusters_filter_output_dir = sys.argv[6]  # Directory containing filtered/processed junction files. To be used to get files and also if we are taking junction pcs

#File names
covariate_output_file = covariate_output_dir + tissue_type + '_' + covariate_regression_method + '_' + str(num_pc) + '_covariates.txt'  # Ouptut file
variance_explained_output_file = covariate_output_dir + tissue_type + '_' + covariate_regression_method + '_' + str(num_pc) + '_variance_explained.txt'  # Output file that will store variance explained for each junction pc (only used if covariate_regression_method == 'junction_pc')
jxn_file = clusters_filter_output_dir + tissue_type + '_hg19_filtered_xt_reclustered_gene_mapped.txt'  # Filtered/processsed junction file for this tissue
covariate_file_v6p = covariate_directory_v6 + tissue_type + '.covariates.txt'  # Covariate file used in v6p eqtl analysis for this tissue


if covariate_regression_method == 'junction_pc':  # Run PCA on jxn matrix. Use top $num_pc principal component loadings as covariates. Also use non-PC covariates in covariate_file_v6p (race, gender, platform, etc)
    extract_junction_principal_components(tissue_type, covariate_output_file, jxn_file, covariate_file_v6p, num_pc, variance_explained_output_file)
elif covariate_regression_method == 'v6p_eqtl_covariates':  # Use non-PC covariates in covariate_file_v6p (race, gender, platform, etc), as well as top num_pc PCs from covariate file
    extract_v6p_eqtl_covariates(tissue_type, covariate_output_file, covariate_file_v6p, num_pc, jxn_file)

