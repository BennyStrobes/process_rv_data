import sys
import pdb
import numpy as np


def extract_tissues(tissue_list_input_file):
    f = open(tissue_list_input_file)
    arr = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        arr.append(data[0])
    return arr


def merge_clusters(ss1, ss2, ss_to_clusters, clusters, cluster_number, header_line):
    cluster_name_old_1 = ss_to_clusters[ss1]  # Cluster name corresponding to old ss1
    cluster_name_old_2 = ss_to_clusters[ss2]  # Cluster name corresponding to old ss2
    cluster_name = 'cluster' + str(cluster_number)  # New cluster_id

    ss_to_clusters[ss1] = cluster_name  # Map ss1 to new_cluster_id
    ss_to_clusters[ss2] = cluster_name  # Map ss2 to new_cluster_id
    clusters[cluster_name] = []
    clusters[cluster_name].append(header_line)  # Map cluster_id to this jxn

    for old_header_line in clusters[cluster_name_old_1]:  # Remap all jxns in old cluster name corresponding to ss1 to new cluster id
        line_header_info = old_header_line.split(':')
        line_ss1 = line_header_info[0] + '_' + line_header_info[1]
        line_ss2 = line_header_info[0] + '_' + line_header_info[2]
        ss_to_clusters[line_ss1] = cluster_name
        ss_to_clusters[line_ss2] = cluster_name
        clusters[cluster_name].append(old_header_line)
    for old_header_line in clusters[cluster_name_old_2]:  # Remap all jxns in old cluster name corresponding to ss2 to new cluster id
        line_header_info = old_header_line.split(':')
        line_ss1 = line_header_info[0] + '_' + line_header_info[1]
        line_ss2 = line_header_info[0] + '_' + line_header_info[2]
        ss_to_clusters[line_ss1] = cluster_name
        ss_to_clusters[line_ss2] = cluster_name
        clusters[cluster_name].append(old_header_line)

    clusters.pop(cluster_name_old_1)  # Remove old cluster names
    clusters.pop(cluster_name_old_2)  # Remove old cluster names
    return ss_to_clusters, clusters


def initial_pass_for_cluster_ids(input_file, clusters, ss_to_clusters, cluster_number):
    f = open(input_file)
    head_count = 0

    for line in f:
        line = line.rstrip()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        data = line.split()

        header_line = data[0]  # First column contains jxn id info
        header_info = header_line.split(':')  # Convert jxn id info into array

        ss1 = header_info[0] + '_' + header_info[1]  # String corresponding to one of splice sites for this jxn
        ss2 = header_info[0] + '_' + header_info[2]  # String corresponding to other spice site for this jxn

        if ss1 not in ss_to_clusters and ss2 not in ss_to_clusters:  # Neither splice site has been seen before --> create new cluster_id
            cluster_name = 'cluster' + str(cluster_number)
            ss_to_clusters[ss1] = cluster_name
            ss_to_clusters[ss2] = cluster_name
            clusters[cluster_name] = []
            clusters[cluster_name].append(header_line)
            cluster_number = cluster_number + 1
        elif ss1 not in ss_to_clusters and ss2 in ss_to_clusters:  # ss2 has been seen before, but ss1 has not been seen before. Map this jxn to the cluster ss2 is mapped to
            cluster_name = ss_to_clusters[ss2]
            ss_to_clusters[ss1] = cluster_name
            ss_to_clusters[ss2] = cluster_name
            clusters[cluster_name].append(header_line)
        elif ss1 in ss_to_clusters and ss2 not in ss_to_clusters:  # ss1 has been seen before, but ss2 has not been seen before. Map this jxn to the cluster ss1 is mapped to.
            cluster_name = ss_to_clusters[ss1]
            ss_to_clusters[ss1] = cluster_name
            ss_to_clusters[ss2] = cluster_name
            clusters[cluster_name].append(header_line)
        elif ss1 in ss_to_clusters and ss2 in ss_to_clusters:  # ss1 and ss2 have been seen before (most interesting case)
            cluster_name1 = ss_to_clusters[ss1]
            cluster_name2 = ss_to_clusters[ss2]
            if cluster_name1 != cluster_name2:  # The cluster ss1 is previously mapped to and the cluster ss2 is previously mapped to disagree. --> Merge clusters to a new one and delete old
                ss_to_clusters, clusters = merge_clusters(ss1, ss2, ss_to_clusters, clusters, cluster_number, header_line)
                cluster_number = cluster_number + 1
            else:  # ss1 and ss2 previously mapped to the same cluster.
                cluster_name = cluster_name1
                ss_to_clusters[ss1] = cluster_name
                ss_to_clusters[ss2] = cluster_name
                clusters[cluster_name].append(header_line)

    return clusters, ss_to_clusters, cluster_number


def final_pass_and_print(clusters, ss_to_clusters, input_file, output_file):
    '''
    Take first pass to update tissue_cluster_counts (ie keep track of how many times a cluster is called in this tissue)
    '''
    f = open(input_file)
    t = open(output_file, 'w')
    head_count = 0

    tissue_cluster_counts = {}  # Map from cluster_id to how many jxns that cluster_id contains in this tissue

    for line in f:
        line = line.rstrip()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        data = line.split()
        header_line = data[0]  # First column of jxn file contains info on jxn
        header_info = header_line.split(':')  # Split info into array

        ss1 = header_info[0] + '_' + header_info[1]  # First splice site for jxn
        ss2 = header_info[0] + '_' + header_info[2]  # Second splice site for jxn

        cluster_name1 = ss_to_clusters[ss1]  # Cluster_id corresponding to ss1
        cluster_name2 = ss_to_clusters[ss2]  # Cluster_id corresponding to ss2
        if cluster_name1 != cluster_name2:  # Make sure cluster_name1 == cluster_name2 --> else there is an error
            print('EROROOROR')
            pdb.set_trace()
        if len(np.unique(clusters[cluster_name1])) == 1:  # Filter out jxns that belong to a cluster with only one jxn (across tissues)
            continue
        # Update tissue_cluster_counts
        if cluster_name1 not in tissue_cluster_counts:  # Cluster_name has yet to be seen in this tissue. Initialize key
            tissue_cluster_counts[cluster_name1] = 1
        else:  # Cluster_name has been seen before in this tissue. Add count
            tissue_cluster_counts[cluster_name1] = tissue_cluster_counts[cluster_name1] + 1
    f.close()
    '''
    Take second pass to filter any jxns that have tissue_cluster_counts == 1
    '''
    f = open(input_file)
    head_count = 0

    for line in f:
        line = line.rstrip()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            t.write(line + '\n')
            continue
        data = line.split()
        header_line = data[0]  # First column of jxn file contains info on jxn
        header_info = header_line.split(':')  # Split into array
        ss1 = header_info[0] + '_' + header_info[1]
        ss2 = header_info[0] + '_' + header_info[2]
        cluster_name1 = ss_to_clusters[ss1]
        cluster_name2 = ss_to_clusters[ss2]
        if cluster_name1 != cluster_name2:
            print('EROROOROR')
            pdb.set_trace()
        if len(np.unique(clusters[cluster_name1])) == 1:
            continue
        if tissue_cluster_counts[cluster_name1] == 1:  # Filter out any jxn that has no other jxns mapped to this cluster in this tissue
            continue
        new_header_line = header_info[0] + ':' + header_info[1] + ':' + header_info[2] + ':' + cluster_name1  # Jxnid string for output
        t.write(new_header_line + '\t' + '\t'.join(data[1:]) + '\n')
    t.close()


def run_analysis(tissues, clusters_filter_output_dir, input_suffix, output_suffix):
    clusters = {}  # Initialize object that maps cluster id to a list of all jxns assigned to that cluster id
    ss_to_clusters = {}  # Initialize object that maps splice site to cluster_id
    cluster_number = 0  # initialize integer variable that keeps track of the next cluster_id name to be used

    for tissue in tissues:
        input_file = clusters_filter_output_dir + tissue + input_suffix  # Input jxn file
        clusters, ss_to_clusters, cluster_number = initial_pass_for_cluster_ids(input_file, clusters, ss_to_clusters, cluster_number)  # Update objects 'clusters', 'ss_to_clusters', and 'cluster_number' for this specific tissue

    for tissue in tissues:
        input_file = clusters_filter_output_dir + tissue + input_suffix  # Input jxn file
        output_file = clusters_filter_output_dir + tissue + output_suffix  # Ouput jxn file
        final_pass_and_print(clusters, ss_to_clusters, input_file, output_file)  # Run through jxn file one more time. To make sure a jxn is only printed if it has more than one jxn mapped to its cluster IN THIS TISSUE ALONE

tissue_list_input_file = sys.argv[1]
clusters_filter_output_dir = sys.argv[2]  # Input dir and output dir


tissues = extract_tissues(tissue_list_input_file)  # get array of tissue types


input_suffix = '_hg19_filtered.txt'  # Suffix of input files
output_suffix = '_hg19_filtered_xt_reclustered.txt'  # Suffix of outputfiles

run_analysis(tissues, clusters_filter_output_dir, input_suffix, output_suffix)


'''
ERROR CHECKING RUN:
1. A given splice site always maps to the same cluster across tissues
'''
