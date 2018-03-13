import numpy as np
import os
import sys
import pdb
import gzip

def get_cluster_to_gene_mapping(file_name):
	f = open(file_name)
	dicti = {} #mapping from cluster_name to gene_name
	count = 0
	for line in f:
		if count == 0:
			count = count +1
			continue
		line = line.rstrip()
		data = line.split()
		header = data[0]
		header_data = header.split(':')
		cluster_name = header_data[3]
		genes = header_data[4].split(',')
		gene_dicti = {}
		for gene in genes:
			gene_dicti[gene] =1
		if cluster_name not in dicti:
			dicti[cluster_name] = gene_dicti
		else:
			for gene in genes:
				dicti[cluster_name][gene] = 1
	return dicti

def merge_output_files(total_nodes,tissue_specific_outlier_root,suffix,cluster_to_gene_mapping):
	genes_to_min_array = {}
	for node_num in range(total_nodes):
		file_name = tissue_specific_outlier_root + str(node_num) + suffix
		f = open(file_name)
		count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if count == 0:
				count = count + 1
				header = 'GENE_ID' + '\t' + '\t'.join(data[1:])
				continue
			cluster_id = data[0]
			gene_dicti = cluster_to_gene_mapping[cluster_id]
			counts = np.asarray(data[1:]).astype(float)
			for gene in gene_dicti.keys():
				if gene not in genes_to_min_array:
					genes_to_min_array[gene] = counts
				else:
					old_min = genes_to_min_array[gene]
					new_min = []
					for i,old_pval in enumerate(old_min):
						line_pval = counts[i]
						new_pval = np.nanmin([line_pval,old_pval])
						new_min.append(new_pval)
					genes_to_min_array[gene] = new_min
		f.close()
	t = open(tissue_specific_outlier_root + suffix,'w')
	t.write(header + '\n')
	for gene in genes_to_min_array.keys():
		t.write(gene + '\t' + '\t'.join(np.asarray(genes_to_min_array[gene]).astype(str)) + '\n')
	t.close()


tissue_type= sys.argv[1]
total_nodes= int(sys.argv[2])
tissue_specific_outlier_root= sys.argv[3]
tissue_specific_jxn_file=sys.argv[4]


cluster_to_gene_mapping = get_cluster_to_gene_mapping(tissue_specific_jxn_file)

merge_output_files(total_nodes,tissue_specific_outlier_root,'_emperical_pvalue.txt',cluster_to_gene_mapping)


#merge_output_files(total_nodes,tissue_specific_outlier_root,'_parametric_pvalue.txt',cluster_to_gene_mapping)
