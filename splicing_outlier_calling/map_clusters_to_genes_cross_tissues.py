import numpy as np
import os
import sys
import pdb
import gzip

def extract_tissues(tissue_list_input_file):
    f = open(tissue_list_input_file)
    arr = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        arr.append(data[0])
    return arr


#Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
def add_gene_to_chromosome_object(chromosome, start, end, gene_name):
    for pos in range(start, end+1):
        if chromosome[pos] == 'NULL':  # No genes already mapped to this position
            chromosome[pos] = gene_name
        else:  # at least one gene already mapped to this position
            chromosome[pos] = chromosome[pos] + ',' + gene_name
    return chromosome


#Create an array of lenth(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
def make_chromosome_with_gene_names(chrom_num, gencode_hg19_gene_annotation_file):
    #initialize chromosome array
    chromosome = ['NULL']*259250621
    #loop through gencode file
    f = gzip.open(gencode_hg19_gene_annotation_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#'):  # ignore header lines
            continue
        gene_type = data[13].split('"')[1]  # ie protein_coding,pseudo_gene,etc
        gene_name = data[9].split('"')[1]  # ensamble id
        line_chrom_num = data[0]
        gene_part = data[2]  # gene,UTR,exon,etc
        if gene_type != 'protein_coding':  # limit to only protein coding genes
            continue
        if line_chrom_num != 'chr' + chrom_num:  # limit to chromosome of interest
            continue
        if gene_part != 'gene':  # ignore other parts of the gene as 'gene' encomposes everything
            continue
        start = int(data[3])  # 5' gene start site (this is the min of all UTR,exons,etc)
        end = int(data[4])  # 3' gene end site (this is the max of all UTR,exons,etc)
        #Fill in chromosome object from chromosome[start:end+1] with the addition of this gene name
        chromosome = add_gene_to_chromosome_object(chromosome, start, end, gene_name)
        #NOTE: ENSAMBLE Id's never come up more than once when gene_part == 'gene'
    return chromosome


def get_unique_genes(stringer):
    info = stringer.split(',')
    return ','.join(np.unique(info))


#Loop through jxn file. For each jxn, see what genes overlap is 5' ss and 3' ss. And the union of the genes to cluster_to_genes_mapping for the jxns associated cluster
def align_clusters_to_genes(chromosome, chrom_num, cluster_to_genes_mapping, input_file_name):
    f = open(input_file_name)
    count = 0  # header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if count == 0:  # ignore header
            count = count + 1
            continue
        junction_info = data[0].split(':')
        line_chrom_num = junction_info[0]
        if line_chrom_num != 'chr' + chrom_num:
            continue
        start = int(junction_info[1])
        end = int(junction_info[2])
        cluster_id = junction_info[3]
        genes_at_start = chromosome[start]
        genes_at_end = chromosome[end]
        if genes_at_start == 'NULL' and genes_at_end == 'NULL':  # Both ends of the jxn overlap zero genes
            continue
        elif genes_at_start == 'NULL':
            gene_string = genes_at_end
        elif genes_at_end == 'NULL':
            gene_string = genes_at_start
        else:
            gene_string = genes_at_start + ',' + genes_at_end
        if cluster_id not in cluster_to_genes_mapping:
            cluster_to_genes_mapping[cluster_id] = get_unique_genes(gene_string)
        else:
            cluster_to_genes_mapping[cluster_id] = get_unique_genes(gene_string + ',' + cluster_to_genes_mapping[cluster_id])
    return cluster_to_genes_mapping


#Perform mapping analysis on each chromosome seperately
#At each chromosome update the cluster_to_genes_mapping for clusters on that chromosome
def map_clusters_on_chromosome(chrom_num, cluster_to_genes_mapping, clusters_filter_output_dir, input_suffix, tissues, gencode_hg19_gene_annotation_file):
    #Create an array of lenth(chromosome) [in bp]. Value of an element corresponds to a list of genes that overlap that basepair. This is used for efficient searching.
    chromosome = make_chromosome_with_gene_names(chrom_num, gencode_hg19_gene_annotation_file)
    #Loop through jxn file. For each jxn, see what genes overlap is 5' ss and 3' ss. And the union of the genes to cluster_to_genes_mapping for the jxns associated cluster
    for tissue in tissues:
        input_file_name = clusters_filter_output_dir + tissue + input_suffix
        cluster_to_genes_mapping = align_clusters_to_genes(chromosome, chrom_num, cluster_to_genes_mapping, input_file_name)
    return cluster_to_genes_mapping


def print_helper_map_clusters_to_genes(input_file_name, output_file_name, cluster_to_genes_mapping):
    f = open(input_file_name)
    t = open(output_file_name, 'w')
    count = 0  # header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if count == 0:
            count = count + 1
            t.write(line + '\n')
            continue
        junction_info = data[0].split(':')
        cluster_id = junction_info[3]
        if cluster_id not in cluster_to_genes_mapping:
            continue
        new_jxn_id = data[0] + ':' + cluster_to_genes_mapping[cluster_id]
        t.write(new_jxn_id + '\t' + '\t'.join(data[1:]) + '\n')
    t.close()


def map_clusters_to_genes(tissues, gencode_hg19_file, clusters_filter_output_dir, input_suffix, output_suffix):
    #Main dictionary to keep track of cluster --> genes.
    #Key will be leafcutter cluster id. Value is string gene1,gene2,gene3,... Or 'NULL' if no gene is mapped
    cluster_to_genes_mapping = {}
    #Perform mapping analysis on each chromosome seperately
    #At each chromosome update the cluster_to_genes_mapping for clusters on that chromosome
    for chrom_num in range(1, 23):
        cluster_to_genes_mapping = map_clusters_on_chromosome(str(chrom_num), cluster_to_genes_mapping, clusters_filter_output_dir, input_suffix, tissues, gencode_hg19_file)

    for tissue in tissues:
        input_file = clusters_filter_output_dir + tissue + input_suffix
        output_file = clusters_filter_output_dir + tissue + output_suffix
        print_helper_map_clusters_to_genes(input_file, output_file, cluster_to_genes_mapping)


def report_cluster_info(tissues, clusters_filter_output_dir, output_suffix, output_file):
    clusters = {}
    for tissue in tissues:
        file_name = clusters_filter_output_dir + tissue + output_suffix
        f = open(file_name)
        count = 0
        for line in f:
            if count == 0:
                count = count + 1
                continue
            line = line.rstrip()
            data = line.split()

            junction_info = data[0].split(':')
            cluster_id = junction_info[3]
            jxn_name = junction_info[0] + ':' + junction_info[1] + ':' + junction_info[2]
            genes = junction_info[4].split(',')
            if cluster_id not in clusters:
                clusters[cluster_id] = {}
                clusters[cluster_id]['jxns'] = []
                clusters[cluster_id]['jxns'].append(jxn_name)
                clusters[cluster_id]['genes'] = []
                for gene in genes:
                    clusters[cluster_id]['genes'].append(gene)
            else:
                clusters[cluster_id]['jxns'].append(jxn_name)
                for gene in genes:
                    clusters[cluster_id]['genes'].append(gene)
    t = open(output_file, 'w')
    t.write('cluster_id\tjxns\tgenes\n')
    for cluster in sorted(clusters.keys()):
        genes = np.unique(clusters[cluster]['genes'])
        jxns = np.unique(clusters[cluster]['jxns'])
        t.write(cluster + '\t' + ','.join(jxns) + '\t' + ','.join(genes) + '\n')


tissue_list_input_file = sys.argv[1]
clusters_filter_output_dir = sys.argv[2]  # Input dir and output dir
gencode_hg19_file = sys.argv[3]

tissues = extract_tissues(tissue_list_input_file)  # get array of tissue types

input_suffix = '_hg19_filtered_xt_reclustered.txt'  # Suffix of input files
output_suffix = '_hg19_filtered_xt_reclustered_gene_mapped.txt'  # Suffix of output files

map_clusters_to_genes(tissues, gencode_hg19_file, clusters_filter_output_dir, input_suffix, output_suffix)

# Make output file containing all clusters as well as which junctions are mapped to which cluster, as well as which genes are mapped to which cluster
report_cluster_info(tissues, clusters_filter_output_dir, output_suffix, clusters_filter_output_dir + 'cluster_info.txt')
