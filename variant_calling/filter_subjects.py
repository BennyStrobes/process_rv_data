import numpy as np
import os
import sys
import pdb




def get_rna_seq_individuals(rna_seq_individuals_file):
    f = open(rna_seq_individuals_file)
    individ = {}
    for line in f:
        line = line.rstrip()
        individ[line] = 1
    return individ

def filter_european_individuals(subject_phenotype_file, rna_seq_indi):
    final_indi = {}
    f = open(subject_phenotype_file)
    count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if count == 0:
            count = count + 1
            continue
        race = data[4]
        line_indi = data[0]
        if race != '3':
            continue
        if line_indi not in rna_seq_indi:
            continue
        final_indi[line_indi] = 1
    print('There are ' + str(len(final_indi)) + ' individuals that we have rna seq measurements for AND are of european ancestry')
    return final_indi


def print_helper(final_individuals, output_file):
    t = open(output_file, 'w')
    for indi in sorted(final_individuals.keys()):
        t.write(indi + '\n')
    t.close()

subject_phenotype_file = sys.argv[1]
rna_seq_individuals_file = sys.argv[2]
output_file = sys.argv[3]

rna_seq_individuals = get_rna_seq_individuals(rna_seq_individuals_file)
final_individuals = filter_european_individuals(subject_phenotype_file, rna_seq_individuals)
print_helper(final_individuals, output_file)