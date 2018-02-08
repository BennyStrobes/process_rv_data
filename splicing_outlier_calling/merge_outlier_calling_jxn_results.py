import numpy as np
import os
import sys
import pdb


def go(tissue_type, total_nodes, tissue_specific_outlier_root, version):
    output_file = tissue_specific_outlier_root + '_jxn_level_' + version + '.txt'
    t = open(output_file, 'w')
    for node_number in range(total_nodes):
        input_file = tissue_specific_outlier_root + str(node_number) + '_' + version
        f = open(input_file)
        count = 0
        for line in f:
            line = line.rstrip()
            if count == 0:
                count = count + 1
                if node_number == 0:
                    t.write(line + '\n')
                continue
            t.write(line + '\n')
        f.close()
    t.close()


tissue_type = sys.argv[1]
total_nodes = int(sys.argv[2])
tissue_specific_outlier_root = sys.argv[3]

go(tissue_type, total_nodes, tissue_specific_outlier_root, 'emperical_pvalue.txt')
