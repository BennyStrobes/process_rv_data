import sys
import pdb
import os
import numpy as np 
import time
import gzip
from scipy import stats
import scipy.stats as ss

def go(input_file,output_file):
	f = open(input_file)
	t = open(output_file,'w')
	count = 0
	for line in f:
		line = line.rstrip()
		data =line.split()
		if count == 0:
			count = count + 1
			continue
		info = data[0].split(':')
		cluster_id = info[0] + ':' + info[1] + ':' + info[2]
		counts = np.asarray(data[1:]).astype(float)
		number = len(np.where(counts > 0)[0])
		t.write(cluster_id + '\t' + str(number) + '\n')
	t.close()

input_file = sys.argv[1]
output_file = sys.argv[2]


go(input_file,output_file)
