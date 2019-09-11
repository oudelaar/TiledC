#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 09:26:32 2018

@author: oudelaar
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

###############################################################################

# Specify directories, file names and resolution

my_dir = "/path/" 
my_file = "myfile_5000_iced"

#5000 bp resolution
#chr11:29902951-33226736
#from coordinates.bed file:
#chr11	29900000	29904999	0
#chr11	33225000	33229999	665
#so 666 bins

n_bins = 666

###############################################################################

full_file_name = my_dir + my_file + ".matrix"
my_dir_out = "/path/"
full_file_name_out = my_dir_out + my_file 

###############################################################################

matrix = np.zeros((n_bins, n_bins))

with open(full_file_name) as f:
    for line in f:
        x1, x2, count = line.split()
        bin1 = int(x1)
        bin2 = int(x2)
        matrix[bin1, bin2] = float(count)
        matrix[bin2, bin1] = float(count)

mask = np.tri(matrix.shape[0], k = -1)
matrix_half = np.ma.array(matrix, mask = mask)         
threshold = np.percentile(matrix_half, 98)    
  
plot = plt.imshow(matrix_half, interpolation = "nearest", origin = "upper", vmin = 0.001, 
                  vmax = threshold, cmap=plt.cm.jet)
cbar = plt.colorbar()
plt.axis('off')  
plt.savefig(full_file_name_out + "_jet_" + str(threshold) + "_98.pdf", dpi=1000)
