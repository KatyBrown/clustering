# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 13:04:21 2014

@author: katherineb
"""

import numpy as np
import pandas as pd
from convert_distmat import convert_distmat
def read_distmat(filename):
    '''
    reads a distance matrix from filename
    if distances are a list of pairwise comparisons rather than a matrix, 
    converts to matrix using convert_distmat
    '''
    ids = open(filename).readline().strip().split("\t")
    if len(ids) == 3:
        result = convert_distmat(filename)
        mat = result[0]
        ids = result[1]
    else:
        mat = np.genfromtxt(filename, delimiter = "\t", skip_header=1, usecols = range(1,len(ids)+1))
    return mat, ids


