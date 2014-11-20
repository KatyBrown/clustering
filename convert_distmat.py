# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 14:14:02 2014

@author: katherineb
"""
import numpy as np
import pandas as pd
def convert_distmat(data):
    '''
    reads a distance matrix in "long" format and converts to wide format
    input format:
        a   b   d
        a   c   d
        b   c   d
   
    '''
    mat = pd.read_table(data, delimiter="\s", names=['id1', 'id2', 'dist'])
    mat = pd.pivot(index = mat['id1'], columns = mat['id2'], values = mat['dist'])
    mat.update(mat.transpose())
    if mat.index.all() != mat.columns.all():
        i = list(set(mat.index).difference(set(mat.columns)))[0]
        j = list(set(mat.columns).difference(set(mat.index)))[0]
        print j
        mat.loc[:,i] = mat.loc[i,:] 
        mat.loc[j,:] = mat.loc[:,j]
    mat = mat.replace("NaN", float(0))    
    mat = mat.sort(axis=0)
    mat = mat.sort(axis=1)
    mat2 = mat.values
    ids = list(mat.index)
    return mat2, ids