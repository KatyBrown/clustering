# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 15:31:36 2014

@author: katherineb
"""
import numpy as np
import pandas as pd
import math

class density_peak_clusters(object):
    def __init__(self, data, outfile):
        self.matrix = data[0]
        self.ids = data[1]
        self.outfile = outfile
        self.dc = np.percentile(self.matrix, 2) # dc  = cutoff distance
        self.pis = self.gaussian_k()
        if self.pis != None: # if > 0 cluster centres are identified
            self.deltas = self.delta_all()
            self.nnhds = self.nnhd_all()
            self.ccs = self.cluster_centres()
            self.ccs = self.assign()
            self.dens = self.near_other_clusters()
            self.inbs = self.in_border()
        else:
            self.inbs = None
              
    def get_clusters(self):
        '''
        returns a dataframe of each id and the cluster that the point with that
        id falls into (or 0 if the point is not in a cluster)
        '''
        if self.inbs != None:
            double = np.vstack((self.ids, self.inbs)).transpose()
            df = pd.DataFrame(data=double)
            df.to_csv(self.outfile, sep= "\t", header=False, index=False)
            return df
        else:
            z = np.zeros(len(self.ids), dtype=int)
            double = np.vstack((self.ids, z)).transpose()
            df = pd.DataFrame(data=double)
            df.to_csv(self.outfile, sep= "\t", header=False, index=False)
            return df
        
    def gaussian_k(self):
        """
        estimates pi for all points
        pi = density = number of points that are closer than dc to a point
        returns self.pis
        """
        ids = self.ids
        matrix = self.matrix
        pis = np.zeros(len(ids))
        for id1 in range(len(ids)):
            for id2 in range(id1+1, len(ids)):
                dist = matrix[(id1, id2)]
                pis[id1] = pis[id1] + math.exp(-(dist/self.dc)*(dist/self.dc))
                pis[id2] = pis[id2] + math.exp(-(dist/self.dc)*(dist/self.dc))
        if len(pis[pis > max(pis) * 0.95]) > len(pis[pis < max(pis) * 0.1]):
            return None
        return pis
        
    def calc_delta(self, row, pi):
        """
        calculates delta for a single row
        delta = minimum distance between the point and any other point with 
        higher density
        """
        try:
            return min(row[self.pis > pi])
        except:
            return max(row)
            
    def delta_all(self):
        """
        runs calc_delta on all data points
        generates self.deltas - an ordered list of deltas
        """
        i = 0
        deltas = list()
        for row in self.matrix:
            pi = self.pis[i]
            delta = self.calc_delta(row, pi)
            deltas.append(delta)
            i += 1
        return np.array(deltas)
    
    def nnhd_all(self):
        """
        calculates nnhd for all rows
        nnhd = nearest neighbour of higher density
        generates nnhds - an ordered list of nnhds
        """
        i = 0
        nnhds = list()
        for row in self.matrix:
            delta = self.deltas[i]
            nnhds.append(int(np.where(row == delta)[0]))
            i += 1
        return np.array(nnhds)
    
    def cluster_centres(self):
        ''' 
        identifies the cluster centres - points with high pi and high delta
        returns ccs - the group assignations of the cluster centres with 0s for
        every other point    
        '''
        pis = self.pis
        deltas = self.deltas
        clusters = np.zeros(len(pis), dtype=int)
        pmax = max(pis) * 0.3
        dmax = max(deltas) * 0.15
        inds = np.intersect1d(np.where(pis > pmax)[0], np.where(deltas > dmax)[0]) # cluster centres are points within the top right of the decision graph
        p = 1
        for i in inds:
            clusters[i] = p
            p += 1
        return clusters
    
    def assign(self):
        '''
        assigns the remaining points to clusters based on the cluster of their nnhd
        returns ccs - an ordered list of cluster assignations
        '''
        ccs = self.ccs
        nnhds = self.nnhds
        j = 0
        while True:
            j += 1
            i = 0
            for cc in ccs:
                if cc == 0:
                    nnhd = nnhds[i]
                    if ccs[nnhd] != 0:
                        ccs[i] = ccs[nnhd]
                i += 1
            if j > len(nnhds):
                print "clustering failed" # if more loops than points in the dataset
                break
            if 0 not in ccs: # when all the points are assigned to a cluster, break
                return ccs
        
    def near_other_clusters(self):
        '''
        determines whether each point is within dc of a point in another cluster
        for each cluster returns the border density
        border density = the highest density of a point within dc of a point in 
        another cluster
        '''
        i = 0
        ccs = self.ccs
        dens = [0] * len(list(set(ccs)))
        for row in self.matrix:
            x = np.where(row < self.dc)[0]
            y = np.where(self.ccs == self.ccs[i])[0]     
            z = np.setdiff1d(x,y) # finds the overlap between these
            if len(z) != 0:
                if dens[(self.ccs[i] - 1)] < self.pis[i]:
                    dens[(self.ccs[i] - 1)] = self.pis[i]
            i += 1
        return np.array(dens)

    def in_border(self):
        ''' 
        finds and returns the points with a density below the border density 
        for their cluster
        '''
        inbs = []
        for i in range(len(self.pis)):
            cluster_dens = self.dens[self.ccs[i] - 1]
            if self.pis[i] > cluster_dens:
                inbs.append(self.ccs[i])
            else:
                inbs.append(0)
            i += 1
        return np.array(inbs)