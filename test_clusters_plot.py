# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 12:59:51 2014

@author: katherineb
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 13:39:45 2014

@author: katherineb
"""
from sklearn import datasets
import math
import matplotlib.pyplot as plt
from read_distmat import read_distmat
from class_density_cluster import density_peak_clusters


def points(x, y):
    """
    makes a list of coordinates of positions of points for all
    indices in x and y
    """

    ps = []
    for i in range(len(x)):
        p1 = x[i]
        p2 = y[i]
        p = (p1, p2)
        ps.append(p)

    return ps


def make_dist_mat(p, output):
    """
    calculates the distance between pairs of points
    writes the distance matrix to file output
    """
    out = open(output, "w")
    for i in range(len(p)):
        for j in range(len(p)):
            x1 = p[i][0]
            x2 = p[j][0]
            y1 = p[i][1]
            y2 = p[j][1]
            dist = math.hypot(x2-x1, y2-y1)
            out.write(str(i) + "\t" + str(j) + "\t" + str(dist) + "\n")
    out.close()


def output_raw(data, output):
    """ writes the raw data to file output"""
    out = open(output, "w")
    for line in data:
        out.write(str(line[0])+"\t"+str(line[1])+"\n")
    out.close()


def pointcol(clusters, points, z):

    plt.subplot(4, 2, z)
    cols = plt.cm.spectral
    selectcols = np.linspace(0.05, 0.95, max(clusters))
    colours = [(0.7, 0.7, 0.7, 1.0)]
    for i in selectcols:
        print i
        print cols(i)
        colours.append(cols(i))
    i = 0
    for item in points:
        plt.plot(points[i][0],
                 points[i][1],
                 markerfacecolor=colours[clusters[i]],
                 markeredgecolor=colours[clusters[i]],
                 marker='o',
                 lw=0)
        i += 1

# random clustered datasets with different parameters
data1 = datasets.make_blobs(centers=8, n_samples=300)[0]
data2 = datasets.make_blobs(centers=3, n_samples=100)[0]
data4 = datasets.make_blobs(centers=2, n_samples=1000)[0]
data5 = datasets.make_moons(n_samples=100, noise=.05)[0]
data6 = np.vstack(((np.random.rand(100, 2) * 20) - 10, data1))
data7 = ((np.random.rand(300, 2) * 20) - 10)
data8 = datasets.make_circles(n_samples=500, noise=.05)[0]
datasets = [data1, data2, data3, data4, data5, data6, data7, data8]

i = 1
plt.figure(figsize=(20, 20))
for data in datasets:
    print "dataset " + str(i)
    try:
        X = list(data[:,0])
        Y = list(data[:,1])
    except:
        X = list(data[0])
        Y = list(data[1])
    output_raw(data, "raw_" + str(i) + ".out")
    p = points(X, Y)
    dist = make_dist_mat(p, "dist_" + str(i) + ".out")
    data = read_distmat("dist_" + str(i) + ".out")
    outfile = "ccs_" + str(i) + ".out"
    a = density_peak_clusters(data, outfile).get_clusters()
    clusters = a[1]
    pointcol(clusters, p, i)
    i += 1
