#!/usr/bin/env python

import numpy as np
import nibabel as nib
from scipy.sparse import coo_matrix, cs_graph_components

# Borrowed from nipy.graph.graph
# https://github.com/nipy/nipy/blob/master/nipy/algorithms/graph/graph.py
def graph_3d_grid(xyz, k=18):
    """ Utility that computes the six neighbors on a 3d grid

    Parameters
    ----------
    xyz: array of shape (n_samples, 3); grid coordinates of the points
    k: neighboring system, equal to 6, 18, or 26

    Returns
    -------
    i, j, d 3 arrays of shape (E),
            where E is the number of edges in the resulting graph
            (i, j) represent the edges, d their weights
    """
    if np.size(xyz) == 0:
        return None
    lxyz = xyz - xyz.min(0)
    m = 3 * lxyz.max(0).sum() + 2

    # six neighbours
    n6 = [np.array([1, m, m ** 2]), np.array([m ** 2, 1, m]),
         np.array([m, m ** 2, 1])]

    # eighteen neighbours
    n18 = [np.array([1 + m, 1 - m, m ** 2]),
           np.array([1 + m, m - 1, m ** 2]),
           np.array([m ** 2, 1 + m, 1 - m]),
           np.array([m ** 2, 1 + m, m - 1]),
           np.array([1 - m, m ** 2, 1 + m]),
           np.array([m - 1, m ** 2, 1 + m])]

    # twenty-six neighbours
    n26 = [np.array([1 + m + m ** 2, 1 - m, 1 - m ** 2]),
           np.array([1 + m + m ** 2, m - 1, 1 - m ** 2]),
           np.array([1 + m + m ** 2, 1 - m, m ** 2 - 1]),
           np.array([1 + m + m ** 2, m - 1, m ** 2 - 1])]

    # compute the edges in each possible direction
    def create_edges(lxyz, nn, l1dist=1, left=np.array([]), right=np.array([]),
                     weights=np.array([])):
        q = 0
        for nn_row in nn:
            v1 = np.dot(lxyz, nn_row)
            o1 = np.argsort(v1)
            sv1 = v1[o1]
            nz = np.squeeze(np.nonzero(sv1[: - 1] - sv1[1:] == - l1dist))
            o1z, o1z1 = o1[nz], o1[nz + 1]
            left = np.hstack((left, o1z, o1z1))
            right = np.hstack((right, o1z1, o1z))
            q += 2 * np.size(nz)
        weights = np.hstack((weights, np.sqrt(l1dist) * np.ones(q)))
        return left, right, weights

    i, j, d = create_edges(lxyz, n6, 1.)
    if k >= 18:
        i, j, d = create_edges(lxyz, n18, 2, i, j, d)
    if k == 26:
        i, j, d = create_edges(lxyz, n26, 3, i, j, d)
    i, j = i.astype(np.int), j.astype(np.int)

    # reorder the edges to have a more standard order
    order = np.argsort(i + j * (len(i) + 1))
    i, j, d = i[order], j[order], d[order]
    return i, j, d


def fstat_given_pval(fstats, pvals, given_p):
    A       = np.ones((len(fstats),2))
    A[:,0]  = pvals
    m,c     = np.linalg.lstsq(A, fstats)[0]
    ret_f   = m*given_p + c
    return ret_f

def cluster_correct_perms(Fperms, fstats, pvals, mask, clust_thr=0.05, vox_thr=0.05):
    """
    Cluster correct the MDMR output using a permutation-based test.
    
    Note that the voxel-level threshold is converted from a p-value into an F-statistic.
    This is done by regressing the p-values onto the f-stats, and then predicting the 
    f-statistic based on the desired p-value. This f-statistic is then used as the threshold
    throughout.
        
    Parameters
    ----------
    Fperms : list of strings
        A length `N` list of file paths of the nifti files of subjects
    mask_file : string
        Path to a mask file in nifti format (this will serve as a prior mask)
    dtype : string
        Numpy data type, should be one of 'float16', 'float32', or 'float64'
    
    Returns
    -------
    joint_mask : string
        Path to joint mask file in nifti format
        (this will be the current directory followed by join_mask.nii.gz)
    
    
    
    Return a table with indices, cluster sizes, cluster pvals
    along with the voxelwise cluster data.
    Note that the table will be for the cluster corrected output.
    """
    
    ###
    # Setup
    
    # Coordinates within mask
    xyz     = np.where(mask > 0)
    xyz     = np.array(xyz).T
    
    # Determine f-stat threshold given a p-value
    f_thr   = fstat_given_pval(fstats, pvals, vox_thr)
    
    # Only retain coordinates greater than threshold
    xyz_th  = xyz[fstats > f_thr]
    
    # Get clusters above threshold
    clust   = clusterize(xyz_th)
    
    ###
    
    
    ###
    # Cluster Info
    
    # Label
    labels  = np.array([ l for l in np.unique(clust) if l > -2 ])
    
    # Size
    clust_sizes = np.array([ (clust==l).sum() for l in labels ])
    
    ###
    
    
    ###
    # Cluster each permutation and get maximum cluster size
    
    def cluster_perm_data(perm_fs, f_thr):
        xyz_th  = xyz[perm_fs > f_thr]
        labels  = clusterize(xyz_th)
        return max_clust_size(labels)
    
    max_sizes = np.apply_along_axis(cluster_perm_data, 1, Fperms, f_thr)
    
    ###
    
    
    ###
    # Threshold clusters of significance
    
    # Get each clusters p-values
    clust_pvals = np.array([ (cz>=max_sizes).mean() for cz in clust_sizes ])
    
    # Threshold
    subset      = clust_pvals < clust_thr
    labels_thr  = labels[subset]
    clust_thr   = np.zeros_like(clust)
    for i,l in enumerate(labels_thr):
        clust_thr[clust==l] = i+1
    
    ###
    
    
    ###
    # Summarize corrected clusters
    
    nclust      = subset.sum()
    clust_inds  = np.array(range(nclust)) + 1
    table       = np.hstack((clust_inds, clust_sizes[subset], clust_pvals[subset]))
    
    ###
    
    return (table, clust_thr)


def cluster_file(img_file, thr, mask_file):
    mask = nib.load(mask_file).get_data()
    img  = nib.load(img_file).get_data()
    return cluster_data(img, thr, mask)

def cluster_data(img, thr, mask):
    """docstring for cluster_data"""
    xyz     = np.where(mask > 0)
    xyz_a   = np.array(xyz).T
    xyz_th  = xyz_a[img[xyz] > thr]
    labels  = clusterize(xyz_th)
    return labels

def clusterize(xyz_th, k=26):
    """docstring for clusterize"""
    nvoxs       = xyz_th.shape[0]
    i,j,d       = graph_3d_grid(xyz_th, k=k)    
    adj         = coo_matrix((d, (i,j)), shape=(nvoxs,nvoxs))
    nc, labels  = cs_graph_components(adj)
    return labels

def cluster_sizes(labels):
    return np.array([ (labels==l).sum() for l in np.unique(labels) if l > -2 ])

def max_clust_size(labels):
    return cluster_sizes(labels).max()

def fsl_cluster_sizes(img_file, thr):
    from subprocess import Popen, PIPE
    cmd = "cluster -i %s -t %f | awk '{print $2}' | grep -v '^Index'" % (img_file, thr)
    p   = Popen(cmd, shell=True, stdout=PIPE)
    stdout, _ = p.communicate()
    clust_sizes = np.array([ int(i) for i in stdout.split("\n") if i != '' ])
    return clust_sizes
