import os
import time
from os.path import expanduser
import matplotlib
import numpy as np
import nibabel as nb
import pandas as pd
import nilearn.image as image
import scipy as sp
import imp
from os.path import expanduser


from nilearn import datasets
from nilearn.input_data import NiftiMasker
from nilearn.plotting import plot_roi, show
from nilearn.image.image import mean_img
from nilearn.image import resample_img
from matplotlib import pyplot as plt


from sklearn import cluster, datasets
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler, normalize

home2 = expanduser("~")


#CPAC = imp.load_source('CPAC', home2 + '/C-PAC')
##bascutils = imp.load_source('bascutils', home2 + '/C-PAC/CPAC/basc/utils.py')
#basc = imp.load_source('basc', home2 + '/C-PAC/CPAC/basc/basc.py')
#utilsutils = imp.load_source('utilsutils', home2 + '/C-PAC/CPAC/utils/utils.py')

#import CPAC
#from CPAC.basc.utils import standard_bootstrap, adjacency_matrix, cluster_timeseries, cluster_matrix_average, individual_stability_matrix
#from CPAC.basc.basc import group_stability_matrix, individual_group_clustered_maps, individual_stability_matrix, nifti_individual_stability, ndarray_to_vol, create_basc
#from CPAC.utils.utils import safe_shape

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def timeseries_bootstrap(tseries, block_size):
    """
    Generates a bootstrap sample derived from the input time-series.  Utilizes Circular-block-bootstrap method described in [1]_.

    Parameters
    ----------
    tseries : array_like
        A matrix of shapes (`M`, `N`) with `M` timepoints and `N` variables
    block_size : integer
        Size of the bootstrapped blocks

    Returns
    -------
    bseries : array_like
        Bootstrap sample of the input timeseries


    References
    ----------
    .. [1] P. Bellec; G. Marrelec; H. Benali, A bootstrap test to investigate
       changes in brain connectivity for functional MRI. Statistica Sinica,
       special issue on Statistical Challenges and Advances in Brain Science,
       2008, 18: 1253-1268.

    Examples
    --------

    >>> x = np.arange(50).reshape((5,10)).T
    >>> sample_bootstrap(x,3)
    array([[ 7, 17, 27, 37, 47],
           [ 8, 18, 28, 38, 48],
           [ 9, 19, 29, 39, 49],
           [ 4, 14, 24, 34, 44],
           [ 5, 15, 25, 35, 45],
           [ 6, 16, 26, 36, 46],
           [ 0, 10, 20, 30, 40],
           [ 1, 11, 21, 31, 41],
           [ 2, 12, 22, 32, 42],
           [ 4, 14, 24, 34, 44]])

    """
    import numpy as np

    k = np.ceil(float(tseries.shape[0])/block_size)
    r_ind = np.floor(np.random.rand(1,k)*tseries.shape[0])

    blocks = np.dot(np.arange(0,block_size)[:,np.newaxis], np.ones([1,k]))
    block_offsets = np.dot(np.ones([block_size,1]), r_ind)

    block_mask = (blocks + block_offsets).flatten('F')[:tseries.shape[0]]

    block_mask = np.mod(block_mask, tseries.shape[0])

    return tseries[block_mask.astype('int'), :]

def standard_bootstrap(dataset):
    """
    Generates a bootstrap sample from the input dataset

    Parameters
    ----------
    dataset : array_like
        A matrix of where dimension-0 represents samples

    Returns
    -------
    bdataset : array_like
        A bootstrap sample of the input dataset

    Examples
    --------
    """
    n = dataset.shape[0]
    b = np.random.random_integers(0, high=n-1, size=n)
    return dataset[b]

def cluster_timeseries(X, n_clusters, similarity_metric = 'k_neighbors', affinity_threshold = 0.0, neighbors = 10):
    """
    Cluster a given timeseries

    Parameters
    ----------
    X : array_like
        A matrix of shape (`N`, `M`) with `N` samples and `M` dimensions
    n_clusters : integer
        Number of clusters
    similarity_metric : {'k_neighbors', 'correlation', 'data'}
        Type of similarity measure for spectral clustering.  The pairwise similarity measure
        specifies the edges of the similarity graph. 'data' option assumes X as the similarity
        matrix and hence must be symmetric.  Default is kneighbors_graph [1]_ (forced to be
        symmetric)
    affinity_threshold : float
        Threshold of similarity metric when 'correlation' similarity metric is used.

    Returns
    -------
    y_pred : array_like
        Predicted cluster labels

    Examples
    --------


    References
    ----------
    .. [1] http://scikit-learn.org/dev/modules/generated/sklearn.neighbors.kneighbors_graph.html


    if similarity_metric == 'correlation':
        # Calculate empirical correlation matrix between samples
        Xn = X - X.mean(1)[:,np.newaxis]
        Xn = Xn/np.sqrt( (Xn**2.).sum(1)[:,np.newaxis] )
        C_X = np.dot(Xn, Xn.T)
        C_X[C_X < affinity_threshold] = 0
        from scipy.sparse import lil_matrix
        C_X = lil_matrix(C_X)
    elif similarity_metric == 'data':
        C_X = X
    elif similarity_metric == 'k_neighbors':
        from sklearn.neighbors import kneighbors_graph
        C_X = kneighbors_graph(X, n_neighbors=neighbors)
        C_X = 0.5 * (C_X + C_X.T)
    else:
        raise ValueError("Unknown value for similarity_metric: '%s'." % similarity_metric)

    #sklearn code is not stable for bad clusters which using correlation as a stability metric
    #tends to give for more info see:
    #http://scikit-learn.org/dev/modules/clustering.html#spectral-clustering warning
    #from sklearn import cluster
    #algorithm = cluster.SpectralClustering(k=n_clusters, mode='arpack')
    #algorithm.fit(C_X)
    #y_pred = algorithm.labels_.astype(np.int)

    from python_ncut_lib import ncut, discretisation
    eigen_val, eigen_vec = ncut(C_X, n_clusters)
    eigen_discrete = discretisation(eigen_vec)

    #np.arange(n_clusters)+1 isn't really necessary since the first cluster can be determined
    #by the fact that the each cluster is a disjoint set
    y_pred = np.dot(eigen_discrete.toarray(), np.diag(np.arange(n_clusters))).sum(1)

    """
#    # sampledata=generate_blobs()
#    #X= bg_func
#
#    # normalize dataset for easier parameter selection
#    X = pd.DataFrame(X)
#    X = StandardScaler().fit_transform(X)
#    spectral = cluster.SpectralClustering(n_clusters=n_clusters, eigen_solver='arpack', random_state = 5, affinity="nearest_neighbors", n_neighbors = 10, assign_labels='discretize')
#
#
#
#    #t0 = time.time()
#    spectral.fit(X)
#    #t1 = time.time()
#    if hasattr(spectral, 'labels_'):
#        y_pred = spectral.labels_.astype(np.int)
#    else:
#        y_pred = spectral.predict(X)



##########################
    #Set affinity type, and affinity threshold
    X = pd.DataFrame(X)
    X_dist = sp.spatial.distance.pdist(X, metric = similarity_metric)
    X_dist = sp.spatial.distance.squareform(X_dist)
    sim_matrix=1-X_dist
    sim_matrix[sim_matrix<affinity_threshold]=0

    spectral = cluster.SpectralClustering(n_clusters, eigen_solver='arpack', random_state = 5, affinity="precomputed", assign_labels='discretize')

    #change distribution to be normalized between 0 and 1 before passing into the spectral clustering algorithm


    #t0 = time.time()
    spectral.fit(sim_matrix)
    #t1 = time.time()

    y_pred = spectral.labels_.astype(np.int)




    return y_pred


def cross_cluster_timeseries(data1, data2, n_clusters, similarity_metric):


    """
    Cluster a timeseries dataset based on its relationship to a second timeseries dataset

    Parameters
    ----------
    data1 : array_like
        A matrix of shape (`N`, `M`) with `N1` samples and `M1` dimensions.
        This is the matrix to receive cluster assignment
    data2 : array_like
        A matrix of shape (`N`, `M`) with `N2` samples and `M2` dimensions.
        This is the matrix with which distances will be calculated to assign clusters to data1
    n_clusters : integer
        Number of clusters
    similarity_metric : {'euclidean', 'correlation', 'minkowski', 'cityblock', 'seuclidean'}
        Type of similarity measure for distance matrix.  The pairwise similarity measure
        specifies the edges of the similarity graph. 'data' option assumes X as the similarity
        matrix and hence must be symmetric.  Default is kneighbors_graph [1]_ (forced to be
        symmetric)
    affinity_threshold : float
        Threshold of similarity metric when 'correlation' similarity metric is used.

    Returns
    -------
    y_pred : array_like
        Predicted cluster labels


    Examples
    --------
    np.random.seed(30)
    offset = np.random.randn(30)
    x1 = np.random.randn(200,30) + 2*offset
    x2 = np.random.randn(100,30) + 44*np.random.randn(30)
    x3 = np.random.randn(400,30)
    sampledata1 = np.vstack((x1,x2,x3))

    np.random.seed(99)
    offset = np.random.randn(30)
    x1 = np.random.randn(200,30) + 2*offset
    x2 = np.random.randn(100,30) + 44*np.random.randn(30)
    x3 = np.random.randn(400,30)
    sampledata2 = np.vstack((x1,x2,x3))

    cross_cluster(sampledata1, sampledata2, 3, 'euclidean')


    References
    ----------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html#scipy.spatial.distance.cdist
    http://scikit-learn.org/stable/modules/clustering.html#spectral-clustering
    """


#    #Simulating data
#    blobs1 = generate_simple_blobs(27)
#    blobs2 = generate_simple_blobs(30)
#
#    blobs2 = blobs2[0:120,:]
#    blobs2[0:50,:] = blobs1[0:50,]
#
#    blobs1=blobs1[0:50,:]
#    blobs2=blobs2[0:40,:]
#
#    data1=blobs1
#    data2=blobs2

    data1_df = pd.DataFrame(data1)
    data2_df = pd.DataFrame(data2)

#    plt.imshow(sim_matrix)
#
#    a=np.hstack(y_pred)
#    plt.hist(a,bins='auto')

    #FIGURE OUT HOW TO CREATE SYMMETRIC AND 0-1 NORMED SIMILARITY MATRIX
    #CUTOFF AFFINITY AT SPECIFIC LEVEL??
    #CUTOFF THE AFFINITY MATRIX RIGHT AFTER THE FIRST CLUSTERING

    dist_btwn_df_1_2 = sp.spatial.distance.cdist(data1_df, data2_df, metric = similarity_metric)
    dist_of_1 = sp.spatial.distance.pdist(dist_btwn_df_1_2, metric = similarity_metric)
    dist_matrix = sp.spatial.distance.squareform(dist_of_1)
    sim_matrix=1-dist_matrix
    sim_matrix[sim_matrix<0.3]=0

#    sim_of_1 = 1 - dist_of_1
#
#
#
#    sim_btwn_df_1_2 = 1 - dist_btwn_df_1_2
#    #set similarities below 0.3 to 0
#    sim_btwn_df_1_2[sim_btwn_df_1_2<0.3] = 0
#
#    #corr = np.corrcoef(dist)
#    dist_of_1 = sp.spatial.distance.pdist(dist_btwn_df_1_2, metric = similarity_metric)
#    sim_of_1 = 1 - dist_of_1
#    dist_matrix = sp.spatial.distance.squareform(dist_of_1)
#    sim_matrix=1-dist_matrix
#    dist_matrix_scaler= StandardScaler().fit_transform(dist_matrix)
#    dist_matrix_norm = normalize(dist_matrix)
#    delta = sqrt(data1_df.shape[0])
#    #y_pred = cluster_timeseries(dist, n_clusters, affinity = 'precomputed')
#    sq_sim_of1 = 1-dist_matrix_scaler
#
#
#    sq_sim_of1 = np.exp(-beta * sq_dist_of1 / sq_dist_of1.std())
#    sq_sim_of1= np.exp(- dist_matrix ** 2 / (2. * delta ** 2))
#
#
#    X = StandardScaler().fit_transform(X)

    spectral = cluster.SpectralClustering(n_clusters, eigen_solver='arpack', random_state = 5, affinity="precomputed", assign_labels='discretize')

    #change distribution to be normalized between 0 and 1 before passing into the spectral clustering algorithm


    #t0 = time.time()
    spectral.fit(sim_matrix)
    #t1 = time.time()

    y_pred = spectral.labels_.astype(np.int)

    return y_pred



def adjacency_matrix(cluster_pred):
    """
    Calculate adjacency matrix for given cluster predictions

    Parameters
    ----------
    cluster_pred : array_like
        A matrix of shape (`N`, `1`) with `N` samples

    Returns
    -------
    A : array_like
        Adjacency matrix of shape (`N`,`N`)

    Examples
    --------
    >>> import numpy as np
    >>> from CPAC.basc import cluster_adjacency_matrix
    >>> x = np.asarray([1, 2, 2, 3, 1])[:,np.newaxis]
    >>> cluster_adjacency_matrix(x).astype('int')
    array([[1, 0, 0, 0, 1],
           [0, 1, 1, 0, 0],
           [0, 1, 1, 0, 0],
           [0, 0, 0, 1, 0],
           [1, 0, 0, 0, 1]])

    """
    x = cluster_pred.copy()

    if(len(x.shape) == 1):
        x = x[:, np.newaxis]

    # Force the cluster indexing to be positive integers
    if(x.min() <= 0):
        x += -x.min() + 1

    A = np.dot(x**-1., x.T) == 1

    return A

def cluster_matrix_average(M, cluster_assignments):
    """
    Calculate the average element value within a similarity matrix for each cluster assignment, a measure
    of within cluster similarity.  Self similarity (diagonal of similarity matrix) is removed.

    Parameters
    ----------
    M : array_like
    cluster_assignments : array_like

    Returns
    -------
    s : array_like

    Examples
    --------
    >>> import numpy as np
    >>> from CPAC import basc
    >>> S = np.arange(25).reshape(5,5)
    >>> assign = np.array([0,0,0,1,1])
    >>> basc.cluster_matrix_average(S, assign)
    array([  6.,   6.,   6.,  21.,  21.])

    """

    if np.any(np.isnan(M)):
        np.save('bad_M.npz', M)
        raise ValueError('M matrix has a nan value')

    cluster_ids = np.unique(cluster_assignments)
    s = np.zeros((cluster_ids.shape[0], cluster_assignments.shape[0]), dtype='float64')
    s_idx = 0
    for cluster_id in cluster_ids:
        s[s_idx, :] = M[:,cluster_assignments == cluster_id].mean(1)
        s_idx += 1
#        k = (cluster_assignments == cluster_id)[:, np.newaxis]
#        print 'Cluster %i size: %i' % (cluster_id, k.sum())
#        K = np.dot(k,k.T)
#        K[np.diag_indices_from(K)] = False
#        if K.sum() == 0: # Voxel with its own cluster
#            s[k[:,0]] = 0.0
#        else:
#            s[k[:,0]] = M[K].mean()


    return s

def individual_stability_matrix(Y1, n_bootstraps, k_clusters, Y2=None, cross_cluster=False, cbb_block_size = None, affinity_threshold = 0.5):
    """
    Calculate the individual stability matrix of a single subject by bootstrapping their time-series

    Parameters
    ----------
    Y1 : array_like
        A matrix of shape (`V`, `N`) with `V` voxels `N` timepoints
    Y2 : array_like
        A matrix of shape (`V`, `N`) with `V` voxels `N` timepoints
        For Cross-cluster solutions- this will be the matrix by which Y1 is clustered
    n_bootstraps : integer
        Number of bootstrap samples
    k_clusters : integer
        Number of clusters
    cbb_block_size : integer, optional
        Block size to use for the Circular Block Bootstrap algorithm
    affinity_threshold : float, optional
        Minimum threshold for similarity matrix based on correlation to create an edge

    Returns
    -------
    S : array_like
        A matrix of shape (`V1`, `V1`), each element v1_{ij} representing the stability of the adjacency of voxel i with voxel j
    """
    if affinity_threshold < 0.0:
        raise ValueError('affinity_threshold %d must be non-negative value' % affinity_threshold)

    #flipped the N and V values bc originally data was being put in transposed
    N1 = Y1.shape[1]
    V1 = Y1.shape[0]

    if(cbb_block_size is None):
        cbb_block_size = int(np.sqrt(N1))

    S = np.zeros((V1, V1))

    if (cross_cluster is True):
        for bootstrap_i in range(n_bootstraps):
            N2 = Y2.shape[1]
            cbb_block_size2 = int(np.sqrt(N2))
            Y_b1 = timeseries_bootstrap(Y1, cbb_block_size)
            Y_b2 = timeseries_bootstrap(Y2, cbb_block_size2)
            S += adjacency_matrix(cross_cluster_timeseries(Y_b1, Y_b2, k_clusters, similarity_metric = 'correlation'))
        S /= n_bootstraps
    else:
        for bootstrap_i in range(n_bootstraps):
            Y_b1 = timeseries_bootstrap(Y1, cbb_block_size)
            S += adjacency_matrix(cluster_timeseries(Y_b1, k_clusters, similarity_metric = 'correlation', affinity_threshold = affinity_threshold)[:,np.newaxis])
        S /= n_bootstraps

    return S

def data_compression(Y1, mask, output_size):
    """
    data : array_like
         A matrix of shape (`V`, `N`) with `V` voxels `N` timepoints
         The functional dataset that needs to be reduced
    output_size : integer
        The number of elements that the data should be reduced to
        
    """
    
    #%%
    #Pseudocode-
    # Option 1- Take in the functional data and mask, and mask here. and convert data to numpy array
    # Option 2- Take in the masked functional data for ROI1, and ROI2, and separately reduce dimensionality.
    #%%
    #Inputs- Data, mask1, mask2, output_size

###################################################################
## Transform nifti files to a data matrix with the NiftiMasker
#from nilearn import input_data
#
## The NiftiMasker will extract the data on a mask. We do not have a
## mask, hence we need to compute one.
##
## This is resting-state data: the background has not been removed yet,
## thus we need to use mask_strategy='epi' to compute the mask from the
## EPI images
#nifti_masker = input_data.NiftiMasker(memory='nilearn_cache',
#                                      mask_strategy='epi', memory_level=1,
#                                      standardize=False)
#
#func_filename = dataset.func[0]
## The fit_transform call computes the mask and extracts the time-series
## from the files:
#fmri_masked = nifti_masker.fit_transform(func_filename)
#
## We can retrieve the numpy array of the mask
#mask = nifti_masker.mask_img_.get_data().astype(bool)


##################################################################
# Perform Ward clustering
# -----------------------
#
# We use spatially-constrained Ward clustering. For this, we need to
# compute from the mask a matrix giving the voxel-to-voxel connectivity

# Compute connectivity matrix: which voxel is connected to which
    from sklearn.feature_extraction import image
    shape = mask.shape
    connectivity = image.grid_to_graph(n_x=shape[0], n_y=shape[1],
                                   n_z=shape[2], mask=mask)


##################################################################
# Then we use FeatureAgglomeration from scikit-learn. Indeed, the voxels
# are the features of the data matrix.
#
# In addition, we use caching. As a result, the clustering doesn't have
# to be recomputed later.

# Computing the ward for the first time, this is long...
    from sklearn.cluster import FeatureAgglomeration
# If you have scikit-learn older than 0.14, you need to import
# WardAgglomeration instead of FeatureAgglomeration
    import time
    start = time.time()
    ward = FeatureAgglomeration(n_clusters=output_size, connectivity=connectivity,
                            linkage='ward')
    ward.fit(Y1.T)
    print("Ward agglomeration compressing voxels into clusters: %.2fs" % (time.time() - start))

## Compute the ward with more clusters, should be faster as we are using
## the caching mechanism
#start = time.time()
#ward = FeatureAgglomeration(n_clusters=2000, connectivity=connectivity,
#                            linkage='ward', memory='nilearn_cache')
#ward.fit(fmri_masked)
#print("Ward agglomeration 2000 clusters: %.2fs" % (time.time() - start))

##################################################################
# Visualize results
# ------------------
#
# First we display the labels of the clustering in the brain.
#
# To visualize results, we need to transform the clustering's labels back
# to a neuroimaging volume. For this, we use the NiftiMasker's
# inverse_transform method.
#from nilearn.plotting import plot_roi, plot_epi, show

## Unmask the labels
#
## Avoid 0 label
#labels = ward.labels_ + 1
#labels_img = nifti_masker.inverse_transform(labels)
#
##from nilearn.image import mean_img
##mean_func_img = mean_img(func_filename)
##
##
##first_plot = plot_roi(labels_img, mean_func_img, title="Ward parcellation",
##                      display_mode='xz')
##
### common cut coordinates for all plots
##cut_coords = first_plot.cut_coords
#
###################################################################
## labels_img is a Nifti1Image object, it can be saved to file with the
## following code:
#labels_img.to_filename('parcellation.nii')


##################################################################
# Second, we illustrate the effect that the clustering has on the
# signal. We show the original data, and the approximation provided by
# the clustering by averaging the signal on each parcel.
#
# As you can see below, this approximation is very good, although there
# are only 2000 parcels, instead of the original 60000 voxels

## Display the original data
#plot_epi(nifti_masker.inverse_transform(fmri_masked[0]),
#         cut_coords=cut_coords,
#         title='Original (%i voxels)' % fmri_masked.shape[1],
#         vmax=fmri_masked.max(), vmin=fmri_masked.min(),
#         display_mode='xz')

# A reduced data can be create by taking the parcel-level average:
# Note that, as many objects in the scikit-learn, the ward object exposes
# a transform method that modifies input features. Here it reduces their
# dimension
    data_reduced = ward.transform(Y1.T)
    return data_reduced
#
## Display the corresponding data compressed using the parcellation
#fmri_compressed = ward.inverse_transform(data_reduced)
#compressed_img = nifti_masker.inverse_transform(fmri_compressed[0])
#
#plot_epi(compressed_img, cut_coords=cut_coords,
#         title='Compressed representation (2000 parcels)',
#         vmax=fmri_masked.max(), vmin=fmri_masked.min(),
#         display_mode='xz')
#
#show()
