import numpy as np

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
    
    """

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
#    from sklearn import cluster
#    algorithm = cluster.SpectralClustering(k=n_clusters, mode='arpack')
#    algorithm.fit(C_X)
#    y_pred = algorithm.labels_.astype(np.int)

    from python_ncut_lib import ncut, discretisation
    eigen_val, eigen_vec = ncut(C_X, n_clusters)
    eigen_discrete = discretisation(eigen_vec)

    #np.arange(n_clusters)+1 isn't really necessary since the first cluster can be determined
    #by the fact that the each cluster is a disjoint set
    y_pred = np.dot(eigen_discrete.toarray(), np.diag(np.arange(n_clusters))).sum(1)
    
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
    s = np.zeros( (cluster_ids.shape[0], cluster_assignments.shape[0]), dtype='float64')
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

def individual_stability_matrix(Y, n_bootstraps, k_clusters, cbb_block_size = None, affinity_threshold = 0.5):
    """
    Calculate the individual stability matrix of a single subject by bootstrapping their time-series
    
    Parameters
    ----------
    Y : array_like
        A matrix of shape (`N`, `V`) with `N` timepoints and `V` voxels
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
        A matrix of shape (`V`, `V`), each element v_{ij} representing the stability of the adjacency of voxel i with voxel j
    """
    if affinity_threshold < 0.0:
        raise ValueError('affinity_threshold %d must be non-negative value' % affinity_threshold)
    
    N = Y.shape[0]
    V = Y.shape[1]
    
    if(cbb_block_size is None):
        cbb_block_size = int(np.sqrt(N))

    S = np.zeros((V,V))
    for bootstrap_i in range(n_bootstraps):
        Y_b = timeseries_bootstrap(Y, cbb_block_size)
        S += adjacency_matrix(cluster_timeseries(Y_b.T, k_clusters, similarity_metric = 'correlation', affinity_threshold = affinity_threshold)[:,np.newaxis])
    S /= n_bootstraps

    return S