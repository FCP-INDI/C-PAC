# CPAC/network_centrality/utils.py
#
# Network centrality utilities

# Method to return recommended block size based on memory restrictions 
def calc_blocksize(timeseries, memory_allocated=None, 
                   include_full_matrix=False, sparsity_thresh=0.0):
    '''
    Method to calculate blocksize to calculate correlation matrix
    as per the memory allocated by the user. By default, the block
    size is 1000 when no memory limit is specified.

    If memory allocated is specified, then block size is calculated
    as memory allocated subtracted by the memory of the timeseries 
    and centrality output, then divided by the size of one correlation 
    map. That is how many correlation maps can we calculate simultaneously 
    in memory?

    Parameters
    ----------
    timeseries : numpy array
       timeseries data: `nvoxs` x `ntpts`
    memory_allocated : float
       memory allocated in GB for degree centrality
    include_full_matrix : boolean
        Boolean indicating if we're using the entire correlation matrix
        in RAM (needed during eigenvector centrality).
        Default is False
    sparsity_thresh : float
        a number between 0 and 1 that represents the number of
        connections to keep during sparsity thresholding.
        Default is 0.0.

    Returns
    -------
    block_size : an integer
      size of block for matrix calculation
    '''

    # Import packages
    import numpy as np

    # Init variables
    block_size = 1000   # default

    nvoxs   = timeseries.shape[0]
    ntpts   = timeseries.shape[1]
    nbytes  = timeseries.dtype.itemsize

    # If we need the full matrix for centrality calculation
    if include_full_matrix:
        memory_for_full_matrix = nvoxs * nvoxs * nbytes
    # Otherwise, we're doing it in blocks
    else:
        memory_for_full_matrix = 0

    # Memory variables
    memory_for_timeseries   = nvoxs * ntpts * nbytes
    memory_for_output       = 2 * nvoxs * nbytes    # bin and wght outputs
    needed_memory = memory_for_timeseries + \
                    memory_for_output + \
                    memory_for_full_matrix

    if memory_allocated:
        available_memory = memory_allocated * 1024.0**3  # assume it is in GB
        ## memory_for_block = # of seed voxels * nvoxs * nbytes
        block_size = int( (available_memory - needed_memory)/(nvoxs*nbytes) )
        # If we're doing degree/sparisty thresholding, calculate block_size
        if sparsity_thresh:
            # k - block_size, v - nvoxs, d - nbytes, m - memory_allocated
            # Solve for k: (-d/2 - 20.5)*k^2 + (41*v + d*v -d/2 - 20.5)*k - m = 0
            coeffs = np.zeros(3)
            coeffs[0] = -nbytes/2 - 20.5
            coeffs[1] = 41*nvoxs + nbytes*nvoxs - nbytes/2 - 20.5
            coeffs[2] = -(available_memory - needed_memory)
            roots = np.roots(coeffs)
            # If roots are complex, then the block_size needed to occupy all of
            # the available memory is bigger than the number of voxels.
            # So set block_size = nvoxs
            if np.iscomplex(roots[0]):
                block_size = nvoxs
            # If the roots are real, test the roots for condition
            else:
                root = roots[np.where(roots <= nvoxs)]
                root = root[np.where(root > 0)]
                if len(root) == 1:
                    block_size = np.floor(root[0])
                else:
                    block_size = 1000

    # Test if calculated block size is beyond max/min limits
    if block_size > nvoxs:
        block_size = nvoxs
    elif block_size < 1:
        memory_usage = (needed_memory + 2.0*nvoxs*nbytes)/1024.0**3
        raise MemoryError('Not enough memory available to perform degree '\
                          'centrality. Need a minimum of %.2fGB' % memory_usage)

    # Convert block_size to an integer before returning
    block_size = int(block_size)

    # Return memory usage and block size
    if sparsity_thresh:
        # Calculate RAM usage by blocking algorithm
        m = (-nbytes/2 - 20.5)*block_size**2 + \
            (41*nvoxs + nbytes*nvoxs - nbytes/2 - 20.5)*block_size
        # Calculate RAM usage by sparse matrix at end
        max_conns = nvoxs**2-nvoxs
        # Max number of connections * i + j (32-bit ints) + w
        m2 = np.round(max_conns*sparsity_thresh)*(4 + 4 + nbytes)
        if m2 > m:
            m = m2
        memory_usage = (needed_memory + m)/1024.0**3
    else:
        memory_usage = (needed_memory + block_size*nvoxs*nbytes)/1024.0**3

    # Print information
    print 'block_size -> %i voxels' % block_size
    print '# of blocks -> %i' % np.ceil(float(nvoxs)/block_size)
    print 'expected usage -> %.2fGB' % memory_usage

    return block_size


# Method to calculate correlation coefficient from (one or two) datasets
def calc_corrcoef(X, Y=None):
    '''
    Method to calculate correlation 
    Each of the columns in X will be correlated 
    with each of the columns in Y. Each column 
    represents a variable, with the rows containing the observations.

    Parameters
    ----------
    X : numpy array
       array of shape x1, x2
    Y : numpy array
      array of shape y1, y2

    Returns
    -------
    r : numpy array
      array containing correlation values of shape x2, y2
    '''

    # Import packages
    import numpy as np

    if Y is None:
        Y = X

    if X.shape[0] != Y.shape[0]:
        raise Exception("X and Y must have the same number of rows.")

    X = X.astype(float)
    Y = Y.astype(float)

    X -= X.mean(axis=0)[np.newaxis,...]
    Y -= Y.mean(axis=0)

    xx = np.sum(X**2, axis=0)
    yy = np.sum(Y**2, axis=0)

    r = np.dot(X.T, Y)/np.sqrt(np.multiply.outer(xx,yy))

    return r


# Method to cluster the data (used in lFCD)
def cluster_data(img, thr, xyz_a, k=26):
    '''docstring for cluster_data'''

    # Import packages
    from scipy.sparse import coo_matrix, cs_graph_components
    import numpy as np

    # Threshold the entire correlation map and find connected components, store this in sparse matrix
    val_idx = img > thr                 # store valid indices
    xyz_th = xyz_a[val_idx]             # find the 3D indices corresponding to the above threshold voxels
    i,j,d = graph_3d_grid(xyz_th, k=k)  # find the connected components for the above threshold voxels
    nvoxs = xyz_th.shape[0]             # store the number of correlated voxels in entire network
    adj = coo_matrix((d, (i,j)), shape=(nvoxs,nvoxs)) # and store the connected nodes and weights in sparse matrix

    # Identify the connected components (clusters) within the graph
    nc, labels = cs_graph_components(adj)

    # Copy the node labels to their voxel equivalents
    lbl_img = np.zeros(img.shape)           # init lbl_img - map to store label data
    # add 2 so that labels corresponding to unconnected voxels (-2)
    # will be zero in lbl_img, and label==0 will now equal 2
    lbl_img[val_idx] = labels + 2 
    return lbl_img


# Convert probability threshold value to correlation threshold
def convert_pvalue_to_r(scans, threshold):
    '''
    Method to calculate correlation threshold from p_value

    Parameters
    ----------
    scans : int
        Total number of scans in the data
    threshold : float
        input p_value

    Returns
    -------
    rvalue : float
        correlation threshold value 
    '''

    # Import packages
    import scipy.stats as s
    import math

    print "p_value ->", threshold
    x = 1-threshold/2
    dof = scans-2
    #Inverse Survival Function (Inverse of SF)
    tvalue = s.t.isf(x, dof)
    rvalue = math.sqrt(math.pow(tvalue, 2)/(dof+ math.pow(tvalue,2)))

    return rvalue


# Borrowed from nipy.graph.graph
# https://github.com/nipy/nipy/blob/master/nipy/algorithms/graph/graph.py
def graph_3d_grid(xyz, k=18):
    '''
    Utility that computes the six neighbors on a 3d grid

    Parameters
    ----------
    xyz: array of shape (n_samples, 3); grid coordinates of the points
    k: neighboring system, equal to 6, 18, or 26

    Returns
    -------
    i, j, d 3 arrays of shape (E),
            where E is the number of edges in the resulting graph
            (i, j) represent the edges, d their weights
    '''

    # Import packages
    import numpy as np

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


# Function to map a centrality matrix to a nifti image
def map_centrality_matrix(centrality_matrix, aff, mask, template_type):
    '''
    Method to map centrality matrix to a nifti image

    Parameters
    ----------
    centrality_matrix : tuple (string, array_like)
        tuple containing matrix name and degree/eigenvector centrality matrix
    aff : ndarray
        Affine matrix of the input data
    mask : ndarray
        Mask or roi data matrix
    template_type : int
        type of template: 0 for mask, 1 for roi

    Returns
    -------
    out_file : string (nifti image)
        nifti image mapped from the centrality matrix

    Raises
    ------
    Exception
    '''

    import nibabel as nib
    import os
    import numpy as np

    try:
        out_file, matrix = centrality_matrix

        out_file = os.path.join(os.getcwd(), out_file + '.nii.gz')
        sparse_m = np.zeros((mask.shape), dtype=float)

        print 'mapping centrality matrix to nifti image...', out_file

        if int(template_type) == 0:
            cords = np.argwhere(mask)
            index=0
            for val in cords:
                x,y,z=val
                sparse_m[x,y,z]= matrix[index]
                index+=1

        elif int(template_type) == 1:
            nodes = np.unique(mask).tolist()
            nodes.sort()
            index = 0
            for n in nodes:
                if n> 0:
                    cords = np.argwhere(mask==n)
                    for val in cords:
                        x,y,z = val
                        if isinstance(matrix[index], list):
                            sparse_m[x,y,z] = matrix[index][0]
                        else:
                            sparse_m[x,y,z]=matrix[index]
                    index+=1

        nifti_img = nib.Nifti1Image(sparse_m, aff)
        nifti_img.to_filename(out_file)

        return out_file
    except:
        print 'Error in mapping centrality matrix to nifti image'
        raise


# Function to actually do the list merging
def merge_lists(deg_list=[],eig_list=[],lfcd_list=[]):
    merged_list = []
    merged_list.extend(deg_list)
    merged_list.extend(eig_list)
    merged_list.extend(lfcd_list)

    return merged_list

