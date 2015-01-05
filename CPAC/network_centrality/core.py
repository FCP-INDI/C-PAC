# coding: utf-8

"""
These are function nodes of sorts that compute a centrality measure.
"""

####
# Degree Centrality
####

import numpy as np
#import pyximport
#pyximport.install(setup_args={'include_dirs': [np.get_include()]})
from CPAC.network_centrality.thresh_and_sum import *


def degree_centrality(corr_matrix, r_value, method, out=None):
    """
    Calculate centrality for the rows in the corr_matrix using
    a specified correlation threshold. The centrality output can 
    be binarized or weighted.
    
    Paramaters
    ---------
    corr_matrix : numpy.ndarray
    r_value : float
    method : str
        Can be 'binarize' or 'weighted'
    out : numpy.ndarray (optional)
        If specified then should have shape of `corr_matrix.shape[0]`
    
    Returns
    -------
    out : numpy.ndarray
    """
    
    if method not in ["binarize", "weighted"]:
        raise Exception("Method must be one of binarize or weighted and not %s" % method)
    
    if corr_matrix.dtype.itemsize == 8:
        dtype   = "double"
        r_value = np.float64(r_value)
    else:
        dtype   = "float"
        r_value = np.float32(r_value)
    
    if out is None:
        out = np.zeros(corr_matrix.shape[0], dtype=corr_matrix.dtype)
    print 'about to call thresh_and_sum'
    func_name   = "centrality_%s_%s" % (method, dtype)
    func        = globals()[func_name]
    func(corr_matrix, out, r_value)
    
    return out


def fast_degree_centrality(m):
    from numpy import linalg as LA
    
    ntpts = m.shape[0]
    nvoxs = m.shape[1]
    
    wts   = np.ones((nvoxs,1))/np.sqrt(nvoxs)   # node weights
    part1 = m.dot(wts)                          # part one of M*v
    part2 = m.T.dot(part1)                      # part two of M*v
    
    return part2



####
# Eigenvector Centrality
####

def eigenvector_centrality(corr_matrix, 
                           r_value=None, 
                           method=None, 
                           to_transform=True, 
                           ret_eigenvalue=False):
    """
    Examples
    --------
    >>> # Simulate Data
    >>> import numpy as np
    >>> ntpts = 100; nvoxs = 1000
    >>> m = np.random.random((ntpts,nvoxs))
    >>> mm = m.T.dot(m) # note that need to generate connectivity matrix here
    >>> # Execute
    >>> from CPAC.network_centrality.core import eigenvector_centrality
    >>> eigenvector = slow_eigenvector_centrality(mm)
    """
    from scipy.sparse import linalg as LA
    from scipy.sparse import csc_matrix
    
    if method not in ["binarize", "weighted"]:
        raise Exception("Method must be one of binarize or weighted and not %s" % method)
    
    # Don't transform if binarize
    if method == "binarize" and to_transform is True:
        to_transform = False
    
    if corr_matrix.dtype.itemsize == 8:
        dtype   = "double"
        r_value = np.float64(r_value)
    else:
        dtype   = "float"
        r_value = np.float32(r_value)
    
    # Create function name by gathering it's parts
    # Threshold? Transform? Method? Datatype?
    func_args = []; func_args.append(corr_matrix)
    func_elems = []
    if r_value is not None:
        func_args.append(r_value)
        func_elems.append("thresh")

    if to_transform is True:
        func_elems.append("transform")

    func_elems.append(method)
    func_elems.append(dtype)
    # Combine to create function name
    func_name = "_".join(func_elems)
    
    # Execute function
    func        = globals()[func_name]
    func(*func_args)
    
    #using scipy method, which is a wrapper to the ARPACK functions
    #http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html
    eigenValue, eigenVector = LA.eigsh(corr_matrix, k=1, which='LM', maxiter=1000)
    
    if ret_eigenvalue:
        return eigenValue, np.abs(eigenVector)
    else:
        return np.abs(eigenVector)


def fast_eigenvector_centrality(m, maxiter=99, verbose=True):
    """
    The output here is based on a transfered correlation matrix of m.
    Where it equals (1+r)/2.
    
    References
    ----------
    .. [1] Wink, A.M., de Munck, J.C., van der Werf, Y.D., van den Heuvel, O.A., Barkhof, F., 2012. Fast Eigenvector Centrality Mapping of Voxel-Wise Connectivity in Functional Magnetic Resonance Imaging: Implementation, Validation, and Interpretation. Brain Connectivity 2, 265–274.
    
    Examples
    --------
    >>> # Simulate Data
    >>> import numpy as np
    >>> ntpts = 100; nvoxs = 1000
    >>> m = np.random.random((ntpts,nvoxs)) # note that don't need to generate connectivity matrix
    >>> # Execute
    >>> from CPAC.network_centrality.core import fast_eigenvector_centrality
    >>> eigenvector = fast_eigenvector_centrality(m)
    """
    from numpy import linalg as LA
    
    ntpts = m.shape[0]
    nvoxs = m.shape[1]
    
    # Initialize eigenvector estimate
    
    vprev = 0                                   # initialize previous ECM estimate
    vcurr = np.ones((nvoxs,1))/np.sqrt(nvoxs)   # initialize estimate with L2-norm == 1
    
    i = 0                                       # reset iteration counter
    dnorm = 1                                   # initial value for difference L2-norm
    cnorm = 0                                   # initial value to estimate L2-norm
    
    # Efficient power iteration
    
    while (i < maxiter) & (dnorm > cnorm):
        vprev = vcurr                           # start with previous estimate
        prevsum = np.sum(vprev)                 # sum of estimate
        vcurr_1 = m.dot(vprev)                  # part one of M*v
        vcurr_2 = m.T.dot(vcurr_1)              # part two of M*v
        vcurr_3 = vcurr_2 + prevsum             # adding sum -- same effect as [M+1]*v
        vcurr   = vcurr_3/LA.norm(vcurr_3,2)    # normalize L2-norm
                
        i += 1
        dnorm = LA.norm(vcurr-vprev, 2)
        cnorm = LA.norm(vcurr,2) * np.spacing(1)
        if verbose:
            print "iteration %02d, || v_i - v_(i-1) || / || v_i * epsilon || = %0.16f / %0.16f" % (i, dnorm, cnorm) #  some stats for the users
    
    if (i >= maxiter) & (dnorm > cnorm):
        print "Error: algorithm did not converge"   # TODO: convert to an exception?
    
    # now the vcurr value will be the ECM
    return vcurr
 
