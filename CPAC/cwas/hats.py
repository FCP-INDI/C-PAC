import numpy as np

# TODO: make better exception class?
def check_rank(x):
    k    = x.shape[1]
    rank = np.linalg.matrix_rank(x)
    if rank < k:
        raise Exception("matrix is rank deficient (rank %i vs cols %i)" % (rank, k))

def add_intercept(x):
    """
    Adds an intercept column to the left of the matrix
    
    Paramaters
    ----------
    x : ndarray
        Design matrix (e.g. with 1st column as your intercept)
    
    Returns
    -------
    x : ndarray
    """
    uno = np.ones((x.shape[0],1))   # intercept
    xx  = np.hstack((uno,x))        # design matrix
    return xx

def hatify(x):
    """
    Distance-based hat matrix
    
    Paramaters
    ----------
    x : ndarray
        Design matrix (e.g. with 1st column as your intercept)
    
    Notes
    -----
    This function assumes that the input is not rank-deficient.
    
    Returns
    -------
    H : ndarray
        This will be a `x.shape[0]` by `x.shape[0]` matrix.
    """
    Q1, R1 = np.linalg.qr(x)
    H = Q1.dot(Q1.T)
    return H

def permute_design(x, cols, indexperm):
    """docstring for permute_design"""
    Xj          = x.copy()
    Xp          = np.take(Xj[:,cols], indexperm, axis=0)    
    Xj[:,cols]  = Xp    
    return Xj

# make sure this function doesn't overwrite
# the original x
def gen_h(x, cols=None, indexperm=None):
    """
    Permuted hat matrix
    
    Parameters
    ----------
    x : ndarray
        Design matrix (e.g. with 1st column as your intercept)
    cols : list (optional)
        Columns to be permuted (if `indexperm` is specified)
    indexperm : list (optional)
        Re-ordering (permuting) of rows in `x`
    
    Returns
    -------
    H : ndarray
        This will be a `x.shape[0]` by `x.shape[0]` matrix.
    """
    if indexperm is not None:
        x = permute_design(x, cols, indexperm)
    H = hatify(x)
    return H
    
def gen_h2(x, cols, indexperm=None):
    """
    Permuted regressor-specific hat matrix
    
    Parameters
    ----------
    x : ndarray
        Design matrix (e.g. with 1st column as your intercept)
    cols : list
        Columns to be permuted (if `indexperm` is specified)
    indexperm : list (optional)
        Re-ordering (permuting) of rows in `x`
    
    Returns
    -------
    H2 : ndarray
        This will be a `x.shape[0]` by `x.shape[0]` matrix.
    """
    # H
    H = gen_h(x, cols, indexperm)
    # H2
    # take H and subtract it by everything other than the columns of interest
    other_cols = [ i for i in range(x.shape[1]) if i not in cols ]
    Xj = x[:,other_cols]
    H2 = H - hatify(Xj)
    return H2

