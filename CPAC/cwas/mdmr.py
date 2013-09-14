import numpy as np
from hats import *

def permuted_index(n, strata=None):
    if strata is None:
        perms = np.random.permutation(n)
    else:
        perms = np.array(range(n))
        elems = np.unique(strata)
        for elem in elems:
            inds = perms[strata == elem]
            if len(inds) > 1:
                perms[strata == elem] = np.random.permutation(inds)
    return perms

def gen_perms(nperms, nobs, strata=None):
    perms = np.zeros((nperms, nobs), dtype=np.int)
    for i in range(nperms):
        perms[i,:] = permuted_index(nobs, strata)
    return perms

def add_original_index(perms):
    nobs = perms.shape[1]
    perms = np.vstack((range(nobs), perms))
    return perms

def gower_center(yDis):
    n = yDis.shape[0]
    I = np.eye(n,n)
    uno = np.ones((n,1))
    
    A = -0.5*(yDis**2)
    C = I - (1.0/n)*uno.dot(uno.T)
    G = C.dot(A).dot(C)
    
    return G

def gower_center_many(dmats):
    nobs    = np.sqrt(dmats.shape[0])
    ntests  = dmats.shape[1]
    Gs      = np.zeros_like(dmats)
    
    for i in range(ntests):
        Dmat    = dmats[:,i].reshape(nobs,nobs)
        Gs[:,i] = gower_center(Dmat).flatten()
    
    return Gs

def gen_h2_perms(x, cols, perms):
    nperms  = perms.shape[0]
    nobs    = perms.shape[1]
    
    H2perms = np.zeros((nobs**2, nperms))
    for i in range(nperms):
        H2 = gen_h2(x, cols, perms[i,:])
        H2perms[:,i] = H2.flatten()
    
    return H2perms

def gen_ih_perms(x, cols, perms):
    nperms  = perms.shape[0]
    nobs    = perms.shape[1]
    I       = np.eye(nobs,nobs)
    
    IHperms = np.zeros((nobs**2, nperms))
    for i in range(nperms):
        IH = I - gen_h(x, cols, perms[i,:])
        IHperms[:,i] = IH.flatten()
    
    return IHperms

def calc_ssq_fast(Hs, Gs, transpose=True):
    if transpose:
        ssq = Hs.T.dot(Gs)
    else:
        ssq = Hs.dot(Gs)
    return ssq

def ftest_fast(Hs, IHs, Gs, df_among, df_resid, **ssq_kwrds):
    SS_among = calc_ssq_fast(Hs, Gs, **ssq_kwrds)
    SS_resid = calc_ssq_fast(IHs, Gs, **ssq_kwrds)
    F = (SS_among/df_among)/(SS_resid/df_resid)
    return F

def calc_ssq_slow(H, G):
    ssq = np.trace(H.dot(G).dot(H))
    return ssq

def ftest_slow(H, IH, G, df_among, df_resid):
    SS_among = calc_ssq_slow(H, G)
    SS_resid = calc_ssq_slow(IH, G)
    F = (SS_among/df_among)/(SS_resid/df_resid)
    return F

def mdmr_single(yDis, *args, **kwrds):
    nobs = yDis.shape[0]
    ys   = yDis.reshape(nobs**2, 1)
    return mdmr(ys, *args, **kwrds)

def fperms_to_pvals(fstats, F_perms):
    nperms,ntests = F_perms.shape
    pvals = np.zeros(ntests)
    for i in range(ntests):
        j        = (F_perms[:,i] >= fstats[i]).sum().astype('float')
        pvals[i] = j/nperms
    return pvals

def mdmr(ys, x, cols, perms, strata=None, debug_output=False):
    """
    Multivariate Distance Matrix Regression
    
    Parameters
    ----------
    ys : ndarray
    x : ndarray
    perms : integer or ndarray
    strata : list or ndarray
    
    Returns
    --------
    ps : float
    Fs : float
    Fperms : ndarray
    perms : ndarray
    
    References
    -----------
    .. [1] Anderson, M. J. 2002. DISTML v.2: a FORTRAN computer program to calculate a distance-based multivariate analysis for a linear model. Dept. of Statistics University of Auckland. (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
    .. [2] Anderson, M. J. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology 26: 32-46.
    .. [3] Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed. Elsevier Science BV, Amsterdam.
    .. [4] McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to community data: a comment on distance-based redundancy analysis. Ecology 290-297.
    .. [5] Neter, J., M. H. Kutner, C. J. Nachtsheim, and W. Wasserman. 1996. Applied linear statistical models. 4th ed. Irwin, Chicago, Illinois.
    """
    check_rank(x)
    
    ntests  = ys.shape[1]
    nobs    = x.shape[0]
    if nobs != np.sqrt(ys.shape[0]):
        raise Exception("# of observations incompatible between x and ys")
    
    ## Distance matrix => Gower's centered matrix
    # G is similar to matrix of inner products from distances used in Partha Niyogi's
    # multidimensional scaling
    Gs = gower_center_many(ys)
    
    # Degrees of freedom
    df_among = len(cols)
    df_resid = nobs - x.shape[1]
    df_total = nobs - 1
    
    # Permutations
    if type(perms) is int:
        perms = gen_perms(perms, nobs, strata)
    perms  = add_original_index(perms)
    nperms = perms.shape[0]
    
    # Permuted versions of H2 and IH
    H2perms = gen_h2_perms(x, cols, perms)
    IHperms = gen_ih_perms(x, cols, perms)
    
    # Permutations of Fstats
    F_perms = ftest_fast(H2perms, IHperms, Gs,
                         df_among, df_resid)
    
    # F-values
    Fs = F_perms[0,:]
    
    # Significance
    ps = fperms_to_pvals(Fs, F_perms)
    
    if debug_output:
        return (ps, Fs, F_perms, perms, Gs, H2perms, IHperms, df_among, df_resid)
    else:
        return (ps, Fs, F_perms, perms)
