import numpy as np

def norm_cols(X):
    """
    Theano expression which centers and normalizes columns of X `||x_i|| = 1`
    """
    Xc = X - X.mean(0)
    return Xc/np.sqrt( (Xc**2.).sum(0) )



def calc_cwas(subjects_data, regressor, iter):
    """
    Performs Connectome-Wide Association Studies (CWAS) [1]_ for every voxel.  Implementation based on
    [2]_.
    
    Parameters
    ----------
    subjects_data : ndarray
        Numpy data array of shape (`S`,`T`,`V`), `S` subjects, `T` timepoints, `V` voxels.  The number
        of timepoints `T` can vary between subjects.
    regressor : ndarray
        Matrix of shape (`S`, `R`), `S` subjects and `R` regressors
    iter : integer
        Number of permutations to derive significance tests

    Returns
    -------
    F_set : ndarray
        Pseudo-F statistic calculated for every voxel
    p_set : ndarray
        Significance probabilities of F_set based on permutation tests
    
    Notes
    -----
    The distance matrix can potentially take up a great deal of memory and therefore is not
    returned.
    
    References
    ----------
    .. [1] Shehzad Z, Reiss PT, Adelstein JS, Emerson JW, Chabernaud C, Mennes M, Di Martino A, Kelly C, Castellanos FX, Milham MP. (June 2011). Connectome-Wide Association Studies (CWAS): A Multivariate Distance-Based Approach. Poster to be presented at the Annual Meeting of the Organization for Human Brain Mapping, Quebec City.
    .. [2] Xiao-Wei Song, Zhang-Ye Dong, Xiang-Yu Long, Su-Fang Li, Xi-Nian Zuo, Chao-Zhe Zhu, Yong He, Chao-Gan Yan, Yu-Feng Zang. (2011) REST: A Toolkit for Resting-State Functional Magnetic Resonance Imaging Data Processing. PLoS ONE 6(9): e25031. doi:10.1371/journal.pone.0025031
    
    """

    nSubjects = subjects_data.shape[0]
    nVoxels = subjects_data[0].shape[1]
    #Number of timepoints may be consistent between subjects

    subjects_data_n = np.zeros_like(subjects_data)
    for i in range(nSubjects):
        subjects_data_n[i] = norm_cols(subjects_data[i])

    # Distance matrices for every voxel
    D = np.zeros((nVoxels, nSubjects, nSubjects))
    
    F_set = np.zeros(nVoxels)
    p_set = np.zeros(nVoxels)
    
    # For a particular voxel v, its spatial correlation map for every subject
    S = np.zeros((nSubjects, nVoxels))
    
    for i in range(nVoxels):
        for i_s in range(nSubjects):
            # Correlate voxel i with every other voxel in the same subject
            S[i_s,:] = subjects_data_n[i_s][:,i][np.newaxis, :].dot(subjects_data_n[i_s])
        # Remove auto-correlation column to prevent infinity in Fischer z transformation
        S0 = np.delete(S,i,1)
        
        S0 = 0.5*np.log((1+S0)/(1-S0))
        
        #Normalize the rows
        S0 = norm_cols(S0.T).T
        D[i,:,:] = 1-S0.dot(S0.T)
        F_set[i], p_set[i] = y_mdmr(D[i], regressor, iter)
    
    return F_set, p_set

def y_mdmr(yDis, x, iter):
    """
    Multivariate Distance Matrix Regression
    
    Parameters
    ----------
    yDis : ndarray
    x : ndarray
    iter : integer
    
    Returns
    --------
    F : float
    p : float
    
    Notes
    -----
    Implementation based on work by YAN Chao-Gan (ycg.yan@gmail.com) and References section.
    
    References
    -----------
    .. [1] Anderson, M. J. 2002. DISTML v.2: a FORTRAN computer program to calculate a distance-based multivariate analysis for a linear model. Dept. of Statistics University of Auckland. (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
    .. [2] Anderson, M. J. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology 26: 32-46.
    .. [3] Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed. Elsevier Science BV, Amsterdam.
    .. [4] McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to community data: a comment on distance-based redundancy analysis. Ecology 290-297.
    .. [5] Neter, J., M. H. Kutner, C. J. Nachtsheim, and W. Wasserman. 1996. Applied linear statistical models. 4th ed. Irwin, Chicago, Illinois.
    """
    n = yDis.shape[0]
    A = -0.5*(yDis**2)
    I = np.eye(n,n)
    uno = np.ones((n,1))
    C = I - (1.0/n)*uno.dot(uno.T)
    
    # G is similar to matrix of inner products from distances used in Partha Niyogi's
    # multidimensional scaling

    G = C.dot(A).dot(C)
    xx = np.hstack((uno, x))
    Q1, R1 = np.linalg.qr(xx)

    H = Q1.dot(Q1.T)
    m = xx.shape[1]
    
    df_among = m-1
    df_resid = n-m
    df_total = n-1
    
    SS_among = np.trace(H.dot(G).dot(H))
    SS_resid = np.trace((I-H).dot(G).dot(I-H))
    
    SS_total = np.trace(G)
    
    MS_among = SS_among/df_among
    MS_resid = SS_resid/df_resid
    
    F = MS_among/MS_resid
    
    if iter > 0:
        F_perm = np.zeros(iter-1)
        for i in range(iter-1):
            IndexRandPerm = np.random.permutation(n)
            # Would it be better to permute the regressors x rather than G?
            G_perm = G.take(IndexRandPerm, axis=0)
            G_perm = G_perm.take(IndexRandPerm, axis=1)
            MS_perm = np.trace(H.dot(G_perm).dot(H))/df_among
            MSE_perm = np.trace((I-H).dot(G_perm).dot(I-H))/df_resid
            F_perm[i] = MS_perm/MSE_perm
        j = (F_perm >= F).sum().astype('float')
        p = (j+1.0)/iter
    else:
        print 'No permutation tests done.'
        p = np.NAN
    
    return F, p