import numpy as np
from mdmr import *

def norm_cols(X):
    """
    Theano expression which centers and normalizes columns of X `||x_i|| = 1`
    """
    Xc = X - X.mean(0)
    return Xc/np.sqrt( (Xc**2.).sum(0) )



def calc_cwas(subjects_data, regressor, cols, iter, strata=None):
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
    nVoxels = subjects_data.shape[2]
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
        
        p_set[i], F_set[i], _, _ = mdmr(D[i].reshape(nSubjects**2,1), regressor, cols, iter, strata)
    
    return F_set, p_set
