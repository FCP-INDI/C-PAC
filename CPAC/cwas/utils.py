import numpy as np
from mdmr import *
from subdist import *

def calc_cwas(subjects_data, regressor, cols, iter, voxel_range, strata=None):
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
    voxel_range : tuple
        (start, end) tuple specify the range of voxels (inside the mask) to perform cwas on.    
    strata : None or list
        todo
        
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
    
    D            = calc_subdists(subjects_data, voxel_range)
    F_set, p_set = calc_mdmrs(D, regressor, cols, iter, strata)
    
    return F_set, p_set

def calc_subdists(subjects_data, voxel_range):
    nSubjects   = len(subjects_data)
    vox_inds    = range(*voxel_range)
    nVoxels     = len(vox_inds)
    #Number of timepoints may be consistent between subjects
    
    subjects_normed_data = norm_subjects(subjects_data)
    
    # Distance matrices for every voxel
    D = np.zeros((nVoxels, nSubjects, nSubjects))
        
    # For a particular voxel v, its spatial correlation map for every subject
    S = np.zeros((nSubjects, 1, nVoxels))
    
    for i in range(nVoxels):
        S    = ncor_subjects(subjects_normed_data, [vox_inds[i]])
        S0   = np.delete(S[:,0,:], vox_inds[i], 1)    # remove autocorrelations
        S0   = fischers_transform(S0)
        D[i] = compute_distances(S0)
    
    return D

def calc_mdmrs(D, regressor, cols, iter, strata=None):
    nVoxels = D.shape[0]
    nSubjects = D.shape[1]
    
    F_set = np.zeros(nVoxels)
    p_set = np.zeros(nVoxels)
    
    for i in range(nVoxels):
        p_set[i], F_set[i], _, _ = mdmr(D[i].reshape(nSubjects**2,1), regressor, cols, iter, strata)
    
    return F_set, p_set
