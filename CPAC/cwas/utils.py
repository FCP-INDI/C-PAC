import numpy as np
from mdmr import *
from subdist import *
from memlimit import *

def calc_cwas(subjects_data, regressor, cols, nperms, voxel_range, memlimit=4, strata=None, dtype='float64'):
    """
    Performs Connectome-Wide Association Studies (CWAS) [1]_ for every voxel.  For similar implementations see 
    [2]_ and [3]_.
    
    Parameters
    ----------
    subjects_data : ndarray or list
        Numpy data array of shape (`S`,`T`,`V`), `S` subjects, `T` timepoints, `V` voxels.  The number
        of timepoints `T` can vary between subjects in which case data is a list of length `S` with each
        element being a numpy array of shape (`T`, `V`).
    regressor : ndarray
        Matrix of shape (`S`, `R`), `S` subjects and `R` regressors
    nperms : integer
        Number of permutations to derive significance tests
    voxel_range : tuple
        (start, end) tuple specify the range of voxels (inside the mask) to perform cwas on.    
    memlimit : int
        TODO
    strata : None or list
        TODO
    dtype : str
        TODO
        
    Returns
    -------
    F_set : ndarray
        Pseudo-F statistic calculated for every voxel
    p_set : ndarray
        Significance probabilities of F_set based on permutation tests
    Fperms : ndarray
        TODO
    
    Notes
    -----
    The distance matrix can potentially take up a great deal of memory and therefore is not
    returned.
    
    References
    ----------
    .. [1] Shehzad Z, Reiss PT, Adelstein JS, Emerson JW, Chabernaud C, Mennes M, Di Martino A, Kelly C, Castellanos FX, Milham MP. (June 2011). Connectome-Wide Association Studies (CWAS): A Multivariate Distance-Based Approach. Poster to be presented at the Annual Meeting of the Organization for Human Brain Mapping, Quebec City.
    .. [2] Xiao-Wei Song, Zhang-Ye Dong, Xiang-Yu Long, Su-Fang Li, Xi-Nian Zuo, Chao-Zhe Zhu, Yong He, Chao-Gan Yan, Yu-Feng Zang. (2011) REST: A Toolkit for Resting-State Functional Magnetic Resonance Imaging Data Processing. PLoS ONE 6(9): e25031. doi:10.1371/journal.pone.0025031
    .. [3] Shehzad Z. connectir R package.
    
    """
    
    # TODO:
    # - add a test for varying the dtype?
    voxel_block = voxel_blocks_for_subdists(memlimit, subjects_data, dtype)
    D           = calc_subdists(subjects_data, voxel_range, voxel_block, dtype)
    
    # TODO:
    # - calculate the voxel block for mdmr but need to edit mdmr to give you that blessing
    # - also pass on the dtype variable
    voxel_block = voxel_blocks_for_mdmr(memlimit, D, nperms, dtype)
    F_set, p_set, Fperms = calc_mdmrs(D, regressor, cols, nperms, strata, voxel_block)
    
    return F_set, p_set, Fperms


###
# Computing Distances
###

def voxel_blocks_for_subdists(memlimit, subjects_data, dtype):
    """docstring for voxel_blocks_for_subdists"""
    nsubjs  = len(subjects_data)
    nvoxs   = subjects_data[0].shape[0]
    ntpts   = [ sdata.shape[1] for sdata in subjects_data ]
    
    # First determine the things with fixed requirements
    mem_func_data   = subjs_func_memory(nvoxs, ntpts, dtype)
    mem_dmat        = dmat_memory(nvoxs, nsubjs, dtype)
    
    # Second determine the flexible things in life (# of voxels, of course)
    resid_memlimit  = memlimit - mem_func_data - mem_dmat
    nvoxs           = nvoxs_with_conn_map(resid_memlimit, nvoxs, nsubjs, dtype)
    
    # TODO: better Exception name?
    if nvoxs == 0:
        raise Exception("Memory limit of %i leaves no voxels to process" % memlimit)
    
    return nvoxs

def split_list_into_groups(l, n):
    """
    l: your list
    n: number of elements within each group when splitting the list
    """
    # Borrowed from Steg at 
    # http://stackoverflow.com/questions/1624883/alternative-way-to-split-a-list-into-groups-of-n
    return [ l[i:i+n] for i in range(0, len(l), n) ]    

def calc_subdists(subjects_data, voxel_range, voxel_block=1, dtype='float64', dfun=compute_distances):
    """
    Calculates the distances at each voxel between connectivity patterns of 
    all possible pairs of participants.
    
    Parameters
    ----------
    subjects_data : ndarray
        Numpy data array of shape (`S`,`T`,`V`), `S` subjects, `T` timepoints, `V` voxels.  The number
        of timepoints `T` can vary between subjects in which case data is a list of length `S` with each
        element being a numpy array of shape (`T`, `V`).
    voxel_range : tuple
        (start, end) tuple specify the range of voxels (inside the mask) to perform cwas on.
    voxel_block : int (default = 1)
        Voxels connectivity patterns to be calculated at once for all participants.
    dtype : str
        TODO
    dfun : function (default=compute_distances)
        Function used to compute the distances. Default function computes pearson correlation between 
        participant connectivity patterns. This function must take a `S` x `V` matrix of connectivity 
        patterns and return a `S` x `S` distance matrix.
    
    Returns
    -------
    D : ndarray
        Distance matrix calculated for every voxel (`V` voxels x `S` subjects x `S` subjects)
    
    """
    nSubjects   = len(subjects_data)
    nVoxels     = len(range(*voxel_range))
    
    vox_inds    = split_list_into_groups(range(*voxel_range), voxel_block)
    subvox_inds = split_list_into_groups(range(nVoxels), voxel_block)
    nGroups     = len(vox_inds)
    
    # Norm the subject's functional data time-series
    subjects_normed_data = norm_subjects(subjects_data)
    
    # Distance matrices for every voxel
    D = np.zeros((nVoxels, nSubjects, nSubjects), dtype=dtype)
    
    for i in range(nGroups):
        S    = ncor_subjects(subjects_normed_data, vox_inds[i])
        S0   = replace_autocorrelations(S)
        S0   = fischers_transform(S0)
        for ij,j in enumerate(subvox_inds[i]):
            D[j] = compute_distances(S0[:,ij,:])
    
    return D

def classic_calc_subdists(subjects_data, voxel_range):
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


###
# Computing MDMR
###

def voxel_blocks_for_mdmr(memlimit, D, nperms, dtype):
    """docstring for voxel_blocks_for_mdmr"""
    nvoxs   = D.shape[0]
    nsubjs  = D.shape[1]
    
    # First determine the things with fixed requirements
    mem_dmat        = dmat_memory(nvoxs, nsubjs, dtype)
    mem_perms       = perm_mats_memory(nperms, nsubjs, dtype) # includes H2/IH
    mem_fperms      = fperms_memory(nperms, nvoxs, dtype)
    
    # Second determine the flexible things in life (# of voxels, of course)
    resid_memlimit  = memlimit - mem_dmat - mem_perms - mem_fperms
    nvoxs           = nvoxs_per_mdmr(memlimit, nperms, nsubjs, nvoxs, dtype)
    
    # TODO: better Exception name?
    if nvoxs == 0:
        raise Exception("Memory limit of %i leaves no voxels to process, sadness" % memlimit)
    
    return nvoxs

def calc_mdmrs(D, regressor, cols, perms, strata=None, voxel_block=1):
    nVoxels     = D.shape[0]
    nSubjects   = D.shape[1]
    vox_inds    = split_list_into_groups(range(nVoxels), voxel_block)
    
    perms, H2perms, IHperms = gen_perm_mats(regressor, cols, perms, strata)
    nperms      = perms.shape[0]
    
    Fs          = np.zeros(nVoxels)
    ps          = np.zeros(nVoxels)
    Fperms      = np.zeros((nperms, nVoxels))
    
    for i,inds in enumerate(vox_inds):
        ps[inds], Fs[inds], Fperms[:,inds], _ = mdmr(D[inds], regressor, cols, 
                                                     perms, strata, H2perms, 
                                                     IHperms)
    
    return Fs, ps, Fperms


def classic_calc_mdmrs(D, regressor, cols, iter, strata=None):
    nVoxels = D.shape[0]
    nSubjects = D.shape[1]
    
    F_set = np.zeros(nVoxels)
    p_set = np.zeros(nVoxels)
    
    for i in range(nVoxels):
        p_set[i], F_set[i], _, _ = mdmr(D[i].reshape(nSubjects**2,1), regressor, cols, iter, strata)
    
    return F_set, p_set
