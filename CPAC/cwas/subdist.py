import numpy as np

# TODO:
# - less transposes when computing the distances?
# - better approach to fix autocorrelation issue

def norm_cols(X):
    """
    Theano expression which centers and normalizes columns of X `||x_i|| = 1`
    """
    Xc = X - X.mean(0)
    return Xc/np.sqrt( (Xc**2.).sum(0) )

def norm_subjects(subjects_data):
    nSubjects = len(subjects_data)
    subjects_data_n = [None] * nSubjects
    for i in range(nSubjects):
        subjects_data_n[i] = norm_cols(subjects_data[i])
    return subjects_data_n

# add a possible 2nd column?
def ncor(normed_data, vox_inds):
    #if type(vox_inds) != list or type(vox_inds) != tuple:
    #    vox_inds = [vox_inds]
    return normed_data[:,vox_inds].T.dot(normed_data)

def ncor_subjects(subjects_normed_data, vox_inds):
    nSubjects = len(subjects_normed_data)
    nVoxels   = subjects_normed_data[0].shape[1]
    nSeeds    = len(vox_inds)
    
    S         = np.zeros((nSubjects, nSeeds, nVoxels))
    for i in range(nSubjects):
        S[i] = ncor(subjects_normed_data[i], vox_inds)
    
    ## Prevent infinity for Fischer's Tranfrom of autocorrelation (1s)
    #for j in range(nSeeds):
    #    S[:,j,vox_inds[j]] = S[:,j,vox_inds[j]] - 1e-9
    
    return S

def fischers_transform(S):
    return np.arctanh(S)

def compute_distances(S0):
    S0   = norm_cols(S0.T).T
    dmat = 1 - S0.dot(S0.T)
    return dmat
