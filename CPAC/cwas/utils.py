import theano.tensor as T
import numpy as np

def TnormCols(X):
    """
    Theano expression which centers and normalizes columns of X `||x_i|| = 1`
    """
    Xc = X - X.mean(0)
    return Xc/T.sqrt( (Xc**2.).sum(0) )

def TzscorrCols(Xn):
    """
    Theano expression which returns Fisher transformed correlation values between columns of a
    normalized input, `X_n`.  Diagonal is set to zero.
    """
    C_X = T.dot(Xn.T, Xn)-T.eye(Xn.shape[1])
    return 0.5*T.log((1+C_X)/(1-C_X))

def calc_distance_matrix():
    R = np.zeros((nVoxels, nSubjects, nSubjects))
    S = np.zeros((nSubjects, nVoxels))
    
    for i in range(nVoxels):
        for i_s in range(subjects_data.shape[0]):
            S[i_s,:] = tdot(subjects_data[i_s,:,i][np.newaxis, :], subjects_data[i_s])
        S[:, i] = 0
        S = 0.5*np.log((1+S)/(1-S))
        S = tnormcols(S.T).T
        R[i,:,:] = tdot(S,S.T)