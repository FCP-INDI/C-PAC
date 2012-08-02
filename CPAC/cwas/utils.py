import theano
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

def calc_distance_matrix(subjects_data):
    
    # Compile theano functions
    X,Y = T.dmatrices('X','Y')
    tdot = theano.function([X,Y], T.dot(X,Y))
    tnormcols = theano.function([X], TnormCols(X))

    nSubjects = subjects_data.shape[0]
    nTimepoints = subjects_data.shape[1]
    nVoxels = subjects_data.shape[2]
    
    subjects_data_n = np.zeros_like(subjects_data)
    for i in range(nSubjects):
        subjects_data_n[i] = tnormcols(subjects_data[i])
    
    # Distance matrices for every voxel
    D = np.zeros((nVoxels, nSubjects, nSubjects))
    
    # For a particular voxel v, its spatial correlation map for every subject
    S = np.zeros((nSubjects, nVoxels))
    
    for i in range(nVoxels):
        for i_s in range(nSubjects):
            # Correlate voxel i with every other voxel in the same subject
            S[i_s,:] = tdot(subjects_data_n[i_s,:,i][np.newaxis, :], subjects_data_n[i_s])
        # Remove auto-correlation column to prevent infinity in Fischer z transformation
        S0 = np.delete(S,i,1)
        
        S0 = 0.5*np.log((1+S0)/(1-S0))
        
        #Normalize the rows
        S0 = tnormcols(S0.T).T
        D[i,:,:] = 1-tdot(S0,S0.T)
    
    return D

def y_mdmr(yDis, x, iter):
    n = yDis.shape[0]
    A = -0.5*(yDis**2)
    I = np.eye(n,n)
    uno = np.ones((n,1))
    C = I - (1.0/n)*np.dot(uno,uno.T)
    G = np.dot(np.dot(C, A),C)
    
    np.linalg.qr()