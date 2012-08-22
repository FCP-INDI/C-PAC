import numpy as np

def calc_compcor_components(data, nComponents, wm_mask, csf_mask):
    import scipy
    
    wmcsf_mask = (csf_mask + wm_mask).astype('bool')
    
    print 'Detrending and centering data'
    Y = scipy.signal.detrend(data[wmcsf_mask], axis=1, type='linear').T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    
    print 'Calculating SVD decomposition of Y*Y\''
    U, S, Vh = np.linalg.svd(np.dot(Yc, Yc.T))
    
    return U[:,:nComponents]
