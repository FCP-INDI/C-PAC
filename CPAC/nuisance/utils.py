import numpy as np

def calc_compcor_components(data, nComponents, wm_sigs, csf_sigs):
    import scipy.signal as signal
    
    wmcsf_sigs = np.vstack((wm_sigs, csf_sigs))

    # filter out any voxels whose variance equals 0
    print 'Removing zero variance components'
    wmcsf_sigs = wmcsf_sigs[wmcsf_sigs.std(1)!=0,:]

    if wmcsf_sigs.shape.count(0):
        print 'No wm or csf signals left after removing those with zero variance'
        raise IndexError
    
    print 'Detrending and centering data'
    Y = signal.detrend(wmcsf_sigs, axis=1, type='linear').T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    Yc = Yc / np.tile(np.array(Y.std(0)).reshape(1,Y.shape[1]), (Y.shape[0],1))
    
    print 'Calculating SVD decomposition of Y*Y\''
    U, S, Vh = np.linalg.svd(Yc)
    
    return U[:,:nComponents]

def erode_mask(data):
    mask = data != 0
    eroded_mask = np.zeros_like(data, dtype='bool')
    max_x, max_y, max_z = data.shape
    x,y,z = np.where(data != 0)
    for i in range(x.shape[0]):
        if (max_x-1) == x[i] or \
           (max_y-1) == y[i] or \
           (max_z-1) == z[i] or \
           x[i] == 0 or \
           y[i] == 0 or \
           z[i] == 0:
            eroded_mask[x[i],y[i],z[i]] = False
        else:
            eroded_mask[x[i],y[i],z[i]] = mask[x[i], y[i], z[i]] * \
                                          mask[x[i] + 1, y[i], z[i]] * \
                                          mask[x[i], y[i] + 1, z[i]] * \
                                          mask[x[i], y[i], z[i] + 1] * \
                                          mask[x[i] - 1, y[i], z[i]] * \
                                          mask[x[i], y[i] - 1, z[i]] * \
                                          mask[x[i], y[i], z[i] - 1]

    eroded_data = np.zeros_like(data)
    eroded_data[eroded_mask] = data[eroded_mask]
    
    return eroded_data
