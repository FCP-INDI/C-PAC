import numpy as np


def calc_compcor_components(data, nComponents, wm_sigs, csf_sigs):

    import scipy.signal as signal
    
    wmcsf_sigs = np.vstack((wm_sigs, csf_sigs))

    # filter out any voxels whose variance equals 0
    print 'Removing zero variance components'
    wmcsf_sigs = wmcsf_sigs[wmcsf_sigs.std(1)!=0,:]

    if wmcsf_sigs.shape.count(0):
        err = "\n\n[!] No wm or csf signals left after removing those " \
              "with zero variance.\n\n"
        raise Exception(err)
    
    print 'Detrending and centering data'
    Y = signal.detrend(wmcsf_sigs, axis=1, type='linear').T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    Yc = Yc / np.tile(np.array(Y.std(0)).reshape(1,Y.shape[1]), (Y.shape[0],1))
    
    print 'Calculating SVD decomposition of Y*Y\''
    U, S, Vh = np.linalg.svd(Yc)
    
    return U[:, :nComponents]


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


def create_despike_regressor_matrix(frames_excluded, total_vols):
    """Create a Numpy array describing which volumes are to be regressed out
    during nuisance regression, for de-spiking.

    :param frames_excluded: 1D file of the volume indices to be excluded. This
    is a 1D text file of integers separated by commas.
    :param total_vols: integer value of the length of the time series (number
    of volumes).
    :return: Numpy array consisting of a row for every volume, and a column
    for every volume being regressed out, with a 1 where they match.
    """

    with open(frames_excluded, 'r') as f:
        excl_vols = f.readlines()

    if len(excl_vols) > 0:
        excl_vols = sorted([int(x) for x in excl_vols[0].split(',') if x != ''])
    else:
        return None

    reg_matrix = np.zeros((total_vols, len(excl_vols)), dtype=int)

    i = 0
    for vol in excl_vols:
        reg_matrix[vol][i] = 1
        i += 1

    return reg_matrix
