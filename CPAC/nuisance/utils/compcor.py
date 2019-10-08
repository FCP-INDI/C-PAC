

def calc_compcor_components(data_file, nComponents, mask_file):

    import scipy.signal as signal
    import nibabel as nb
    from CPAC.utils import safe_shape

    try:
        nii = nb.load(data_file).get_data().astype(np.float64)
    except:
        print('Unable to load data from {0}'.format(data_file)
        raise()
    
    # Check and define regressors which are provided from files
    if mask_file is not None:
        binary_mask = nb.load(mask_file).get_data().astype(np.float64) # change to int16
        # check the dimensions

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

