import os
import scipy.signal as signal
import nibabel as nb
import numpy as np
from CPAC.utils import safe_shape

def calc_compcor_components(data_filename, num_components, mask_filename):

    if num_components < 1:
        raise ValueError('Improper value for num_components ({0}), should be >= 1.'.format(num_components))

    try:
        image_data = nb.load(data_filename).get_data().astype(np.float64)
    except:
        print('Unable to load data from {0}'.format(data_filename))
        raise

    try:
        binary_mask = nb.load(mask_filename).get_data().astype(np.int16)
    except:
        print('Unable to load data from {0}'.format(mask_filename))

    if not safe_shape(image_data, binary_mask):
        raise ValueError('The data in {0} and {1} do not have a consistent shape'.format(data_filename, mask_filename))
    
    # make sure that the values in binary_mask are binary
    binary_mask[binary_mask > 0] = 1
    binary_mask[binary_mask != 1] = 0

    # reduce the image data to only the voxels in the binary mask
    image_data = image_data[binary_mask==1, :]
    
    # filter out any voxels whose variance equals 0
    print 'Removing zero variance components'
    image_data = image_data[image_data.std(1)!=0,:]

    if image_data.shape.count(0):
        err = "\n\n[!] No wm or csf signals left after removing those " \
              "with zero variance.\n\n"
        raise Exception(err)
    
    print 'Detrending and centering data'
    Y = signal.detrend(image_data, axis=1, type='linear').T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    Yc = Yc / np.tile(np.array(Yc.std(0)).reshape(1,Yc.shape[1]), (Yc.shape[0],1))
    
    print 'Calculating SVD decomposition of Y*Y\''
    U, S, Vh = np.linalg.svd(Yc)

    # write out the resulting regressor file
    regressor_file = os.path.join(os.getcwd(), 'compcor_regressors.1D')
    np.savetxt(regressor_file, U[:, :num_components], delimiter='\t', fmt='%16g')

    return regressor_file
