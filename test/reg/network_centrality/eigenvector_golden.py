# eigenvector_golden.py
#
# Author: Daniel Clark, 2016

'''
Calculate eigenvector centrality of nifti image without block sizes
or memory-limiting methods
'''


def eigen_centrality(nii_path, mask_path, thresh_type, thresh_val):
    '''
    Function to compute the eigenvector centrality of a functional
    image without memory-limiting (RAM) considerations

    Parameters
    ----------
    nii_path : string
        filepath to the nifti functional file
    mask_path : string
        filepath to the nifti mask file
    thresh_type : string
        type of thresholding (either 'correlation' or 'sparsity')
    thresh_val : float
        threshold to cutoff in similarity matrix
    '''

    # Import packages
    import nibabel as nib
    import numpy as np
    from scipy.sparse import linalg

    # Init variables
    data = nib.load(nii_path)
    data_arr = data.get_data().astype('float32')
    data_aff = data.get_affine()
    mask = nib.load(mask_path)
    mask_arr = mask.get_data().astype('bool')

    # Get data mask where no variance
    datmask = data_arr.var(axis=3).astype('bool')
    mask_arr = mask_arr & datmask

    # Extract time series (V x T)
    time_series = data_arr[mask_arr]

    # Transpose (numpy math likes T x V for de-meaning) and normalize
    time_series = time_series.T
    ts_demeaned = time_series - time_series.mean(0)
    ts_normd = ts_demeaned / np.sqrt((ts_demeaned**2.0).sum(0))
    # Get info from timeseries
    num_tpts = ts_normd.shape[0]
    num_voxs = ts_normd.shape[1]
    num_conns = (num_voxs**2-num_voxs)/2.0

    # Calculate similarity matrix and threshold
    sim_mat = np.dot(ts_normd.T, ts_normd)
    if thresh_type == 'sparsity':
        thresh_idx = int(thresh_val*num_conns)
        uptri = np.triu(sim_mat, k=1)
        sort_arr = sim_mat[np.where(uptri)]
        sort_arr.sort()
        thresh_val = sort_arr[-thresh_idx]

    # Threshold similarity matrix
    mat_mask = sim_mat >= thresh_val
    bin_mat = mat_mask.astype('float32')
    sim_mat[np.logical_not(mat_mask)] = 0

    # Calculate eigenvectors
    eigen_val, bin_eigen_vec = linalg.eigsh(bin_mat, k=1, which='LM', maxiter=1000)
    eigen_val, wght_eigen_vec = linalg.eigsh(sim_mat, k=1, which='LM', maxiter=1000)

    # Map eigenvector back to nifti
    coords = np.argwhere(mask_arr)
    bin_out_arr = np.zeros(mask_arr.shape)
    wght_out_arr = np.zeros(mask_arr.shape)
    for idx, xyz in enumerate(coords):
        x, y, z = xyz
        bin_out_arr[x, y, z] = bin_eigen_vec[idx]
        wght_out_arr[x, y, z] = wght_eigen_vec[idx]

    # Out nifti images
    bin_out_nii = nib.Nifti1Image(bin_out_arr, data_aff)
    wght_out_nii = nib.Nifti1Image(wght_out_arr, data_aff)
    bin_out_nii.to_filename('eigenvector_%s_binarize.nii.gz' % thresh_type)
    wght_out_nii.to_filename('eigenvector_%s_weighted.nii.gz' % thresh_type)


if __name__ == '__main__':
    eigen_centrality('/home/dclark/work-dir/residual_antswarp_flirt.nii.gz',
                     '/home/dclark/work-dir/benchmark_centrality_mask.nii.gz',
                     'correlation', 0.6)
