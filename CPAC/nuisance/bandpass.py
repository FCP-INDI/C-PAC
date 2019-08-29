import os
import numpy as np
import nibabel as nb

from scipy.fftpack import fft, ifft


def bandpass_voxels(realigned_file, bandpass_freqs, sample_period=None):
    """
    Performs ideal bandpass filtering on each voxel time-series.
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    bandpass_freqs : tuple
        Tuple containing the bandpass frequencies. (LowCutoff_HighPass HighCutoff_LowPass)
    sample_period : float, optional
        Length of sampling period in seconds.  If not specified,
        this value is read from the nifti file provided.
        
    Returns
    -------
    bandpassed_file : string
        Path of filtered output (nifti file).
    
    """

    def ideal_bandpass(data, sample_period, bandpass_freqs):
        # Derived from YAN Chao-Gan 120504 based on REST.
        sample_freq = 1. / sample_period
        sample_length = data.shape[0]

        data_p = np.zeros(int(2**np.ceil(np.log2(sample_length))))
        data_p[:sample_length] = data

        LowCutoff, HighCutoff = bandpass_freqs

        if (LowCutoff is None):  # No lower cutoff (low-pass filter)
            low_cutoff_i = 0
        elif (LowCutoff > sample_freq / 2.):
            # Cutoff beyond fs/2 (all-stop filter)
            low_cutoff_i = int(data_p.shape[0] / 2)
        else:
            low_cutoff_i = np.ceil(
                LowCutoff * data_p.shape[0] * sample_period).astype('int')

        if (HighCutoff > sample_freq / 2. or HighCutoff is None):
            # Cutoff beyond fs/2 or unspecified (become a highpass filter)
            high_cutoff_i = int(data_p.shape[0] / 2)
        else:
            high_cutoff_i = np.fix(
                HighCutoff * data_p.shape[0] * sample_period).astype('int')

        freq_mask = np.zeros_like(data_p, dtype='bool')
        freq_mask[low_cutoff_i:high_cutoff_i + 1] = True
        freq_mask[
            data_p.shape[0] -
            high_cutoff_i:data_p.shape[0] + 1 - low_cutoff_i
        ] = True

        f_data = fft(data_p)
        f_data[freq_mask != True] = 0.
        data_bp = np.real_if_close(ifft(f_data)[:sample_length])

        return data_bp

    nii = nb.load(realigned_file)
    data = nii.get_data().astype('float64')
    mask = (data != 0).sum(-1) != 0
    Y = data[mask].T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))

    if not sample_period:
        hdr = nii.get_header()
        sample_period = float(hdr.get_zooms()[3])
        # Sketchy check to convert TRs in millisecond units
        if sample_period > 20.0:
            sample_period /= 1000.0

    Y_bp = np.zeros_like(Y)
    for j in range(Y.shape[1]):
        Y_bp[:, j] = ideal_bandpass(Yc[:, j], sample_period, bandpass_freqs)

    data[mask] = Y_bp.T
    img = nb.Nifti1Image(data, header=nii.get_header(),
                         affine=nii.get_affine())
    bandpassed_file = os.path.join(os.getcwd(),
                                   'bandpassed_demeaned_filtered.nii.gz')
    img.to_filename(bandpassed_file)

    return bandpassed_file
