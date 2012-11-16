import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def compute_fisher_z_score(correlation_file):

    """
    Computes the fisher z transform of the input correlation map
    If the correlation map contains data for multiple ROIs then 
    the function returns z score for each ROI as a seperate nifti 
    file


    Parameters
    ----------

    correlation_file: string
        Input correlations file
    

    Returns
    -------

    out_file : list (nifti files)
        list of z_scores for mask or ROI
    """

    import nibabel as nb
    import numpy as np
    import os

    corr_img = nb.load(correlation_file)
    corr_data = corr_img.get_data()

    hdr = corr_img.get_header()

    corr_data = np.log((1 + corr_data)/(1 - corr_data)) / 2.0

    dims = corr_data.shape

    out_file = []

    if len(dims) == 5:

        x, y, z, one, roi_number = dims

        corr_data = np.reshape(corr_data, (x*y*z, roi_number), order='F')

        for i in range(0, roi_number):

            sub_data = np.reshape(corr_data[:, i], (x, y, z), order='F')

            sub_img = nb.Nifti1Image(sub_data, header=corr_img.get_header(), affine=corr_img.get_affine())

            sub_z_score_file = os.path.join(os.getcwd(), 'z_score_ROI_%d.nii.gz' % (i+1))

            sub_img.to_filename(sub_z_score_file)

            out_file.append(sub_z_score_file)

    else:

        z_score_img = nb.Nifti1Image(corr_data, header=hdr, affine=corr_img.get_affine())

        z_score_file = os.path.join(os.getcwd(), 'z_score.nii.gz')

        z_score_img.to_filename(z_score_file)

        out_file.append(z_score_file)


    return out_file
