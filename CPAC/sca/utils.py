import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def compute_fisher_z_score(correlation_file, timeseries_one_d):

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

    roi_numbers = []
    if '#' in open(timeseries_one_d, 'r').readline().rstrip('\r\n'):
        roi_numbers = open(timeseries_one_d, 'r').readline().rstrip('\r\n').replace('#', '').split('\t')

    corr_img = nb.load(correlation_file)
    corr_data = corr_img.get_data()

    hdr = corr_img.get_header()

    corr_data = np.log((1 + corr_data)/(1 - corr_data)) / 2.0

    dims = corr_data.shape

    out_file = []

    if len(dims) == 5 or len(roi_numbers) > 0:

        if len(dims) == 5:
            x, y, z, one, roi_number = dims

            corr_data = np.reshape(corr_data, (x*y*z, roi_number), order='F')


        for i in range(0, len(roi_numbers)):

            sub_data = corr_data
            if len(dims) == 5:
                sub_data = np.reshape(corr_data[:, i], (x, y, z), order='F')

            sub_img = nb.Nifti1Image(sub_data, header=corr_img.get_header(), affine=corr_img.get_affine())

            sub_z_score_file = os.path.join(os.getcwd(), 'z_score_ROI_number_%s.nii.gz' % (roi_numbers[i]))

            sub_img.to_filename(sub_z_score_file)

            out_file.append(sub_z_score_file)

    else:

        z_score_img = nb.Nifti1Image(corr_data, header=hdr, affine=corr_img.get_affine())

        z_score_file = os.path.join(os.getcwd(), 'z_score.nii.gz')

        z_score_img.to_filename(z_score_file)

        out_file.append(z_score_file)


    return out_file


def map_to_roi(timeseries, maps):
    """
    Method to map roi numbers to correct
    component maps when split,
    """
    
    try:
        import os
        import numpy as np
        testMat = np.loadtxt(timeseries)
        timepoints, rois = testMat.shape
        
        if rois > timepoints:
            raise Exception('The number of timepoints is smaller than the'
                            + ' number of ROIs to run - therefore the'
                            + ' GLM is underspecified and can\'t run.')
        labels = open(timeseries, 'r').readline().split()
        labels = [os.path.join(os.path.dirname(maps[0]),'z_maps_roi_' + str(l)) \
                  for l in labels]
        maps.sort()
        #output = zip(labels, component_maps)
        print "labels, maps", labels, maps
    except Exception:
        print "Error while mapping roi to dual regression split 3d volumes"
    return labels, maps

