# coding: utf-8
import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from CPAC.reho.utils import *


def create_reho():

    """
    Regional Homogeneity(ReHo) approach to fMRI data analysis

    This workflow computes the ReHo map, z-score on map, transforms z_score in MNI space and smooths it

    Parameters
    ----------

    None

    Returns
    -------
    reHo : workflow
        Regional Homogeneity Workflow

    Notes
    -----

    `Source <https://github.com/openconnectome/C-PAC/blob/master/CPAC/reho/reho.py>`_

    Workflow Inputs: ::

        inputspec.rest_res_filt : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.rest_mask : string (existing nifti file)
            Input Whole Brain Mask of EPI 4D Volume

        inputspec.fieldcoeff_file : string (existing nifti file)
            File with warp coefficients/fields. This typically the output given by the --cout parameter of fnirt during registration step

        inputspec.premat : string (existing affine transformation .mat file)
            Specifies an affine transform that should be applied to the data prior to the non-linear warping(example_func2highres.mat).

        inputspec.standard : string (existing nifti file)
            FSL standard nifti file in user specified resolution

        inputspec.cluster_size : integer
            For a brain voxel the number of neighbouring brain voxels to use for KCC.
            Possible values are 27, 19, 7. Recommended value 27

        fwhm_input.fwhm : list (float)
            full width half max for spatial smoothing ReHo Z scores in MNI space 


    Workflow Outputs: ::

        outputspec.raw_reho_map : string (nifti file)

        outputspec.z_score : string (nifti file)

        outputspec.z_score_non_linear_transformed : string (nifti file)

        outputspec.z_score_non_linear_transformed_smooth : string (nifti file)


    ReHo Workflow Procedure:

    1. Generate ReHo map from the input EPI 4D volume, EPI mask and cluster_size
    2. Compute Z score of the ReHo map by subtracting mean and dividing by standard deviation
    3. Transfor Z score of the ReHo map into MNI space by non-linear transform using applywarp
    4. Spatially smooth the MNI Z score of ReHo with the supplied Full Width Half Max 


    Workflow Graph:

    .. image:: ../images/reho.dot.png
        :width: 500

    Detailed Workflow Graph:

    .. image:: ../images/reho_detailed.dot.png
        :width: 500
        
    References
    ---------- 
    .. [1] Zang, Y., Jiang, T., Lu, Y., He, Y.,  Tian, L. (2004). Regional homogeneity approach to fMRI data analysis. NeuroImage, 22(1), 394, 400. doi:10.1016/j.neuroimage.2003.12.030

    Examples
    --------
    >>> from CPAC import reho
    >>> wf = reho.create_reho()
    >>> wf.inputs.inputspec.rest_res_filt = '/home/data/Project/subject/func/rest_res_filt.nii.gz'
    >>> wf.inputs.inputspec.rest_mask = '/home/data/Project/subject/func/rest_mask.nii.gz'
    >>> wf.inputs.inputspec.premat = '/home/data/Project/subject/func/reg/example_func2highres.mat'
    >>> wf.inputs.inputspec.fieldcoeff_file = '/home/data/Project/subject/anat/reg/highres2standard_warp.nii.gz'
    >>> wf.inputs.inputspec.cluster_size = 27
    >>> wf.run()
    """



    def set_gauss_reho(fwhm):

        op_string = ""

        fwhm = float(fwhm)

        sigma = float(fwhm / 2.3548)

        op = "-kernel gauss %f -fmean " % (sigma)
        op_string = op

        return op_string



    reHo = pe.Workflow(name='reHo')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'cluster_size',
                                                'premat',
                                                'rest_res_filt',
                                                'fieldcoeff_file',
                                                'rest_mask',
                                                'standard']),
                        name='inputspec')


    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')

    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'raw_reho_map',
                                                    'z_score',
                                                    'z_score_non_linear_transformed',
                                                    'z_score_non_linear_transformed_smooth']),
                        name='outputspec')

    op_string = pe.MapNode(util.Function(input_names=['mean', 'std_dev'],
                                         output_names=['op_string'],
                           function=getOpString),
                           name='op_string',
                           iterfield=['mean', 'std_dev'])

    raw_reho_map = pe.MapNode(util.Function(input_names=['in_file', 'mask_file', 'cluster_size'],
                                   output_names=['out_file'],
                     function=compute_reho),
                     name='reho_map',
                     iterfield=['in_file',
                                'mask_file'])

    mean = pe.MapNode(interface=fsl.ImageStats(),
                      name='mean',
                      iterfield=['in_file',
                                 'mask_file'])
    mean.inputs.op_string = '-k %s -m'


    standard_deviation = pe.MapNode(interface=fsl.ImageStats(),
                     name='standard_deviation',
                     iterfield=['in_file',
                                'mask_file'])
    standard_deviation.inputs.op_string = '-k %s -s'


    z_score = pe.MapNode(interface=fsl.MultiImageMaths(),
                   name='z_score',
                   iterfield=['in_file',
                              'operand_files',
                              'op_string'])



    #Registering Z-transformed ReHo to standard space
    z_score_non_linear_mni = pe.MapNode(interface=fsl.ApplyWarp(),
                      name='z_score_non_linear_mni',
                      iterfield=['in_file',
                                 'premat'])

    spatial_smooth = pe.MapNode(interface=fsl.ImageMaths(),
                        name='spatial_smooth',
                        iterfield=['in_file'])



    reHo.connect(inputNode, 'rest_res_filt',
                    raw_reho_map, 'in_file')
    reHo.connect(inputNode, 'rest_mask',
                    raw_reho_map, 'mask_file')
    reHo.connect(inputNode, 'cluster_size',
                    raw_reho_map, 'cluster_size')

    reHo.connect(raw_reho_map, 'out_file',
                    mean, 'in_file')
    reHo.connect(inputNode, 'rest_mask',
                    mean, 'mask_file')
    reHo.connect(raw_reho_map, 'out_file',
                    standard_deviation, 'in_file')
    reHo.connect(inputNode, 'rest_mask',
                    standard_deviation, 'mask_file')
    reHo.connect(mean, 'out_stat',
                    op_string, 'mean')
    reHo.connect(standard_deviation, 'out_stat',
                    op_string, 'std_dev')
    reHo.connect(raw_reho_map, 'out_file',
                    z_score, 'in_file')
    reHo.connect(op_string, 'op_string',
                    z_score, 'op_string')
    reHo.connect(inputNode, 'rest_mask',
                    z_score, 'operand_files')
    reHo.connect(z_score, 'out_file',
                    z_score_non_linear_mni, 'in_file')
    reHo.connect(inputNode, 'standard',
                    z_score_non_linear_mni, 'ref_file')
    reHo.connect(inputNode, 'fieldcoeff_file',
                    z_score_non_linear_mni, 'field_file')
    reHo.connect(inputNode, 'premat',
                    z_score_non_linear_mni, 'premat')
    reHo.connect(z_score_non_linear_mni, 'out_file',
                 spatial_smooth, 'in_file')
    reHo.connect(inputnode_fwhm, ('fwhm', set_gauss_reho),
                 spatial_smooth, 'op_string')


    reHo.connect(z_score, 'out_file',
                 outputNode, 'z_score')

    reHo.connect(raw_reho_map, 'out_file',
                 outputNode, 'raw_reho_map')

    reHo.connect(z_score_non_linear_mni, 'out_file',
                 outputNode, 'z_score_non_linear_transformed')

    reHo.connect(spatial_smooth, 'out_file',
                 outputNode, 'z_score_non_linear_transformed_smooth')


    return reHo


