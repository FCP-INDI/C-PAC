# coding: utf-8
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from CPAC.reho.utils import *


def create_reho(wf_name):

    """
    Regional Homogeneity(ReHo) approach to fMRI data analysis

    This workflow computes the ReHo map, z-score on map

    Parameters
    ----------

    None

    Returns
    -------
    reHo : workflow
        Regional Homogeneity Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/reho/reho.py>`_

    Workflow Inputs: ::

        inputspec.rest_res_filt : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.rest_mask : string (existing nifti file)
            Input Whole Brain Mask of EPI 4D Volume

        inputspec.cluster_size : integer
            For a brain voxel the number of neighbouring brain voxels to use for KCC.
            Possible values are 27, 19, 7. Recommended value 27


    Workflow Outputs: ::

        outputspec.raw_reho_map : string (nifti file)

        outputspec.z_score : string (nifti file)


    ReHo Workflow Procedure:

    1. Generate ReHo map from the input EPI 4D volume, EPI mask and cluster_size
    2. Compute Z score of the ReHo map by subtracting mean and dividing by standard deviation

    .. exec::
        from CPAC.reho import create_reho
        wf = create_reho()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/reho.dot'
        )

    Workflow Graph:

    .. image:: ../../images/generated/reho.png
        :width: 500

    Detailed Workflow Graph:

    .. image:: ../../images/generated/reho_detailed.png
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
    >>> wf.inputs.inputspec.cluster_size = 27
    >>> wf.run()
    """

    reHo = pe.Workflow(wf_name=wf_name)
    inputNode = pe.Node(util.IdentityInterface(fields=['cluster_size',
                                                       'rest_res_filt',
                                                       'rest_mask']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['raw_reho_map']),
                         name='outputspec')

    reho_imports = ['import os', 'import sys', 'import nibabel as nb',
                    'import numpy as np',
                    'from CPAC.reho.utils import f_kendall']
    raw_reho_map = pe.Node(util.Function(input_names=['in_file', 'mask_file',
                                                      'cluster_size'],
                                         output_names=['out_file'],
                                         function=compute_reho,
                                         imports=reho_imports),
                           name='reho_map')

    reHo.connect(inputNode, 'rest_res_filt', raw_reho_map, 'in_file')
    reHo.connect(inputNode, 'rest_mask', raw_reho_map, 'mask_file')
    reHo.connect(inputNode, 'cluster_size', raw_reho_map, 'cluster_size')
    reHo.connect(raw_reho_map, 'out_file', outputNode, 'raw_reho_map')

    return reHo


def reho(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "ReHo",
     "config": ["regional_homogeneity"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-filtered_bold", "desc-cleaned_bold", "desc-brain_bold",
                 "desc-preproc_bold", "bold"],
                "space-bold_desc-brain_mask"],
     "outputs": ["reho"]}
    '''

    cluster_size = cfg.regional_homogeneity['cluster_size']

    # Check the cluster size is supported
    if cluster_size not in [7, 19, 27]:
        err_msg = 'Cluster size specified: %d, is not ' \
                  'supported. Change to 7, 19, or 27 and try ' \
                  'again' % cluster_size
        raise Exception(err_msg)

    reho = create_reho(f'reho_{pipe_num}')
    reho.inputs.inputspec.cluster_size = cluster_size


    node, out = strat_pool.get_data(["desc-filtered_bold",
                                     "desc-cleaned_bold",
                                     "desc-brain_bold",
                                     "desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, reho, 'inputspec.rest_res_filt')

    node, out_file = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out_file, reho, 'inputspec.rest_mask')

    outputs = {
        'reho': (reho, 'outputspec.raw_reho_map')
    }

    return (wf, outputs)
