import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec
from nipype.interfaces import fsl

from nilearn import input_data, masking, image, datasets
from nilearn.image import resample_to_img, concat_imgs
from nilearn.input_data import NiftiMasker, NiftiLabelsMasker

from CPAC.utils.function import Function

import os
import copy
import numpy as np
import nibabel as nb


def select(input_list):
    out_file = input_list[0]
    return out_file


def create_randomise(name='randomise', working_dir=None, crash_dir=None):
    """
    
    Parameters
    ----------
        
    Returns
    -------
    workflow : nipype.pipeline.engine.Workflow
        Randomise workflow.
        
    Notes
    -----
    
    Workflow Inputs::
    
        
    Workflow Outputs::

    
    References
    ----------
    
    """

    if not working_dir:
        working_dir = os.path.join(os.getcwd(), 'Randomise_work_dir')
    if not crash_dir:
        crash_dir = os.path.join(os.getcwd(), 'Randomise_crash_dir')

    wf = pe.Workflow(name=name)
    wf.base_dir = working_dir
    wf.config['execution'] = {'hash_method': 'timestamp',
                              'crashdump_dir': os.path.abspath(crash_dir)}

    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'subjects',
            'design_matrix_file',
            'constrast_file',
            'f_constrast_file',
            'permutations',
        ]),
        name='inputspec'
    )

    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'index_file'
            'threshold_file'
            'localmax_txt_file'
            'localmax_vol_file'
            'max_file'
            'mean_file'
            'pval_file'
            'size_file'
        ]),
        name='outputspec'
    )

    merge = pe.Node(interface=fsl.Merge(), name='fsl_merge')
    merge.inputs.dimension = 't'
    merge.inputs.merged_file = "randomise_merged.nii.gz"

    wf.connect(inputspec, 'subjects', merge, 'in_files')

    mask = pe.Node(interface=fsl.maths.MathsCommand(), name='fsl_maths')
    mask.inputs.args = '-abs -Tmin -bin'
    mask.inputs.out_file = "randomise_mask.nii.gz"

    wf.connect(merge, 'merged_file', mask, 'in_file')

    randomise = pe.Node(interface=fsl.Randomise(), name='randomise')
    wf.connect(mask, 'out_file', randomise, 'mask')
    randomise.inputs.base_name = "randomise"
    randomise.inputs.demean = True
    randomise.inputs.tfce = True
    wf.connect([
        (merge, randomise, [('merged_file', 'in_file')]),
        (inputspec, randomise, [
            ('design_matrix_file', 'design_mat'),
            ('constrast_file', 'tcon'),
            ('f_constrast_file', 'fcon'),
            ('permutations', 'num_perm'),
        ]),
    ])

    select_t_corrected = pe.Node(Function(input_names=["input_list"],
                                          output_names=['out_file'],
                                          function=select),
                                 name='select_t_cor')

    wf.connect(randomise, "t_corrected_p_files",
               select_t_corrected, "input_list")

    thresh = pe.Node(interface=fsl.Threshold(),
                     name='fsl_threshold_contrast')
    thresh.inputs.thresh = 0.95
    thresh.inputs.out_file = 'rando_pipe_thresh_tstat.nii.gz'

    wf.connect(select_t_corrected, "out_file", thresh, "in_file")

    thresh_bin = pe.Node(interface=fsl.maths.MathsCommand(),
                         name='fsl_threshold_bin_contrast')
    thresh_bin.inputs.args = '-bin'
    wf.connect(thresh, "out_file", thresh_bin, "in_file")

    select_t_stat = pe.Node(Function(input_names=["input_list"],
                                     output_names=['out_file'],
                                     function=select),
                            name='select_item_t_stat')

    wf.connect(randomise, "tstat_files", select_t_stat, "input_list")

    apply_mask = pe.Node(interface=fsl.ApplyMask(),
                         name='fsl_applymask_contrast')
    wf.connect(select_t_stat, 'out_file', apply_mask, 'in_file')
    wf.connect(thresh_bin, 'out_file', apply_mask, 'mask_file')

    cluster = pe.Node(interface=fsl.Cluster(),
                      name='cluster_contrast')
    cluster.inputs.threshold = 0.0001
    cluster.inputs.out_index_file = "cluster_index_contrast"
    cluster.inputs.out_localmax_txt_file = "lmax_contrast.txt"
    cluster.inputs.out_size_file = "cluster_size_contrast"
    cluster.inputs.out_threshold_file = "randomise_out_contrast"
    cluster.inputs.terminal_output = 'file'

    wf.connect(apply_mask, 'out_file', cluster, 'in_file')

    wf.connect([
        (cluster, outputspec, [
            ('index_file', 'index_file')
            ('threshold_file', 'threshold_file')
            ('localmax_txt_file', 'localmax_txt_file')
            ('localmax_vol_file', 'localmax_vol_file')
            ('max_file', 'max_file')
            ('mean_file', 'mean_file')
            ('pval_file', 'pval_file')
            ('size_file', 'size_file')
        ])
    ])

    return wf
