# -*- coding: utf-8 -*-

# coding: utf-8
#import os
#import sys
#import re
#import commands
import nipype.pipeline.engine as pe
#import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from CPAC.series_mod.utils import compute_ROI_corr
from CPAC.series_mod.utils import compute_MI  

def create_ROI_corr():

    """
    Simple Network Correlation Matrix calculation for ROIs in the mask file

    Parameters
    ----------

    None

    Returns
    -------
    corr : workflow
        Correlation Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_

    Workflow Inputs: ::

        inputspec.in_file : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.mask_file : string (existing nifti file)
            Input Whole Brain Mask of EPI 4D Volume

    Workflow Outputs: ::

        outputspec.corr_matrix : ROI_number * ROI_number array


    Corr Workflow Procedure:

    1. Generate Correlation Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_ROI_corr()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.run()
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    from CPAC import series_mod
    wf = series_mod.create_ROI_corr()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.run()

    """



    ROI_corr = pe.Workflow(name='ROI_corr')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'corr_mat']),
                        name='outputspec')


    corr_mat = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['corr_mat'],
                     function=compute_ROI_corr),
                     name='corr_mat')


    ROI_corr.connect(inputNode, 'in_file',
                    corr_mat, 'in_file')
    ROI_corr.connect(inputNode, 'mask_file',
                    corr_mat, 'mask_file')  
                    
    ROI_corr.connect(corr_mat, 'corr_mat',
                 outputNode, 'corr_mat')



    return ROI_corr
    
    
def create_MI():

    """
    Simple Network Correlation Matrix calculation for ROIs in the mask file

    Parameters
    ----------

    None

    Returns
    -------
    corr : workflow
        Correlation Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_

    Workflow Inputs: ::

        inputspec.in_file : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.bins : integer
            How many states do you want to be the data transformed to

    Workflow Outputs: ::

        outputspec.MI_matrix : [ROI_number * ROI_number] array

    Corr Workflow Procedure:

    1. Generate MI Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_MI()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.inputs.inputspec.bins = 4
    >>> wf.run()
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    from CPAC import series_mod
    wf = series_mod.create_MI()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.inputs.inputspec.bins = 10
    wf.run()

    """



    MI = pe.Workflow(name='MI_comp')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file',
                                                'bins'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'MI_mat']),
                        name='outputspec')


    MI_mat = pe.Node(util.Function(input_names=['in_file', 'mask_file','bins'],
                                   output_names=['MI_mat'],
                     function=compute_MI),
                     name='MI_mat')


    MI.connect(inputNode, 'in_file',
                    MI_mat, 'in_file')
    MI.connect(inputNode, 'mask_file',
                    MI_mat, 'mask_file') 
    MI.connect(inputNode, 'bins',
                    MI_mat, 'bins')
                    
    MI.connect(MI_mat, 'MI_mat',
                 outputNode, 'MI_mat')



    return MI