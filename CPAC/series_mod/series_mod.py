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
from CPAC.series_mod.utils import compute_ROI_pcorr  
from CPAC.series_mod.utils import compute_MI  
from CPAC.series_mod.utils import compute_TE
#from CPAC.series_mod.utils import compute_ApEn


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
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
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



    corr_matNode = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['corr_mat'],
                     function=compute_ROI_corr),
                     name='corr_calc')


    ROI_corr.connect(inputNode, 'in_file',
                    corr_matNode, 'in_file')
    ROI_corr.connect(inputNode, 'mask_file',
                    corr_matNode, 'mask_file')  
                    
    ROI_corr.connect(corr_matNode, 'corr_mat',
                 outputNode, 'corr_mat')



    return ROI_corr


def create_ROI_pcorr():

    """
    Simple Network Partial Correlation Matrix calculation for ROIs in the mask file

    Parameters
    ----------

    None

    Returns
    -------
    pcorr : workflow
        Partial Correlation Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_

    Workflow Inputs: ::

        inputspec.in_file : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.mask_file : string (existing nifti file)
            Input Whole Brain Mask of EPI 4D Volume

    Workflow Outputs: ::

        outputspec.pcorr_matrix : ROI_number * ROI_number array


    Corr Workflow Procedure:

    1. Generate Partial Correlation Matrix from the input EPI 4D volume and EPI mask.


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
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
    wf.run()

    """



    ROI_pcorr = pe.Workflow(name='ROI_pcorr')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'pcorr_mat']),
                        name='outputspec')


    pcorr_matNode = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['pcorr_mat'],
                     function=compute_ROI_pcorr),
                     name='pcorr_calc')


    ROI_pcorr.connect(inputNode, 'in_file',
                    pcorr_matNode, 'in_file')
    ROI_pcorr.connect(inputNode, 'mask_file',
                    pcorr_matNode, 'mask_file')  
                    
    ROI_pcorr.connect(pcorr_matNode, 'pcorr_mat',
                 outputNode, 'pcorr_mat')



    return ROI_pcorr    
    
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
    >>> wf.run()
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    from CPAC import series_mod
    wf = series_mod.create_MI()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
    wf.run()

    """



    MI = pe.Workflow(name='MI_comp')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'MI_mat']),
                        name='outputspec')


    MI_mat = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['MI_mat'],
                     function=compute_MI),
                     name='MI_mat')


    MI.connect(inputNode, 'in_file',
                    MI_mat, 'in_file')
    MI.connect(inputNode, 'mask_file',
                    MI_mat, 'mask_file') 
                    
    MI.connect(MI_mat, 'MI_mat',
                 outputNode, 'MI_mat')

    return MI
    
    
def create_TE():

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

        outputspec.TE_matrix : [ROI_number * ROI_number] array

    Corr Workflow Procedure:

    1. Generate TE Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_TE()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.run()
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    wf = series_mod.create_TE()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
    wf.run()

    """



    TE = pe.Workflow(name='TE_comp')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'TE_mat']),
                        name='outputspec')


    TE_mat = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['MI_mat'],
                     function=compute_TE),
                     name='MI_mat')


    TE.connect(inputNode, 'in_file',
                    TE_mat, 'in_file')
    TE.connect(inputNode, 'mask_file',
                    TE_mat, 'mask_file') 
                    
    TE.connect(TE_mat, 'TE_mat',
                 outputNode, 'TE_mat')

    return TE
    
#def create_ApEn():
#
#    """
#    ApEn calculation for fMRI file
#
#    Parameters
#    ----------
#
#    None
#
#    Returns
#    -------
#    corr : workflow
#        Correlation Workflow
#
#    Notes
#    -----
#
#    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_
#
#    Workflow Inputs: ::
#
#        inputspec.in_file : string (existing nifti file)
#            Input EPI 4D Volume
#
#	inputspec.m_param : parameter m, time series size
#
#        inputspec.r_param : parameter r - factor?¿
#
#    Workflow Outputs: ::
#
#        outputspec.result_vector : ApEn value
#        
#    References
#    ---------- 
#
#    Examples
#    --------
#    >>> from CPAC import series_mod
#    >>> wf = series_mod.create_ApEn()
#    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
#    >>> wf.inputs.inputspec.m_param = 30
#    >>> wf.inputs.inputspec.m_param = 4
#    >>> wf.run()
#
#    """
#
#
#
#    ApEn = pe.Workflow(name='ApEn')
#    inputNode = pe.Node(util.IdentityInterface(fields=[
#                                                'in_file',
#                                                'm_param',
#						'r_param'
#                                                ]),
#                        name='inputspec')
#
#
#    outputNode = pe.Node(util.IdentityInterface(fields=[
#                                                    'result_vector']),
#                        name='outputspec')
#
#
#
#    ApEn_calc_Node = pe.Node(util.Function(input_names=['in_file', 'm_param','r_param'],
#                                   output_names=['result_vector'],
#                     function=compute_ApEn),
#                     name='ApEn_calc')
#
#
#    ApEn.connect(inputNode, 'in_file',
#                    ApEn_calc_Node, 'in_file')
#    ApEn.connect(inputNode, 'm_param',
#                    ApEn_calc_Node, 'm_param')  
#    ApEn.connect(inputNode, 'r_param',
#                    ApEn_calc_Node, 'r_param')  
#                    
#    ApEn.connect(ApEn_calc_Node, 'result_vector',
#                 outputNode, 'result_vector')
#
#
#
#    return ApEn    