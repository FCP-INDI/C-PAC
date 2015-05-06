# -*- coding: utf-8 -*-

# coding: utf-8
#import os
#import sys
#import re
#import commands
import nipype.pipeline.engine as pe
#import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from CPAC.series_mod.utils import * #this import is wrong... ?¿


def create_corr():

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

        inputspec.rest_res_filt : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.rest_mask : string (existing nifti file)
            Input Whole Brain Mask of EPI 4D Volume

    Workflow Outputs: ::

        outputspec.corr_matrix : string (nifti file)


    Corr Workflow Procedure:

    1. Generate Correlation Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_corr()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.inputs.inputspec.TR = 3
    >>> wf.run()
    """




    corr = pe.Workflow(name='corr')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask',
                                                'TR'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'corr_mat']),
                        name='outputspec')


    corr_mat = pe.Node(util.Function(input_names=['in_file', 'mask','TR'],
                                   output_names=['out_file'],
                     function=compute_corr),
                     name='corr_mat')


    corr.connect(inputNode, 'in_file',
                    corr_mat, 'in_file')
    corr.connect(inputNode, 'mask',
                    corr_mat, 'mask')
    corr.connect(inputNode, 'TR',
                    corr_mat, 'TR')                   
    corr.connect(corr_mat, 'out_file',
                 outputNode, 'corr_mat')



    return corr

def create_te():

    """
    Simple Pairwise Transfer Entropy Matrix calculation for ROIs in the mask file

    Parameters
    ----------

    None

    Returns
    -------
    te : workflow
        Transfer Entropy Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_

    Workflow Inputs: ::

        inputspec.rest_res_filt : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.rest_mask : string (existing nifti file)
            Input Whole Brain Mask of EPI 4D Volume

    Workflow Outputs: ::

        outputspec.te_matrix : string (nifti file)


    Corr Workflow Procedure:

    1. Generate Pairwise Transfer Entropy Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_te()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.inputs.inputspec.TR = 3
    >>> wf.run()
    """




    te = pe.Workflow(name='te')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask',
                                                'TR'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'te_mat']),
                        name='outputspec')


    te_mat = pe.Node(util.Function(input_names=['in_file', 'mask','TR'],
                                   output_names=['out_file'],
                     function=compute_te),
                     name='te_mat')


    te.connect(inputNode, 'in_file',
                    te_mat, 'in_file')
    te.connect(inputNode, 'mask',
                    te_mat, 'mask')
    te.connect(inputNode, 'TR',
                    te_mat, 'TR')                
    te.connect(corr_mat, 'out_file',
                 outputNode, 'te_mat')
                 
    return te
