#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from nilearn import input_data, masking, image, datasets
from nilearn.image import resample_to_img, concat_imgs
from nilearn.input_data import NiftiMasker, NiftiLabelsMasker

from CPAC.utils.function import Function

import os
import copy
import numpy as np
import nibabel as nb

from CPAC.isc.isc import (
    isc,
    isc_significance,
    isc_permutation,
)

from CPAC.isc.isfc import (
    isfc,
    isfc_significance,
    isfc_permutation,
)


def _permutations(perm):
    return range(perm)


def load_data(subjects):
    subject_ids = list(subjects.keys())
    subject_files = list(subjects.values())

    if subject_files[0].endswith('.csv'):
        data = np.array([
            np.genfromtxt(img).T
            for img in subject_files
        ])
        voxel_masker = None

    else:
        images = [nb.load(img) for img in subject_files]
        voxel_masker = NiftiMasker()
        voxel_masker.fit(images)

        data = np.array([
            voxel_masker.transform(img).T
            for img in images
        ])

    # Reshape to voxel x time x subject
    data = np.moveaxis(data, 0, -1).copy(order='C')

    data_file = os.path.abspath('./data.npy')
    np.save(data_file, data)
    return data_file, voxel_masker


def save_data_isc(ISC, p, voxel_masker=None):

    ISC = np.load(ISC)
    p = np.load(p)

    if voxel_masker:
        ISC_image = voxel_masker.inverse_transform(ISC.T)
        p_image = voxel_masker.inverse_transform(p.T)
        final = concat_imgs([ISC_image, p_image])
        final_file = os.path.abspath('./data.nii.gz')
        nb.save(final, final_file)
    else:
        final = np.array([
            ISC, p
        ])
        final_file = os.path.abspath('./data.csv')
        np.savetxt(final_file, final, delimiter=',', fmt='%1.10f')

    return final_file


def save_data_isfc(ISFC, p):
    final = np.array([
        np.load(ISFC), np.load(p)
    ])
    final_file = os.path.abspath('./result.npy')
    np.save(final_file, final)
    return final_file


def node_isc(D, std=None, collapse_subj=True):
    D = np.load(D)

    ISC, ISC_mask = isc(D, std, collapse_subj=True)
    
    f = os.path.abspath('./isc.npy')
    np.save(f, ISC)

    f_mask = os.path.abspath('./isc_mask.npy')
    np.save(f_mask, ISC_mask)
    
    return f, f_mask


def node_isc_significance(ISC, min_null, max_null, two_sided=False):
    ISC = np.load(ISC)
    p = isc_significance(ISC, min_null, max_null, two_sided)
    f = os.path.abspath('./isc-p.npy')
    np.save(f, p)
    return f


def node_isc_permutation(permutation, D, masked, collapse_subj=True, random_state=0):
    D = np.load(D)
    masked = np.load(masked)
    permutation, min_null, max_null = isc_permutation(permutation,
                                                      D,
                                                      masked,
                                                      collapse_subj,
                                                      random_state)
    return permutation, min_null, max_null


def node_isfc(D, std=None, collapse_subj=True):
    D = np.load(D)

    ISFC, ISFC_mask = isfc(D, std, collapse_subj)

    f = os.path.abspath('./isfc.npy')
    np.save(f, ISFC)

    f_mask = os.path.abspath('./isfc_mask.npy')
    np.save(f_mask, ISFC_mask)
    
    return f, f_mask


def node_isfc_significance(ISFC, min_null, max_null, two_sided=False):
    ISFC = np.load(ISFC)
    p = isfc_significance(ISFC, min_null, max_null, two_sided)
    f = os.path.abspath('./isfc-p.npy')
    np.save(f, p)
    return f


def node_isfc_permutation(permutation, D, masked, collapse_subj=True, random_state=0):
    D = np.load(D)
    masked = np.load(masked)
    permutation, min_null, max_null = isfc_permutation(permutation,
                                                       D,
                                                       masked,
                                                       collapse_subj,
                                                       random_state)
    return permutation, min_null, max_null


def create_isc(name='isc'):
    """
    Inter-Subject Correlation
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
        
    Returns
    -------
    workflow : nipype.pipeline.engine.Workflow
        ISFC workflow.
        
    Notes
    -----
    
    Workflow Inputs::
    
        
    Workflow Outputs::

    
    References
    ----------
    .. [1] Simony, E., Honey, C. J., Chen, J., Lositsky, O., Yeshurun,
           Y., Wiesel, A., & Hasson, U. (2016). Dynamic reconfiguration of the
           default mode network during narrative comprehension.
           Nature Communications, 7(May 2015), 1-13.
           https://doi.org/10.1038/ncomms12141
    
    """

    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'subjects',
            'permutations',
            'collapse_subj',
            'std',
            'two_sided',
            'random_state'
        ]),
        name='inputspec'
    )

    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'result'
        ]),
        name='outputspec'
    )

    data_node = pe.Node(Function(input_names=['subjects'],
                                 output_names=['D', 'voxel_masker'],
                                 function=load_data,
                                 as_module=True),
                        name='data')

    save_node = pe.Node(Function(input_names=['ISC', 'p', 'voxel_masker'],
                                 output_names=['result'],
                                 function=save_data_isc,
                                 as_module=True),
                        name='save')

    isc_node = pe.Node(Function(input_names=['D',
                                             'std',
                                             'collapse_subj'],
                                output_names=['ISC', 'masked'],
                                function=node_isc,
                                as_module=True),
                       name='ISC')

    permutations_node = pe.MapNode(Function(input_names=['permutation',
                                                         'D',
                                                         'masked',
                                                         'collapse_subj',
                                                         'random_state'],
                                            output_names=['permutation',
                                                          'min_null',
                                                          'max_null'],
                                            function=node_isc_permutation,
                                            as_module=True),
                                   name='ISC_permutation', iterfield='permutation')

    significance_node = pe.Node(Function(input_names=['ISC',
                                                      'min_null',
                                                      'max_null',
                                                      'two_sided'],
                                         output_names=['p'],
                                         function=node_isc_significance,
                                         as_module=True),
                                name='ISC_p')

    wf = pe.Workflow(name=name)

    wf.connect([
        (inputspec, data_node, [('subjects', 'subjects')]),
        (inputspec, isc_node, [('collapse_subj', 'collapse_subj')]),
        (inputspec, isc_node, [('std', 'std')]),
        (data_node, isc_node, [('D', 'D')]),

        (isc_node, significance_node, [('ISC', 'ISC')]),

        (data_node, permutations_node, [('D', 'D')]),
        (isc_node, permutations_node, [('masked', 'masked')]),
        (inputspec, permutations_node, [('collapse_subj', 'collapse_subj')]),
        (inputspec, permutations_node, [(('permutations', _permutations), 'permutation')]),
        (inputspec, permutations_node, [('random_state', 'random_state')]),

        (permutations_node, significance_node, [('min_null', 'min_null')]),
        (permutations_node, significance_node, [('max_null', 'max_null')]),
        (inputspec, significance_node, [('two_sided', 'two_sided')]),

        (isc_node, save_node, [('ISC', 'ISC')]),
        (significance_node, save_node, [('p', 'p')]),
        (data_node, save_node, [('voxel_masker', 'voxel_masker')]),

        (save_node, outputspec, [('result', 'result')]),
    ])

    return wf


def create_isfc(name='isfc'):
    """
    Inter-Subject Functional Correlation
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
        
    Returns
    -------
    workflow : nipype.pipeline.engine.Workflow
        ISFC workflow.
        
    Notes
    -----
    
    Workflow Inputs::
    
        
    Workflow Outputs::

    
    References
    ----------
    .. [1] Simony, E., Honey, C. J., Chen, J., Lositsky, O., Yeshurun,
           Y., Wiesel, A., & Hasson, U. (2016). Dynamic reconfiguration of the
           default mode network during narrative comprehension.
           Nature Communications, 7(May 2015), 1-13.
           https://doi.org/10.1038/ncomms12141
    
    """

    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'subjects',
            'permutations',
            'collapse_subj',
            'std',
            'two_sided',
            'random_state'
        ]),
        name='inputspec'
    )

    data_node = pe.Node(Function(input_names=['subjects'],
                                 output_names=['D', 'voxel_masker'],
                                 function=load_data,
                                 as_module=True),
                        name='data')

    save_node = pe.Node(Function(input_names=['ISFC', 'p'],
                                 output_names=['result'],
                                 function=save_data_isfc,
                                 as_module=True),
                        name='save')

    outputspec = pe.Node(
        util.IdentityInterface(fields=['result']),
        name='outputspec'
    )

    isfc_node = pe.Node(Function(input_names=['D',
                                             'std',
                                             'collapse_subj'],
                                output_names=['ISFC', 'masked'],
                                function=node_isfc,
                                as_module=True),
                       name='ISFC')

    permutations_node = pe.MapNode(Function(input_names=['permutation',
                                                         'D',
                                                         'masked',
                                                         'collapse_subj',
                                                         'random_state'],
                                            output_names=['permutation',
                                                          'min_null',
                                                          'max_null'],
                                            function=node_isfc_permutation,
                                            as_module=True),
                                   name='ISFC_permutation', iterfield='permutation')

    significance_node = pe.Node(Function(input_names=['ISFC',
                                                      'min_null',
                                                      'max_null',
                                                      'two_sided'],
                                         output_names=['p'],
                                         function=node_isfc_significance,
                                         as_module=True),
                                name='ISFC_p')

    wf = pe.Workflow(name=name)

    wf.connect([
        (inputspec, data_node, [('subjects', 'subjects')]),
        (inputspec, isfc_node, [('collapse_subj', 'collapse_subj')]),
        (inputspec, isfc_node, [('std', 'std')]),
        (data_node, isfc_node, [('D', 'D')]),

        (isfc_node, significance_node, [('ISFC', 'ISFC')]),

        (data_node, permutations_node, [('D', 'D')]),
        (isfc_node, permutations_node, [('masked', 'masked')]),
        (inputspec, permutations_node, [('collapse_subj', 'collapse_subj')]),
        (inputspec, permutations_node, [(('permutations', _permutations), 'permutation')]),
        (inputspec, permutations_node, [('random_state', 'random_state')]),

        (permutations_node, significance_node, [('min_null', 'min_null')]),
        (permutations_node, significance_node, [('max_null', 'max_null')]),
        (inputspec, significance_node, [('two_sided', 'two_sided')]),

        (isfc_node, save_node, [('ISFC', 'ISFC')]),
        (significance_node, save_node, [('p', 'p')]),

        (save_node, outputspec, [('result', 'result')]),
    ])

    return wf
