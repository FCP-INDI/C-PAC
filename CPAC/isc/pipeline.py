#!/usr/bin/env python
# -*- coding: utf-8 -*-

from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from nilearn import input_data, masking, image, datasets
from nilearn.image import resample_to_img, concat_imgs
from nilearn.input_data import NiftiMasker, NiftiLabelsMasker

from CPAC.utils.interfaces.function import Function

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
    subject_files = list(subjects[i] for i in subject_ids)

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
    return subject_ids, data_file, voxel_masker


def save_data_isc(subject_ids, ISC, p, out_dir, collapse_subj=True, voxel_masker=None):

    ISC = np.load(ISC)
    p = np.load(p)

    subject_ids_file = os.path.abspath('./subject_ids.txt')
    np.savetxt(subject_ids_file, np.array(subject_ids), fmt="%s")

    if voxel_masker:
        corr_file, p_file = os.path.abspath('./correlations.nii.gz'), os.path.abspath('./significance.nii.gz')
        corr_out, p_out = os.path.join(out_dir, 'correlations.nii.gz'), os.path.join(out_dir, 'significance.nii.gz')

        ISC_image = voxel_masker.inverse_transform(ISC.T)
        p_image = voxel_masker.inverse_transform(p.T)

        nb.save(ISC_image, corr_file)
        nb.save(p_image, p_file)

        os.makedirs(out_dir)
        nb.save(ISC_image, corr_out)
        nb.save(p_image, p_out)

    else:
        corr_file, p_file = os.path.abspath('./correlations.csv'), os.path.abspath('./significance.csv')
        corr_out, p_out = os.path.join(out_dir, 'correlations.csv'), os.path.join(out_dir, 'significance.csv')

        if ISC.ndim == 1:
            ISC = np.expand_dims(ISC, axis=0)

        if p.ndim == 1:
            p = np.expand_dims(p, axis=0)

        np.savetxt(corr_file, ISC, delimiter=',', fmt='%1.10f')
        np.savetxt(p_file, p, delimiter=',', fmt='%1.10f')

        os.makedirs(out_dir)
        np.savetxt(corr_out, ISC, delimiter=',', fmt='%1.10f')
        np.savetxt(p_out, p, delimiter=',', fmt='%1.10f')

    return subject_ids_file, corr_file, p_file


def save_data_isfc(subject_ids, ISFC, p, out_dir, collapse_subj=True):

    subject_ids_file = os.path.abspath('./subject_ids.txt')
    np.savetxt(subject_ids_file, np.array(subject_ids), fmt="%s")

    os.makedirs(out_dir)

    corr_file = os.path.abspath('./correlations.npy')
    corr_out = os.path.join(out_dir, 'correlations.npy')

    corr = np.load(ISFC)
    np.save(corr_file, corr if collapse_subj else np.moveaxis(corr, -1, 0))
    np.save(corr_out, corr if collapse_subj else np.moveaxis(corr, -1, 0))    

    p_file = os.path.abspath('./significance.npy')
    p_out = os.path.join(out_dir, 'significance.npy')

    p = np.load(p)
    np.save(p_file, p if collapse_subj else np.moveaxis(p, -1, 0))
    np.save(p_out, p if collapse_subj else np.moveaxis(p, -1, 0))

    return subject_ids_file, corr_file, p_file


def node_isc(D, std=None, collapse_subj=True):
    D = np.load(D)

    ISC, ISC_mask = isc(D, std, collapse_subj)
    
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


def create_isc(name='isc', output_dir=None, working_dir=None, crash_dir=None):
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

    if not output_dir:
        output_dir = os.path.join(os.getcwd(), 'ISC_output_dir')
    if not working_dir:
        working_dir = os.path.join(os.getcwd(), 'ISC_work_dir')
    if not crash_dir:
        crash_dir = os.path.join(os.getcwd(), 'ISC_crash_dir')

    wf = pe.Workflow(name=name)
    wf.base_dir = working_dir
    wf.config['execution'] = {'hash_method': 'timestamp',
                              'crashdump_dir': os.path.abspath(crash_dir)}

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
            'correlations',
            'significance'
        ]),
        name='outputspec'
    )

    data_node = pe.Node(Function(input_names=['subjects'],
                                 output_names=['subject_ids', 'D', 'voxel_masker'],
                                 function=load_data,
                                 as_module=True),
                        name='data')

    save_node = pe.Node(Function(input_names=['subject_ids', 'ISC', 'p', 'out_dir', 'collapse_subj', 'voxel_masker'],
                                 output_names=['subject_ids', 'correlations', 'significance'],
                                 function=save_data_isc,
                                 as_module=True),
                        name='save')
    save_node.inputs.out_dir = output_dir

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
        (inputspec, save_node, [('collapse_subj', 'collapse_subj')]),
        (significance_node, save_node, [('p', 'p')]),
        (data_node, save_node, [('subject_ids', 'subject_ids')]),
        (data_node, save_node, [('voxel_masker', 'voxel_masker')]),

        (save_node, outputspec, [('subject_ids', 'subject_ids')]),
        (save_node, outputspec, [('correlations', 'correlations')]),
        (save_node, outputspec, [('significance', 'significance')]),
    ])

    return wf


def create_isfc(name='isfc', output_dir=None, working_dir=None,
                crash_dir=None):
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


    if not output_dir:
        output_dir = os.path.join(os.getcwd(), 'ISC_output_dir')
    if not working_dir:
        working_dir = os.path.join(os.getcwd(), 'ISC_work_dir')
    if not crash_dir:
        crash_dir = os.path.join(os.getcwd(), 'ISC_crash_dir')

    wf = pe.Workflow(name=name)
    wf.base_dir = working_dir
    wf.config['execution'] = {'hash_method': 'timestamp',
                              'crashdump_dir': os.path.abspath(crash_dir)}

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
                                 output_names=['subject_ids', 'D', 'voxel_masker'],
                                 function=load_data,
                                 as_module=True),
                        name='data')

    save_node = pe.Node(Function(input_names=['subject_ids', 'ISFC', 'p', 'out_dir', 'collapse_subj'],
                                 output_names=['subject_ids', 'correlations', 'significance'],
                                 function=save_data_isfc,
                                 as_module=True),
                        name='save')
    save_node.inputs.out_dir = output_dir

    outputspec = pe.Node(
        util.IdentityInterface(fields=['correlations', 'significance']),
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

        (data_node, save_node, [('subject_ids', 'subject_ids')]),
        (inputspec, save_node, [('collapse_subj', 'collapse_subj')]),
        (isfc_node, save_node, [('ISFC', 'ISFC')]),
        (significance_node, save_node, [('p', 'p')]),

        (save_node, outputspec, [('subject_ids', 'subject_ids')]),
        (save_node, outputspec, [('correlations', 'correlations')]),
        (save_node, outputspec, [('significance', 'significance')]),
    ])

    return wf
