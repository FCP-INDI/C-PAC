#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from CPAC.utils.function import Function

import os
import copy
import numpy as np
import nibabel as nb

from isfc import (
    isc,
    isc_significance,
    isc_permutation,
    isfc,
    isfc_significance,
    isfc_permutation,
)


def _permutations(perm):
    return range(perm)


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
            'D',
            'permutations',
            'collapse_subj',
            'two_sided',
            'random_state'
        ]),
        name='inputspec'
    )

    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'ISC', 'significance'
        ]),
        name='outputspec'
    )

    isc_node = pe.Node(Function(input_names=['D',
                                             'collapse_subj'],
                                output_names=['ISC'],
                                function=isc,
                                as_module=True),
                       name='ISC')

    permutations_node = pe.MapNode(Function(input_names=['permutation',
                                                         'D',
                                                         'collapse_subj',
                                                         'random_state'],
                                            output_names=['permutation',
                                                          'min_null',
                                                          'max_null'],
                                            function=isc_permutation,
                                            as_module=True),
                                   name='ISC_permutation', iterfield='permutation')

    significance_node = pe.Node(Function(input_names=['ISC',
                                                      'min_null',
                                                      'max_null',
                                                      'two_sided'],
                                         output_names=['p'],
                                         function=isc_significance,
                                         as_module=True),
                                name='ISC_p')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputspec, isc_node, [('D', 'D')]),
        (inputspec, isc_node, [('collapse_subj', 'collapse_subj')]),

        (isc_node, significance_node, [('ISC', 'ISC')]),

        (inputspec, permutations_node, [('D', 'D')]),
        (inputspec, permutations_node, [('collapse_subj', 'collapse_subj')]),
        (inputspec, permutations_node, [(('permutations', _permutations), 'permutation')]),
        (inputspec, permutations_node, [('random_state', 'random_state')]),

        (permutations_node, significance_node, [('min_null', 'min_null')]),
        (permutations_node, significance_node, [('max_null', 'max_null')]),
        (inputspec, significance_node, [('two_sided', 'two_sided')]),

        (isc_node, outputspec, [('ISC', 'ISC')]),
        (significance_node, outputspec, [('p', 'significance')]),
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
            'D',
            'permutations',
            'collapse_subj',
            'two_sided',
            'random_state'
        ]),
        name='inputspec'
    )

    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'ISFC', 'significance'
        ]),
        name='outputspec'
    )

    isfc_node = pe.Node(Function(input_names=['D',
                                             'collapse_subj'],
                                output_names=['ISFC'],
                                function=isfc,
                                as_module=True),
                       name='ISFC')

    permutations_node = pe.MapNode(Function(input_names=['permutation',
                                                         'D',
                                                         'collapse_subj',
                                                         'random_state'],
                                            output_names=['permutation',
                                                          'min_null',
                                                          'max_null'],
                                            function=isfc_permutation,
                                            as_module=True),
                                   name='ISFC_permutation', iterfield='permutation')

    significance_node = pe.Node(Function(input_names=['ISFC',
                                                      'min_null',
                                                      'max_null',
                                                      'two_sided'],
                                         output_names=['p'],
                                         function=isfc_significance,
                                         as_module=True),
                                name='ISFC_p')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputspec, isfc_node, [('D', 'D')]),
        (inputspec, isfc_node, [('collapse_subj', 'collapse_subj')]),

        (isfc_node, significance_node, [('ISFC', 'ISFC')]),

        (inputspec, permutations_node, [('D', 'D')]),
        (inputspec, permutations_node, [('collapse_subj', 'collapse_subj')]),
        (inputspec, permutations_node, [(('permutations', _permutations), 'permutation')]),
        (inputspec, permutations_node, [('random_state', 'random_state')]),

        (permutations_node, significance_node, [('min_null', 'min_null')]),
        (permutations_node, significance_node, [('max_null', 'max_null')]),
        (inputspec, significance_node, [('two_sided', 'two_sided')]),

        (isfc_node, outputspec, [('ISFC', 'ISFC')]),
        (significance_node, outputspec, [('p', 'significance')]),
    ])

    return wf