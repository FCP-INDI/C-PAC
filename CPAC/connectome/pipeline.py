#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

import os
import numpy as np
import nibabel as nb
from nilearn.connectome import ConnectivityMeasure


def create_connectome_study(configuration, name='connectivity'):
    """
    Functional Connectivity
    
    This workflow computes different kinds of functional connectivity matrices.
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
        
    Returns
    -------
    workflow : nipype.pipeline.engine.Workflow
        Connectivity workflow.
        
    Notes
    -----
    
    Workflow Inputs::
    
        inputspec.functional : string
            Path of the subject's functional nifti file.
        inputspec.metric : string
            Metric to compute connectivity: correlation,
                                            partial correlation (or partial),
                                            tangent
        
    Workflow Outputs::

        outputspec.connectivity : string
           Path to numpy file with correlation matrix
    
    References
    ----------
    .. [1] G. Varoquaux et al. “Detection of brain functional-connectivity difference in post-stroke patients using group-level covariance modeling, MICCAI 2010.
    
    """

    inputspec = pe.Node(util.IdentityInterface(fields=[]),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=[]),
                         name='outputspec')

    wf = pe.Workflow(name=name)

    {
    'classifier': [{'parameters': None, 'type': 'knn'},
                    {'parameters': None, 'type': 'random_forest'},
                    {'parameters': None, 'type': 'naive_bayes'},
                    {'parameters': None, 'type': 'ridge'},
                    {'feature_extraction': [{'parameters': {'percentile': 10},
                                             'type': 'anova'}],
                     'parameters': {'penalty': 'l1'},
                     'type': 'svc'},
                    {'feature_extraction': [{'parameters': {'percentile': 10},
                                             'type': 'anova'}],
                     'parameters': {'penalty': 'l2'},
                     'type': 'svc'},
                    {'parameters': {'penalty': 'l1'}, 'type': 'svc'},
                    {'parameters': {'penalty': 'l2'}, 'type': 'svc'},
                    {'parameters': {'penalty': 'l1'},
                     'type': 'logistic'},
                    {'parameters': {'penalty': 'l2'}, 'type': 'logistic'}],

    'connectivity': ['partial', 'tangent', 'correlation'],

    'roi': ['aal', 'harvard_oxford', 'power', 'basc', 'ward', 'k_means', 'group_ica', 'dict_learning']}


    conn = pe.Node(FunctionalConnectivity(),
                   name='connectivity')
    
    conn.inputs.vectorize = True
    conn.iterables = ('metric', [
        'correlation',
        'partial',
        'tangent',
    ])

    wf.connect(inputspec, 'functional',
               conn, 'functional')

    wf.connect(conn, 'connectivity',
               outputspec, 'connectivity')

    return wf
