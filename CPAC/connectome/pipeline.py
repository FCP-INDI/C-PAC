#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

import os
import copy
import numpy as np
import nibabel as nb
from collections import OrderedDict as Odict

from CPAC.connectome.cross_validation import (
    CVInterface, JoinerInterface,
    Roi, Connectivity, Classifier
)

from CPAC.connectome.report import ReportInterface


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

    # inputspec = pe.Node(util.IdentityInterface(fields=[]),
    #                     name='inputspec')

    # outputspec = pe.Node(util.IdentityInterface(fields=[]),
    #                      name='outputspec')

    wf = pe.Workflow(name=name)

    from nilearn import datasets
    adhd_dataset = datasets.fetch_adhd(n_subjects=5)
    func_filenames = adhd_dataset.func
    cv = pe.Node(CVInterface(), name='cv')
    cv.inputs.X = func_filenames 
    cv.inputs.y = [1, 2, 2, 1, 2]
    cv.inputs.folds = int(configuration['folds'])

    cv_fields = ['X_train', 'y_train', 'X_valid', 'y_valid', 'fold', 'config']
    cv_fields_connection = list(zip(cv_fields, cv_fields))

    roi = pe.MapNode(
        interface=Roi(),
        name='roi',
        iterfield=cv_fields
    )
    roi.iterables = ('node_config', configuration['roi'])

    connectivity = pe.MapNode(
        interface=Connectivity(),
        name='connectivity',
        iterfield=cv_fields
    )
    connectivity.iterables = ('node_config', configuration['connectivity'])

    classifier = pe.MapNode(
        interface=Classifier(),
        name='classifier',
        iterfield=cv_fields
    )
    classifier.iterables = ('node_config', configuration['classifier'])

    joiner_classifier = pe.JoinNode(
        JoinerInterface(),
        name='joiner_classifier',
        joinsource=classifier,
        joinfield=cv_fields
    )

    joiner_connectivity = pe.JoinNode(
        JoinerInterface(),
        name='joiner_connectivity',
        joinsource=connectivity,
        joinfield=cv_fields
    )

    joiner_roi = pe.JoinNode(
        JoinerInterface(),
        name='joiner_roi',
        joinsource=roi,
        joinfield=cv_fields
    )

    report = pe.Node(interface=ReportInterface(),
                     name="report")
    
    wf.connect([
        (cv, roi, cv_fields_connection),
        (roi, connectivity, cv_fields_connection),
        (connectivity, classifier, cv_fields_connection),
        (classifier, joiner_classifier, cv_fields_connection),
        (joiner_classifier, joiner_connectivity, cv_fields_connection),
        (joiner_connectivity, joiner_roi, cv_fields_connection),
        (joiner_roi, report, cv_fields_connection),
        (cv, report, [('label_encoder', 'label_encoder')]),
    ])


    return wf
