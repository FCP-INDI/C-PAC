#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

import os
import numpy as np
import nibabel as nb

from CPAC.connectome.cross_validation import CVInterface
from CPAC.connectome.classifiers import (
    SVCClassifierInterface
)
from CPAC.connectome.connectivity import (
    CorrelationConnectivityInterface,
    PartialCorrelationConnectivityInterface,
    TangentConnectivityInterface
)
from CPAC.connectome.rois import (
    DictLearningRoiInterface
)


def generate_layer(wf, ifaces_ref, configs, leafs):

    nodes = []

    for leaf, leaf_train, leaf_valid in leafs:

        for config in configs:

            iface = ifaces_ref[config['type']]()
            params = config['parameters'] if 'parameters' in config else {}

            name = "%s_%s" % (
                config['type'],
                '_'.join(
                    '%s_%s' % (k, params[k])
                    for k in sorted(params.keys())
                ))
            name = name.strip('_')
            if leaf != '':
                name = '%s_%s' % (leaf, name)

            node_train = pe.MapNode(interface=iface, name='%s_train' % name, iterfield=['X', 'y', 'fold'])
            node_valid = pe.MapNode(interface=iface, name='%s_valid' % name, iterfield=['X', 'y', 'fold', 'model'])
            # TODO setup parameters

            wf.connect(node_train, 'model', node_valid, 'model')

            wf.connect([
                (leaf_train, node_train, [
                    ('X', 'X'),
                    ('y', 'y'),
                    ('fold', 'fold')
                ]),
                (leaf_valid, node_valid, [
                    ('X', 'X'),
                    ('y', 'y'),
                    ('fold', 'fold')
                ]),
            ])

            nodes += [(name, node_train, node_valid)]

    return nodes


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

    cv = pe.Node(CVInterface(), name='cv')
    cv.inputs.folds = 2
    from nilearn import datasets
    adhd_dataset = datasets.fetch_adhd(n_subjects=5)
    func_filenames = adhd_dataset.func
    cv.inputs.X = func_filenames
    cv.inputs.y = [1, 2, 2, 1, 2]

    input_train = pe.MapNode(util.IdentityInterface(fields=['X', 'y', 'fold']), name='input_train', iterfield=['X', 'y', 'fold'])
    input_valid = pe.MapNode(util.IdentityInterface(fields=['X', 'y', 'fold']), name='input_valid', iterfield=['X', 'y', 'fold'])

    wf.connect([
        (cv, input_train, [
            ('train_X', 'X'),
            ('train_y', 'y'),
            ('fold', 'fold')
        ]),
        (cv, input_valid, [
            ('train_X', 'X'),
            ('train_y', 'y'),
            ('fold', 'fold')
        ]),
    ])

    input_leafs = [('', input_train, input_valid)]

    roi_ifaces = {
        # 'aal': DictLearningRoiInterface,
        # 'harvard_oxford': DictLearningRoiInterface,
        # 'power': DictLearningRoiInterface,
        # 'basc': DictLearningRoiInterface,
        # 'ward': DictLearningRoiInterface,
        # 'k_means': DictLearningRoiInterface,
        # 'group_ica': DictLearningRoiInterface,
        'dict_learning': DictLearningRoiInterface
    }

    connectivity_ifaces = {
        'correlation': CorrelationConnectivityInterface,
        'partial correlation': PartialCorrelationConnectivityInterface,
        'partial': PartialCorrelationConnectivityInterface,
        'tangent': TangentConnectivityInterface,
    }

    classifier_ifaces = {
        'svc': SVCClassifierInterface
    }

    roi_nodes          = generate_layer(wf, roi_ifaces,          configuration['roi'],          input_leafs)
    connectivity_nodes = generate_layer(wf, connectivity_ifaces, configuration['connectivity'], roi_nodes)
    classifier_nodes   = generate_layer(wf, classifier_ifaces,   configuration['classifier'],   connectivity_nodes)

    return wf
