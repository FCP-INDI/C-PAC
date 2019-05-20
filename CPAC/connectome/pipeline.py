#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from CPAC.utils.interfaces.function import Function

import os
import numpy as np
from nilearn.connectome import ConnectivityMeasure


def compute_correlation(time_series, method):

    data = np.genfromtxt(time_series).T

    if method == 'PearsonCorr':
        method = 'correlation'
    elif method == 'PartialCorr':
        method = 'partial correlation'

    connectivity_measure = ConnectivityMeasure(kind=method)
    connectome = connectivity_measure.fit_transform([data])[0]

    file = os.path.abspath('./%s_connectome.npy' % (method.replace(" ", "-")))
    np.save(file, connectome)
    return file


def create_connectome(name='connectome'):

    wf = pe.Workflow(name=name)

    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'time_series',
            'method'
        ]),
        name='inputspec'
    )

    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'connectome',
        ]),
        name='outputspec'
    )

    node = pe.Node(Function(input_names=['time_series', 'method'],
                            output_names=['connectome'],
                            function=compute_correlation,
                            as_module=True),
                   name='connectome')

    wf.connect([
        (inputspec, node, [('time_series', 'time_series')]),
        (inputspec, node, [('method', 'method')]),
        (node, outputspec, [('connectome', 'connectome')]),
    ])

    return wf
    