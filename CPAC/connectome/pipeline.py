#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions for creating connectome connectivity matrices."""
import os
import re
import numpy as np
from nilearn.connectome import ConnectivityMeasure
from nipype.interfaces import utility as util
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.schema import valid_options
from CPAC.utils.interfaces.function import Function


def compute_correlation(time_series, method):
    """Function to create a numpy array file [1] containing a
    correlation matrix for a given time series and method.

    Parameters
    ----------
    time_series : str
        path to time series output

    method : str
        correlation matrix method. See https://fcp-indi.github.io/docs/nightly/user/tse#configuring-roi-time-series-extraction
        for options and how to configure.

    References
    ----------
    .. [1] The NumPy community (2021). NPY format. https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html#npy-format
    """  # pylint: disable=line-too-long  # noqa E501
    desc_regex = re.compile(r'\_desc.*(?=\_)')
    method_string = method.replace(' ', '+')

    data = np.genfromtxt(time_series).T

    if method == 'Pearson':
        method = 'correlation'
    elif method == 'Partial':
        method = 'partial correlation'
    # elif method == 'TangentEmbed':  # "Skip tangent embedding for now"
    #     method = 'tangent'

    connectivity_measure = ConnectivityMeasure(kind=method)
    connectome = connectivity_measure.fit_transform([data])[0]

    existing_desc = desc_regex.search(time_series)
    if hasattr(existing_desc, 'group'):
        matrix_file = time_series.replace(
            existing_desc.group(),
            '+'.join([
                existing_desc.group(),
                method_string
            ])
        ).replace('_timeseries.1D', '_connectome.npy')
    else:
        matrix_file = os.path.abspath(os.path.join(
            os.path.dirname(time_series),
            time_series.replace(
                '_timeseries.1D',
                f'_desc-{method_string}_connectome.npy'
            )
        ))

    np.save(matrix_file, connectome)
    return matrix_file


def create_connectome(name='connectome'):

    wf = pe.Workflow(name=name)

    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'timeseries',
            'measure'
        ]),
        name='inputspec'
    )

    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'connectome',
        ]),
        name='outputspec'
    )

    node = pe.Node(Function(input_names=['timeseries', 'measure'],
                            output_names=['connectome'],
                            function=compute_correlation,
                            as_module=True),
                   name='connectome')

    wf.connect([
        (inputspec, node, [('timeseries', 'timeseries')]),
        (inputspec, node, [('measure', 'measure')]),
        (node, outputspec, [('connectome', 'connectome')]),
    ])

    return wf


def timeseries_connectivity_matrix(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "timeseries_connectivity_matrix",
     "config": ["connectivity_matrix"],
     "switch": "None",
     "option_key": "using",
     "option_val": ["AFNI", "Nilearn"],
     "inputs": ["timeseries"],
     "outputs": ["connectome"]}
    '''
    # pylint: disable=invalid-name,unused-argument
    node, out = strat_pool.get_data(["timeseries"])
    outputs = {}
    for measure in cfg['connectivity_matrix', 'measure']:
        if measure in valid_options['connectivity_matrix']['measure']:
            timeseries_correlation = pe.Node(Function(
                input_names=['time_series', 'measure'],
                output_names=['connectome'],
                function=compute_correlation,
                as_module=True
            ), name=f'timeseries_{measure}_{pipe_num}')

            timeseries_correlation.inputs.measure = measure

            wf.connect(node, out,
                       timeseries_correlation, 'timeseries')
            outputs[f'desc-{measure}_connectome'] = (
                timeseries_correlation, 'outputspec.connectome'
            )

    return (wf, outputs)
