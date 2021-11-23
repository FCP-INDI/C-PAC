#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions for creating connectome connectivity matrices."""
from warnings import warn
import numpy as np
from nilearn.connectome import ConnectivityMeasure
from nilearn.input_data import NiftiLabelsMasker
from nipype import logging
from nipype.interfaces import utility as util
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.interfaces.function import Function
from CPAC.utils.interfaces.netcorr import NetCorr, strip_afni_output_header

logger = logging.getLogger('nipype.workflow')
connectome_methods = {
    'afni': {'Pearson': '',
             'Partial': '-part_corr'},
    'nilearn': {'Pearson': 'correlation',
                'Partial': 'partial correlation'}
}


def connectome_name(timeseries, atlas_name, method):
    """Helper function to create connectome file filename

    Parameters
    ----------
    timeseries : str
        path to input timeseries

    atlas_name : str
        atlas name

    method : str
        BIDS entity value for `desc-` key

    Returns
    -------
    str
    """
    method = ''.join(word.capitalize() for word in method.split(' '))
    new_filename_parts = [part for part in timeseries.split('_')[:-1][::-1] if
                          not part.startswith('space-')]
    atlas_index = len(new_filename_parts) - 1
    if any(filename_part.startswith('desc-') for filename_part in
            new_filename_parts):
        for i, filename_part in enumerate(new_filename_parts):
            if filename_part.startswith('desc-'):
                new_filename_parts[-i] = (
                    filename_part + method.capitalize())
                atlas_index = -(i - 1)
                break
    new_filename_parts.insert(atlas_index, f'atlas-{atlas_name}')
    return '_'.join([*new_filename_parts[::-1], 'connectome.tsv'])


def get_connectome_method(method, tool):
    """Helper function to get tool's method string

    Parameters
    ----------
    method : str

    tool : str

    Returns
    -------
    str or NotImplemented

    Examples
    --------
    >>> get_connectome_method('Pearson', 'AFNI')
    ''
    >>> get_connectome_method('Pearson', 'Nilearn')
    'correlation'
    >>> get_connectome_method('Spearman', 'AFNI')
    NotImplemented
    """
    cm_method = connectome_methods[tool.lower()].get(method, NotImplemented)
    if cm_method is NotImplemented:
        warning_message = (
            f'{method} has not yet been implemented for {tool} in C-PAC.')
        if logger:
            logger.warning(NotImplementedError(warning_message))
        else:
            warn(warning_message, category=Warning)
    return cm_method


def compute_connectome_nilearn(in_rois, in_file, method, atlas_name):
    """Function to compute a connectome matrix using Nilearn

    Parameters
    ----------
    in_rois : Niimg-like object
        http://nilearn.github.io/manipulating_images/input_output.html#niimg-like-objects
        Region definitions, as one image of labels.

    in_file : str
        path to timeseries image

    method: str
        'Pearson' or 'Partial'

    atlas_name: str

    Returns
    -------
    numpy.ndarray or NotImplemented
    """
    output = connectome_name(in_file, atlas_name, f'Nilearn{method}')
    method = get_connectome_method(method, 'Nilearn')
    if method is NotImplemented:
        return NotImplemented
    masker = NiftiLabelsMasker(labels_img=in_rois,
                               standardize=True,
                               verbose=True)
    timeser = masker.fit_transform(in_file)
    correlation_measure = ConnectivityMeasure(kind=method)
    corr_matrix = correlation_measure.fit_transform([timeser])[0]
    np.fill_diagonal(corr_matrix, 0)
    np.savetxt(output, corr_matrix, delimiter='\t')
    return output


def create_connectome_afni(name, method, pipe_num):
    wf = pe.Workflow(name=name)
    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'in_rois',  # parcellation
            'in_file',  # timeseries,
            'mask',
            'method',
            'atlas_name'
        ]),
        name='inputspec'
    )
    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'out_file',
        ]),
        name='outputspec'
    )

    timeseries_correlation = pe.Node(NetCorr(), name=name)
    if method:
        timeseries_correlation.inputs.part_corr = (method == 'Partial')

    strip_header_node = pe.Node(Function(
        input_names=['in_file', 'out_file'], output_names=['out_file'],
        function=strip_afni_output_header),
                                name='netcorrStripHeader'
                                     f'{method}_{pipe_num}')

    name_output_node = pe.Node(Function(input_names=['timeseries',
                                                     'atlas_name',
                                                     'method'],
                                        output_names=['filename'],
                                        function=connectome_name),
                               name=f'connectomeName{method}_{pipe_num}')

    wf.connect([
        (inputspec, timeseries_correlation, [('in_rois', 'in_rois'),
                                             ('in_file', 'in_file'),
                                             ('mask', 'mask')]),
        (inputspec, name_output_node, [('in_file', 'timeseries'),
                                       ('atlas_name', 'atlas_name'),
                                       ('method', 'method')]),
        (name_output_node, strip_header_node, [('filename', 'out_file')]),
        (timeseries_correlation, strip_header_node, [
            ('out_file', 'in_file')]),
        (strip_header_node, outputspec, [('out_file', 'out_file')])])
    return wf


def create_connectome_nilearn(name='connectomeNilearn'):
    wf = pe.Workflow(name=name)
    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'in_rois',  # parcellation
            'in_file',  # timeseries
            'method',
            'atlas_name'
        ]),
        name='inputspec'
    )
    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'out_file',
        ]),
        name='outputspec'
    )
    node = pe.Node(Function(input_names=['in_rois', 'in_file', 'method',
                                         'atlas_name'],
                            output_names=['out_file'],
                            function=compute_connectome_nilearn,
                            as_module=True),
                   name='connectome')
    wf.connect([
        (inputspec, node, [('in_rois', 'in_rois'),
                           ('in_file', 'in_file'),
                           ('method', 'method'),
                           ('atlas_name', 'atlas_name')]),
        (node, outputspec, [('out_file', 'out_file')]),
    ])
    return wf
