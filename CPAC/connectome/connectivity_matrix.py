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


def connectome_name(timeseries, atlas_name, tool, method):
    """Helper function to create connectome file filename

    Parameters
    ----------
    timeseries : str
        path to input timeseries

    atlas_name : str
        atlas name

    tool : str
        connectome tool

    method : str
        BIDS entity value for `desc-` key

    Returns
    -------
    str
    """
    method = ''.join(word.capitalize() for word in [tool, method])
    new_filename_parts = [part for part in timeseries.split('_')[:-1][::-1] if
                          not part.startswith('space-')]
    atlas_index = len(new_filename_parts) - 1
    if any(filename_part.startswith('desc-') for filename_part in
            new_filename_parts):
        for i, filename_part in enumerate(new_filename_parts):
            if filename_part.startswith('desc-'):
                new_filename_parts[-i] = f'desc-{method}'
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
    tool = 'Nilearn'
    output = connectome_name(in_file, atlas_name, tool, method)
    method = get_connectome_method(method, tool)
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
        imports=['import subprocess'],
        function=strip_afni_output_header),
                                name='netcorrStripHeader'
                                     f'{method}_{pipe_num}')

    name_output_node = pe.Node(Function(input_names=['timeseries',
                                                     'atlas_name',
                                                     'tool',
                                                     'method'],
                                        output_names=['filename'],
                                        function=connectome_name),
                               name=f'connectomeName{method}_{pipe_num}')
    name_output_node.inputs.tool = 'Afni'

    wf.connect([
        (inputspec, timeseries_correlation, [('in_rois', 'in_rois'),
                                             ('in_file', 'in_file'),
                                             ('mask', 'mask')]),
        (inputspec, name_output_node, [('in_file', 'timeseries'),
                                       ('atlas_name', 'atlas_name'),
                                       ('method', 'method')]),
        (timeseries_correlation, strip_header_node, [
            ('out_corr_matrix', 'in_file')]),
        (name_output_node, strip_header_node, [('filename', 'out_file')]),
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
<<<<<<< HEAD


def timeseries_connectivity_matrix(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "timeseries_connectivity_matrix",
     "config": ["connectivity_matrix"],
     "switch": "None",
     "option_key": "using",
     "option_val": ["AFNI", "Nilearn"],
     "inputs": ["timeseries", "method", "parcellation"],
     "outputs": ["connectome"]}
    '''  # pylint: disable=invalid-name,unused-argument
    outputs = {}
    for timeseries_analysis in cfg['timeseries_extraction', 'tse_roi_paths']:
        atlas = timeseries_analysis.split('/')[
            -1].split('.')[0].replace('_', '')
        analysis = tse_analyses.get(
            cfg['timeseries_extraction', 'tse_roi_paths', timeseries_analysis]
        )
        if analysis is None:
            continue
        node, out = strat_pool.get_data(f'atlas-{atlas}_desc-{analysis}_'
                                        'timeseries')
        for measure in cfg['connectivity_matrix', 'measure']:
            if measure in valid_options['connectivity_matrix']['measure']:
                if opt in ['AFNI', 'Nilearn']:
                    implementation = _get_method(measure, opt)
                    if implementation is NotImplemented:
                        continue
                    roi_dataflow = create_roi_mask_dataflow(
                            cfg.timeseries_extraction['tse_atlases']['Avg'],
                            f'roi_dataflow_{pipe_num}')
                    roi_dataflow.inputs.inputspec.set(
                        creds_path=cfg.pipeline_setup['input_creds_path'],
                        dl_dir=cfg.pipeline_setup['working_directory']['path']
                    )

                    if opt == 'Nilearn':
                        timeseries_correlation = pe.Node(Function(
                            input_names=['parcellation', 'timeseries'],
                            output_names=['connectomeNilearn'],
                            function=create_connectome_nilearn,
                            as_module=True
                        ), name=f'connectomeNilearn{measure}_{pipe_num}')
                        timeseries_correlation.inputs.measure = measure

                    elif opt == "AFNI":
                        timeseries_correlation = pe.Node(
                            NetCorr(),
                            name=f'connectomeAFNI{measure}_{pipe_num}')
                        if implementation:
                            timeseries_correlation.inputs.part_corr = (
                                measure == 'Partial'
                            )

                    wf.connect(roi_dataflow, 'outputspec.out_file',
                               timeseries_correlation, 'parcellation')

            else:
                timeseries_correlation = pe.Node(Function(
                    input_names=['timeseries', 'method'], #was measure
                    output_names=['connectome'],
                    function=create_connectome,
                    as_module=True
                ), name=f'connectome{measure}_{pipe_num}')

        wf.connect(node, out,
                   timeseries_correlation, 'timeseries')

    return (wf, outputs)
=======
>>>>>>> 7a952bbf9b88b085ce31d16bda59484de062c130
