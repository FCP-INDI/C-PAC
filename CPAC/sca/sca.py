# Copyright (C) 2022-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
from CPAC.pipeline.nodeblock import nodeblock
from nipype.interfaces.afni import preprocess
from CPAC.pipeline import nipype_pipeline_engine as pe
from nipype.interfaces import fsl, utility as util

from CPAC.sca.utils import *
# from CPAC.utils.utils import extract_one_d
from CPAC.utils.datasource import resample_func_roi, \
    create_roi_mask_dataflow, create_spatial_map_dataflow

from CPAC.timeseries.timeseries_analysis import get_roi_timeseries, \
    get_spatial_map_timeseries, resample_function


def create_sca(name_sca='sca'):
    """
    Map of the correlations of the Region of Interest(Seed in native or MNI space) with the rest of brain voxels.
    The map is normalized to contain Z-scores, mapped in standard space and treated with spatial smoothing.

    Parameters
    ----------
    name_sca : a string
        Name of the SCA workflow

    Returns
    -------
    sca_workflow : workflow
        Seed Based Correlation Analysis Workflow

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/sca/sca.py>`_

    Workflow Inputs::
        inputspec.rest_res_filt : string (existing nifti file)
            Band passed Image with Global Signal , white matter, csf and
            motion regression. Recommended bandpass filter (0.001,0.1) )

        inputspec.timeseries_one_d : string (existing nifti file)
            1D 3dTcorr1D compatible timeseries file. 1D file can be timeseries
            from a mask or from a parcellation containing ROIs

    Workflow Outputs::
        outputspec.correlation_file : string (nifti file)
            Correlations of the functional file and the input time series

        outputspec.Z_score : string (nifti file)
            Fisher Z transformed correlations of the seed

    SCA Workflow Procedure:

    1. Compute pearson correlation between input timeseries 1D file and input functional file
       Use 3dTcorr1D to compute that. Input timeseries can be a 1D file containing parcellation ROI's
       or a 3D mask

    2. Compute Fisher Z score of the correlation computed in step above. If a mask is provided then a
       a single Z score file is returned, otherwise z-scores for all ROIs are returned as a list of
       nifti files

    .. exec::
        from CPAC.sca import create_sca
        wf = create_sca()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/sca.dot'
        )

    Workflow:

    .. image:: ../../images/generated/sca.png
        :width: 500

    Detailed Workflow:

    .. image:: ../../images/generated/sca_detailed.png
        :width: 500

    Examples
    --------

    >>> sca_w = create_sca("sca_wf")
    >>> sca_w.inputs.inputspec.functional_file = '/home/data/subject/func/rest_bandpassed.nii.gz'  # doctest: +SKIP
    >>> sca_w.inputs.inputspec.timeseries_one_d = '/home/data/subject/func/ts.1D'  # doctest: +SKIP
    >>> sca_w.run() # doctest: +SKIP

    """

    from CPAC.utils.utils import get_roi_num_list

    sca = pe.Workflow(name=name_sca)
    inputNode = pe.Node(util.IdentityInterface(fields=['timeseries_one_d',
                                                       'functional_file',]),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'correlation_stack',
                                                    'correlation_files',
                                                    'Z_score',
                                                    ]),
                         name='outputspec')

    # 2. Compute voxel-wise correlation with Seed Timeseries
    corr = pe.Node(interface=preprocess.TCorr1D(),
                      name='3dTCorr1D', mem_gb=3.0)
    corr.inputs.pearson = True
    corr.inputs.outputtype = 'NIFTI_GZ'

    sca.connect(inputNode, 'timeseries_one_d',
                corr, 'y_1d')
    sca.connect(inputNode, 'functional_file',
                corr, 'xset')

    # Transform the sub-bricks into volumes
    try:
        concat = pe.Node(interface=preprocess.TCat(), name='3dTCat')
    except AttributeError:
        from nipype.interfaces.afni import utils as afni_utils
        concat = pe.Node(interface=afni_utils.TCat(), name='3dTCat')

    concat.inputs.outputtype = 'NIFTI_GZ'

    # also write out volumes as individual files
    #split = pe.Node(interface=fsl.Split(), name='split_raw_volumes_sca')
    #split.inputs.dimension = 't'
    #split.inputs.out_base_name = 'sca_'

    #get_roi_num_list = pe.Node(util.Function(input_names=['timeseries_file',
    #                                                      'prefix'],
    #                                         output_names=['roi_list'],
    #                                         function=get_roi_num_list),
    #                           name='get_roi_num_list')
    #get_roi_num_list.inputs.prefix = "sca"

    #sca.connect(inputNode, 'timeseries_one_d', get_roi_num_list,
    #            'timeseries_file')

    #rename_rois = pe.MapNode(interface=util.Rename(), name='output_rois',
    #                         iterfield=['in_file', 'format_string'])
    #rename_rois.inputs.keep_ext = True

    #sca.connect(split, 'out_files', rename_rois, 'in_file')
    #sca.connect(get_roi_num_list, 'roi_list', rename_rois, 'format_string')

    sca.connect(corr, 'out_file', concat, 'in_files')
    #sca.connect(concat, 'out_file', split, 'in_file')
    sca.connect(concat, 'out_file',
                outputNode, 'correlation_stack')
    #sca.connect(rename_rois, 'out_file', outputNode,
    #            'correlation_files')

    return sca


def create_temporal_reg(wflow_name='temporal_reg', which='SR'):
    """
    Temporal multiple regression workflow
    Provides a spatial map of parameter estimates corresponding to each
    provided timeseries in a timeseries.txt file as regressors

    Parameters
    ----------

    wflow_name : a string
        Name of the temporal regression workflow

    which: a string
        SR: Spatial Regression, RT: ROI Timeseries

        NOTE: If you set (which = 'RT'), the output of this workflow will be
        renamed based on the header information provided in the
        timeseries.txt file.
        If you run the temporal regression workflow manually, don\'t set
        (which = 'RT') unless you provide a timeseries.txt file with a header
        containing the names of the timeseries.

    Returns
    -------

    wflow : workflow

        temporal multiple regression Workflow



    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/sca/sca.py>`_

    Workflow Inputs::

        inputspec.subject_rest : string (existing nifti file)
            Band passed Image with Global Signal , white matter, csf and
            motion regression. Recommended bandpass filter (0.001,0.1) )

        inputspec.subject_timeseries : string (existing txt file)
            text file containing the timeseries to be regressed on the subjects
            functional file
            timeseries are organized by columns, timepoints by rows

        inputspec.subject_mask : string (existing nifti file)
            path to subject functional mask

        inputspec.demean : Boolean
            control whether to demean model and data

        inputspec.normalize : Boolean
            control whether to normalize the input timeseries to unit standard deviation



    Workflow Outputs::

        outputspec.temp_reg_map : string (nifti file)
            GLM parameter estimate image for each timeseries in the input file

        outputspec.temp_reg_map_zstat : string (nifti file)
            Normalized version of the GLM parameter estimates


    Temporal Regression Workflow Procedure:

    Enter all timeseries into a general linear model and regress these
    timeseries to the subjects functional file to get spatial maps of voxels
    showing activation patterns related to those in the timeseries.

    .. exec::
        from CPAC.sca import create_temporal_reg
        wf = create_temporal_reg()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/create_temporal_regression.dot'
        )

    Workflow:

    .. image:: ../../images/generated/create_temporal_regression.png
        :width: 500

    Detailed Workflow:

    .. image:: ../../images/generated/create_temporal_regression_detailed.png
        :width: 500

    References
    ----------
    `http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/DualRegression/UserGuide <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/DualRegression/UserGuide>`_

    Examples
    --------

    >>> tr_wf = create_temporal_reg('temporal-regression')
    >>> tr_wf.inputs.inputspec.subject_rest = '/home/data/subject/func/rest_bandpassed.nii.gz'  # doctest: +SKIP
    >>> tr_wf.inputs.inputspec.subject_timeseries = '/home/data/subject/func/timeseries.txt'  # doctest: +SKIP
    >>> tr_wf.inputs.inputspec.subject_mask = '/home/data/spatialmaps/spatial_map.nii.gz'  # doctest: +SKIP
    >>> tr_wf.inputs.inputspec.demean = True
    >>> tr_wf.inputs.inputspec.normalize = True
    >>> tr_wf.run() # doctest: +SKIP

    """

    wflow = pe.Workflow(name=wflow_name)

    inputNode = pe.Node(util.IdentityInterface
                        (fields=['subject_rest',
                                 'subject_timeseries',
                                 'subject_mask',
                                 'demean',
                                 'normalize']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface
                         (fields=['temp_reg_map',
                                  'temp_reg_map_files',
                                  'temp_reg_map_z',
                                  'temp_reg_map_z_files']),
                         name='outputspec')

    check_timeseries = pe.Node(util.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=check_ts),
                               name='check_timeseries')

    wflow.connect(inputNode, 'subject_timeseries',
                  check_timeseries, 'in_file')

    temporalReg = pe.Node(interface=fsl.GLM(), name='temporal_regression',
        mem_gb=4.0)
    temporalReg.inputs.out_file = 'temp_reg_map.nii.gz'
    temporalReg.inputs.out_z_name = 'temp_reg_map_z.nii.gz'

    wflow.connect(inputNode, 'subject_rest', temporalReg, 'in_file')
    wflow.connect(check_timeseries, 'out_file', temporalReg, 'design')
    wflow.connect(inputNode, 'demean', temporalReg, 'demean')
    wflow.connect(inputNode, 'normalize', temporalReg, 'des_norm')
    wflow.connect(inputNode, 'subject_mask', temporalReg, 'mask')

    wflow.connect(temporalReg, 'out_file', outputNode, 'temp_reg_map')
    wflow.connect(temporalReg, 'out_z', outputNode, 'temp_reg_map_z')

    '''
    split = pe.Node(interface=fsl.Split(), name='split_raw_volumes')
    split.inputs.dimension = 't'
    split.inputs.out_base_name = 'temp_reg_map_'

    wflow.connect(temporalReg, 'out_file', split, 'in_file')

    split_zstat = pe.Node(interface=fsl.Split(), name='split_zstat_volumes')
    split_zstat.inputs.dimension = 't'
    split_zstat.inputs.out_base_name = 'temp_reg_map_z_'

    wflow.connect(temporalReg, 'out_z',
                  split_zstat, 'in_file')

    if which == 'SR':
        wflow.connect(split, 'out_files',
                      outputNode, 'temp_reg_map_files')
        wflow.connect(split_zstat, 'out_files',
                      outputNode, 'temp_reg_map_z_files')

    elif which == 'RT':
        map_roi_imports = ['import os', 'import numpy as np']

        # get roi order and send to output node for raw outputs
        get_roi_order = pe.Node(util.Function(input_names=['maps',
                                                           'timeseries'],
                                              output_names=['labels',
                                                            'maps'],
                                              function=map_to_roi,
                                              imports=map_roi_imports),
                                name='get_roi_order')

        wflow.connect(split, 'out_files', get_roi_order, 'maps')

        wflow.connect(inputNode, 'subject_timeseries',
                      get_roi_order, 'timeseries')

        rename_maps = pe.MapNode(interface=util.Rename(),
                                 name='rename_maps',
                                 iterfield=['in_file',
                                            'format_string'])
        rename_maps.inputs.keep_ext = True

        wflow.connect(get_roi_order, 'labels', rename_maps, 'format_string')
        wflow.connect(get_roi_order, 'maps', rename_maps, 'in_file')
        wflow.connect(rename_maps, 'out_file',
                      outputNode, 'temp_reg_map_files')

        # get roi order and send to output node for z-stat outputs
        get_roi_order_zstat = pe.Node(util.Function(input_names=['maps',
                                                           'timeseries'],
                                                    output_names=['labels',
                                                                  'maps'],
                                                    function=map_to_roi,
                                                    imports=map_roi_imports),
                                      name='get_roi_order_zstat')

        wflow.connect(split_zstat, 'out_files', get_roi_order_zstat, 'maps')
        wflow.connect(inputNode, 'subject_timeseries',
                      get_roi_order_zstat, 'timeseries')

        rename_maps_zstat = pe.MapNode(interface=util.Rename(),
                                       name='rename_maps_zstat',
                                       iterfield=['in_file',
                                                  'format_string'])
        rename_maps_zstat.inputs.keep_ext = True

        wflow.connect(get_roi_order_zstat, 'labels',
                      rename_maps_zstat, 'format_string')
        wflow.connect(get_roi_order_zstat, 'maps',
                      rename_maps_zstat, 'in_file')

        wflow.connect(rename_maps_zstat, 'out_file',
                      outputNode, 'temp_reg_map_z_files')
    '''

    return wflow


@nodeblock(
    name="SCA_AVG",
    config=["seed_based_correlation_analysis"],
    switch=["run"],
    inputs=["space-template_desc-preproc_bold"],
    outputs=[
        "desc-MeanSCA_timeseries",
        "space-template_desc-MeanSCA_correlations",
        "atlas_name",
    ],
)
def SCA_AVG(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Run Seed-Based Correlation Analysis.'''

    # same workflow, except to run TSE and send it to the resource
    # pool so that it will not get sent to SCA
    resample_functional_roi_for_sca = pe.Node(
        util.Function(input_names=['in_func',
                                   'in_roi',
                                   'realignment',
                                   'identity_matrix'],
                      output_names=['out_func', 'out_roi'],
                      function=resample_func_roi,
                      as_module=True),
        name=f'resample_functional_roi_for_sca_{pipe_num}')

    resample_functional_roi_for_sca.inputs.realignment = \
        cfg.timeseries_extraction['realignment']
    resample_functional_roi_for_sca.inputs.identity_matrix = \
    cfg.registration_workflows['functional_registration'][
        'func_registration_to_template']['FNIRT_pipelines']['identity_matrix']

    roi_dataflow_for_sca = create_roi_mask_dataflow(
        cfg.seed_based_correlation_analysis['sca_atlases']['Avg'],
        f'roi_dataflow_for_sca_{pipe_num}'
    )

    roi_dataflow_for_sca.inputs.inputspec.set(
        creds_path=cfg.pipeline_setup['input_creds_path'],
        dl_dir=cfg.pipeline_setup['working_directory']['path']
    )

    roi_timeseries_for_sca = get_roi_timeseries(
        f'roi_timeseries_for_sca_{pipe_num}')

    node, out = strat_pool.get_data("space-template_desc-preproc_bold")
    # resample the input functional file to roi
    wf.connect(node, out,
                     resample_functional_roi_for_sca, 'in_func')
    wf.connect(roi_dataflow_for_sca, 'outputspec.out_file',
                     resample_functional_roi_for_sca, 'in_roi')

    # connect it to the roi_timeseries
    wf.connect(resample_functional_roi_for_sca, 'out_roi',
                     roi_timeseries_for_sca, 'input_roi.roi')
    wf.connect(resample_functional_roi_for_sca, 'out_func',
                     roi_timeseries_for_sca, 'inputspec.rest')

    sca_roi = create_sca(f'sca_roi_{pipe_num}')

    node, out = strat_pool.get_data("space-template_desc-preproc_bold")
    wf.connect(node, out, sca_roi, 'inputspec.functional_file')

    wf.connect(roi_timeseries_for_sca, 'outputspec.roi_csv',
               #('outputspec.roi_outputs', extract_one_d),
               sca_roi, 'inputspec.timeseries_one_d')

    outputs = {
        'desc-MeanSCA_timeseries':
            (roi_timeseries_for_sca, 'outputspec.roi_csv'),
                                    #('outputspec.roi_outputs',
                                    # extract_one_d)),
        'space-template_desc-MeanSCA_correlations':
            (sca_roi, 'outputspec.correlation_stack'),
        'atlas_name': (roi_dataflow_for_sca, 'outputspec.out_name')
    }

    return (wf, outputs)


@nodeblock(
    name="dual_regression",
    config=["seed_based_correlation_analysis"],
    switch=["run"],
    inputs=["space-template_desc-preproc_bold",
            "space-template_desc-bold_mask"],
    outputs=[
        "space-template_desc-DualReg_correlations",
        "desc-DualReg_statmap",
        "atlas_name",
    ],
)
def dual_regression(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    Run Dual Regression - spatial regression and then temporal regression.
    '''
    resample_spatial_map_to_native_space_for_dr = pe.Node(
        interface=fsl.FLIRT(),
        name=f'resample_spatial_map_to_native_space_for_DR_{pipe_num}'
    )
    resample_spatial_map_to_native_space_for_dr.inputs.set(
        interp='nearestneighbour',
        apply_xfm=True,
        in_matrix_file=
        cfg.registration_workflows['functional_registration'][
            'func_registration_to_template']['FNIRT_pipelines'][
            'identity_matrix']
    )

    spatial_map_dataflow_for_dr = create_spatial_map_dataflow(
        cfg.seed_based_correlation_analysis['sca_atlases']['DualReg'],
        f'spatial_map_dataflow_for_DR_{pipe_num}'
    )

    spatial_map_dataflow_for_dr.inputs.inputspec.set(
        creds_path=cfg.pipeline_setup['input_creds_path'],
        dl_dir=cfg.pipeline_setup['working_directory']['path']
    )

    spatial_map_timeseries_for_dr = get_spatial_map_timeseries(
        f'spatial_map_timeseries_for_DR_{pipe_num}'
    )
    spatial_map_timeseries_for_dr.inputs.inputspec.demean = True

    # resample the input functional file and functional mask
    # to spatial map
    node, out = strat_pool.get_data("space-template_desc-preproc_bold")
    wf.connect(node, out,
               resample_spatial_map_to_native_space_for_dr, 'reference')
    wf.connect(node, out,
               spatial_map_timeseries_for_dr, 'inputspec.subject_rest')

    wf.connect(spatial_map_dataflow_for_dr, 'select_spatial_map.out_file',
               resample_spatial_map_to_native_space_for_dr, 'in_file')

    # connect it to the spatial_map_timeseries
    wf.connect(resample_spatial_map_to_native_space_for_dr, 'out_file',
               spatial_map_timeseries_for_dr, 'inputspec.spatial_map'
    )

    dr_temp_reg = create_temporal_reg(f'temporal_regression_{pipe_num}')
    dr_temp_reg.inputs.inputspec.normalize = \
        cfg.seed_based_correlation_analysis['norm_timeseries_for_DR']
    dr_temp_reg.inputs.inputspec.demean = True

    wf.connect(spatial_map_timeseries_for_dr, 'outputspec.subject_timeseries',
               dr_temp_reg, 'inputspec.subject_timeseries')

    node, out = strat_pool.get_data("space-template_desc-preproc_bold")
    wf.connect(node, out, dr_temp_reg, 'inputspec.subject_rest')

    node, out = strat_pool.get_data("space-template_desc-bold_mask")
    wf.connect(node, out, dr_temp_reg, 'inputspec.subject_mask')

    outputs = {
        'space-template_desc-DualReg_correlations':
            (dr_temp_reg, 'outputspec.temp_reg_map'),
        'desc-DualReg_statmap':
            (dr_temp_reg, 'outputspec.temp_reg_map_z'),
        'atlas_name':
            (spatial_map_dataflow_for_dr, 'select_spatial_map.out_name')
    }

    return (wf, outputs)


@nodeblock(
    name="multiple_regression",
    config=["seed_based_correlation_analysis"],
    switch=["run"],
    inputs=["space-template_desc-preproc_bold",
            "space-template_desc-bold_mask"],
    outputs=[
        "space-template_desc-MultReg_correlations",
        "desc-MultReg_statmap",
        "atlas_name",
    ],
)
def multiple_regression(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Run Multiple Regression.'''

    # same workflow, except to run TSE and send it to the resource
    # pool so that it will not get sent to SCA
    resample_functional_roi_for_multreg = pe.Node(
        resample_function(),
        name=f'resample_functional_roi_for_multreg_{pipe_num}')

    resample_functional_roi_for_multreg.inputs.realignment = \
    cfg.timeseries_extraction['realignment']
    resample_functional_roi_for_multreg.inputs.identity_matrix = \
    cfg.registration_workflows['functional_registration'][
        'func_registration_to_template']['FNIRT_pipelines']['identity_matrix']

    roi_dataflow_for_multreg = create_roi_mask_dataflow(
        cfg.seed_based_correlation_analysis['sca_atlases']['MultReg'],
        f'roi_dataflow_for_mult_reg_{pipe_num}')

    roi_dataflow_for_multreg.inputs.inputspec.set(
        creds_path=cfg.pipeline_setup['input_creds_path'],
        dl_dir=cfg.pipeline_setup['working_directory']['path']
    )

    roi_timeseries_for_multreg = get_roi_timeseries(
        f'roi_timeseries_for_mult_reg_{pipe_num}')

    node, out = strat_pool.get_data("space-template_desc-preproc_bold")
    # resample the input functional file to roi
    wf.connect(node, out, resample_functional_roi_for_multreg, 'in_func')
    wf.connect(roi_dataflow_for_multreg,
                     'outputspec.out_file',
                     resample_functional_roi_for_multreg,
                     'in_roi')

    # connect it to the roi_timeseries
    wf.connect(resample_functional_roi_for_multreg,
                     'out_roi',
                     roi_timeseries_for_multreg,
                     'input_roi.roi')
    wf.connect(resample_functional_roi_for_multreg,
                     'out_func',
                     roi_timeseries_for_multreg,
                     'inputspec.rest')

    sc_temp_reg = create_temporal_reg(
        f'temporal_regression_sca_{pipe_num}',
        which='RT')
    sc_temp_reg.inputs.inputspec.normalize = \
    cfg.seed_based_correlation_analysis['norm_timeseries_for_DR']
    sc_temp_reg.inputs.inputspec.demean = True

    node, out = strat_pool.get_data(["space-template_desc-cleaned_bold",
                                     "space-template_desc-brain_bold",
                                     "space-template_desc-motion_bold",
                                     "space-template_desc-preproc_bold",
                                     "space-template_bold"])
    wf.connect(node, out, sc_temp_reg, 'inputspec.subject_rest')

    wf.connect(roi_timeseries_for_multreg, 'outputspec.roi_csv',
                     #('outputspec.roi_outputs', extract_one_d),
                     sc_temp_reg, 'inputspec.subject_timeseries')

    node, out = strat_pool.get_data('space-template_desc-bold_mask')
    wf.connect(node, out, sc_temp_reg, 'inputspec.subject_mask')

    outputs = {
        'space-template_desc-MultReg_correlations':
            (sc_temp_reg, 'outputspec.temp_reg_map'),
        'desc-MultReg_statmap':
            (sc_temp_reg, 'outputspec.temp_reg_map_z'),
        'atlas_name': (roi_dataflow_for_multreg, 'outputspec.out_name')
    }

    return (wf, outputs)
