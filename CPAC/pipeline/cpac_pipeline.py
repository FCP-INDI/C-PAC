# Copyright (C) 2012-2023  C-PAC Developers

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
import os
import sys
import time
import csv
import shutil
import pickle
import copy
import faulthandler

from time import strftime

import nipype
import yaml
# pylint: disable=wrong-import-order
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nipype_pipeline_engine.plugins import \
    LegacyMultiProcPlugin, MultiProcPlugin
from nipype import config, logging

from indi_aws import aws_utils, fetch_creds

import CPAC

from CPAC.pipeline.check_outputs import check_outputs
from CPAC.pipeline.engine import NodeBlock, initiate_rpool
from CPAC.anat_preproc.anat_preproc import (
    freesurfer_reconall,
    freesurfer_abcd_preproc,
    anatomical_init,
    acpc_align_head,
    acpc_align_head_with_mask,
    acpc_align_brain,
    acpc_align_brain_with_mask,
    registration_T2w_to_T1w,
    non_local_means,
    n4_bias_correction,
    t1t2_bias_correction,
    brain_mask_afni,
    brain_mask_fsl,
    brain_mask_niworkflows_ants,
    brain_mask_unet,
    brain_mask_freesurfer,
    brain_mask_freesurfer_abcd,
    brain_mask_freesurfer_fsl_tight,
    brain_mask_freesurfer_fsl_loose,
    brain_mask_acpc_afni,
    brain_mask_acpc_fsl,
    brain_mask_acpc_niworkflows_ants,
    brain_mask_acpc_unet,
    brain_mask_acpc_freesurfer,
    brain_mask_acpc_freesurfer_abcd,
    correct_restore_brain_intensity_abcd,
    brain_mask_acpc_freesurfer_fsl_tight,
    brain_mask_acpc_freesurfer_fsl_loose,
    brain_extraction_temp,
    brain_extraction,
    anatomical_init_T2,
    acpc_align_head_T2,
    acpc_align_head_with_mask_T2,
    acpc_align_brain_T2,
    acpc_align_brain_with_mask_T2,
    non_local_means_T2,
    n4_bias_correction_T2,
    brain_mask_T2,
    brain_mask_acpc_T2,
    brain_extraction_temp_T2,
    brain_extraction_T2,
)

from CPAC.registration.registration import (
    register_ANTs_anat_to_template,
    overwrite_transform_anat_to_template,
    register_FSL_anat_to_template,
    register_symmetric_ANTs_anat_to_template,
    register_symmetric_FSL_anat_to_template,
    register_ANTs_EPI_to_template,
    register_FSL_EPI_to_template,
    coregistration_prep_vol,
    coregistration_prep_mean,
    coregistration_prep_fmriprep,
    coregistration,
    create_func_to_T1template_xfm,
    create_func_to_T1template_symmetric_xfm,
    warp_wholeheadT1_to_template,
    warp_T1mask_to_template,
    apply_phasediff_to_timeseries_separately,
    apply_blip_to_timeseries_separately,
    warp_timeseries_to_T1template,
    warp_timeseries_to_T1template_deriv,
    warp_sbref_to_T1template,
    warp_bold_mask_to_T1template,
    warp_deriv_mask_to_T1template,
    warp_timeseries_to_EPItemplate,
    warp_bold_mean_to_EPItemplate,
    warp_bold_mask_to_EPItemplate,
    warp_deriv_mask_to_EPItemplate,
    warp_timeseries_to_T1template_abcd,
    single_step_resample_timeseries_to_T1template,
    warp_timeseries_to_T1template_dcan_nhp,
    warp_tissuemask_to_T1template,
    warp_tissuemask_to_EPItemplate
)

from CPAC.seg_preproc.seg_preproc import (
    tissue_seg_fsl_fast,
    tissue_seg_T1_template_based,
    tissue_seg_EPI_template_based,
    tissue_seg_ants_prior,
    tissue_seg_freesurfer
)

from CPAC.func_preproc import (
    calc_motion_stats,
    func_motion_correct,
    func_motion_correct_only,
    func_motion_estimates,
    get_motion_ref,
    motion_estimate_filter
)

from CPAC.func_preproc.func_preproc import (
    func_scaling,
    func_truncate,
    func_despike,
    func_despike_template,
    func_slice_time,
    func_reorient,
    bold_mask_afni,
    bold_mask_fsl,
    bold_mask_fsl_afni,
    bold_mask_anatomical_refined,
    bold_mask_anatomical_based,
    bold_mask_anatomical_resampled,
    bold_mask_ccs,
    bold_masking,
    func_mean,
    func_normalize
)

from CPAC.distortion_correction.distortion_correction import (
    distcor_phasediff_fsl_fugue,
    distcor_blip_afni_qwarp,
    distcor_blip_fsl_topup
)

from CPAC.nuisance.nuisance import (
    choose_nuisance_blocks,
    ICA_AROMA_ANTsreg,
    ICA_AROMA_FSLreg,
    ICA_AROMA_ANTsEPIreg,
    ICA_AROMA_FSLEPIreg,
    erode_mask_T1w,
    erode_mask_CSF,
    erode_mask_GM,
    erode_mask_WM,
    erode_mask_bold,
    erode_mask_boldCSF,
    erode_mask_boldGM,
    erode_mask_boldWM,
    nuisance_regression_template,
    ingress_regressors
)

from CPAC.surface.surf_preproc import surface_postproc
from CPAC.surface.surf_preproc import surface_falff
from CPAC.surface.surf_preproc import surface_alff
from CPAC.surface.surf_preproc import surface_reho
from CPAC.surface.surf_preproc import surface_connectivity_matrix

from CPAC.timeseries.timeseries_analysis import (
    timeseries_extraction_AVG,
    timeseries_extraction_Voxel,
    spatial_regression
)

from CPAC.sca.sca import (
    SCA_AVG,
    dual_regression,
    multiple_regression
)

from CPAC.alff.alff import alff_falff, alff_falff_space_template
from CPAC.reho.reho import reho, reho_space_template
from flowdump import save_workflow_json, WorkflowJSONMeta

from CPAC.utils.workflow_serialization import cpac_flowdump_serializer
from CPAC.vmhc.vmhc import (
    smooth_func_vmhc,
    warp_timeseries_to_sym_template,
    vmhc
)

from CPAC.network_centrality.pipeline import (
    network_centrality
)

from CPAC.pipeline.random_state import set_up_random_state_logger
from CPAC.pipeline.schema import valid_options
from CPAC.utils.trimmer import the_trimmer
from CPAC.utils import Configuration, set_subject
from CPAC.utils.docs import version_report
from CPAC.utils.versioning import REQUIREMENTS
from CPAC.qc.pipeline import create_qc_workflow
from CPAC.qc.xcp import qc_xcp

from CPAC.utils.monitoring import getLogger, log_nodes_cb, log_nodes_initial, \
                                  LOGTAIL, set_up_logger, \
                                  WARNING_FREESURFER_OFF_WITH_DATA
from CPAC.utils.monitoring.draw_gantt_chart import resource_report
from CPAC.utils.utils import (
    check_config_resources,
    check_system_deps,
)

logger = getLogger('nipype.workflow')
faulthandler.enable()

# config.enable_debug_mode()


def run_workflow(sub_dict, c, run, pipeline_timing_info=None, p_name=None,
                 plugin='MultiProc', plugin_args=None, test_config=False):
    '''
    Function to prepare and, optionally, run the C-PAC workflow

    Parameters
    ----------
    sub_dict : dictionary
        subject dictionary with anatomical and functional image paths
    c : Configuration object
        CPAC pipeline configuration dictionary object
    run : boolean
        flag to indicate whether to run the prepared workflow
    pipeline_timing_info : list (optional); default=None
        list of pipeline info for reporting timing information
    p_name : string (optional); default=None
        name of pipeline
    plugin : string (optional); default='MultiProc'
        nipype plugin to utilize when the workflow is ran
    plugin_args : dictionary (optional);
                  default={'status_callback': log_nodes_cb}
        plugin-specific arguments for the workflow plugin

    Returns
    -------
    workflow : nipype workflow
        the prepared nipype workflow object containing the parameters
        specified in the config
    '''
    from CPAC.utils.datasource import bidsier_prefix

    if plugin is not None and not isinstance(plugin, str):
        raise TypeError(
            'CPAC.pipeline.cpac_pipeline.run_workflow requires a '
            'string for the optional "plugin" argument, but a '
            f'{getattr(type(plugin), "__name__", str(type(plugin)))} '
            'was provided.')

    # Assure that changes on config will not affect other parts
    c = copy.copy(c)

    subject_id, p_name, log_dir = set_subject(sub_dict, c)
    c['subject_id'] = subject_id

    set_up_logger(f'{subject_id}_expectedOutputs',
                  filename=f'{bidsier_prefix(c["subject_id"])}_'
                           'expectedOutputs.yml',
                  level='info', log_dir=log_dir, mock=True,
                  overwrite_existing=True)
    if c.pipeline_setup['Debugging']['verbose']:
        set_up_logger('engine', level='debug', log_dir=log_dir, mock=True)

    config.update_config({
        'logging': {
            'log_directory': log_dir,
            'log_to_file': bool(getattr(c.pipeline_setup['log_directory'],
                                        'run_logging', True))
        },
        'execution': {
            'crashfile_format': 'txt',
            'resource_monitor_frequency': 0.2,
            'stop_on_first_crash': c['pipeline_setup', 'system_config',
                                     'fail_fast']}})
    config.enable_resource_monitor()
    logging.update_logging(config)

    # Start timing here
    pipeline_start_time = time.time()
    # at end of workflow, take timestamp again, take time elapsed and check
    # tempfile add time to time data structure inside tempfile, and increment
    # number of subjects

    # Check pipeline config resources
    (
        sub_mem_gb, num_cores_per_sub, num_ants_cores, num_omp_cores
    ) = check_config_resources(c)

    if not plugin_args:
        plugin_args = {}

    plugin_args['memory_gb'] = sub_mem_gb
    plugin_args['n_procs'] = num_cores_per_sub
    plugin_args['raise_insufficient'] = c['pipeline_setup', 'system_config',
                                          'raise_insufficient']
    plugin_args['status_callback'] = log_nodes_cb

    # perhaps in future allow user to set threads maximum
    # this is for centrality mostly
    # import mkl
    os.environ['OMP_NUM_THREADS'] = str(num_omp_cores)
    os.environ['MKL_NUM_THREADS'] = '1'  # str(num_cores_per_sub)
    os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = str(num_ants_cores)

    # TODO: TEMPORARY
    # TODO: solve the UNet model hanging issue during MultiProc
    if "UNet" in c.anatomical_preproc['brain_extraction']['using']:
        c.pipeline_setup['system_config']['max_cores_per_participant'] = 1
        logger.info("\n\n[!] LOCKING CPUs PER PARTICIPANT TO 1 FOR U-NET "
                    "MODEL.\n\nThis is a temporary measure due to a known "
                    "issue preventing Nipype's parallelization from running "
                    "U-Net properly.\n\n")

    # calculate maximum potential use of cores according to current pipeline
    # configuration
    max_core_usage = int(
        c.pipeline_setup['system_config']['max_cores_per_participant']) * \
                     int(c.pipeline_setup['system_config'][
                             'num_participants_at_once'])

    try:
        creds_path = sub_dict['creds_path']
        if creds_path and 'none' not in creds_path.lower():
            if os.path.exists(creds_path):
                input_creds_path = os.path.abspath(creds_path)
            else:
                err_msg = 'Credentials path: "%s" for subject "%s" was not ' \
                          'found. Check this path and try again.' % (
                              creds_path, subject_id)
                raise Exception(err_msg)
        else:
            input_creds_path = None
    except KeyError:
        input_creds_path = None

    # TODO enforce value with schema validation
    try:
        encrypt_data = bool(
            config.pipeline_setup['Amazon-AWS']['s3_encryption'])
    except:
        encrypt_data = False

    information = """
    Environment
    ===========
    {dependency_versions}

    Run command: {run_command}

    C-PAC version: {cpac_version}

    {license_notice}

    Setting maximum number of cores per participant to {cores}
    Setting number of participants at once to {participants}
    Setting OMP_NUM_THREADS to {omp_threads}
    Setting MKL_NUM_THREADS to 1
    Setting ANTS/ITK thread usage to {ants_threads}
    Maximum potential number of cores that might be used during this run: {max_cores}
{random_seed}
"""  # noqa: E501

    execution_info = """

    End of subject workflow {workflow}

    CPAC run complete:

        Pipeline configuration: {pipeline}
        Subject workflow: {workflow}
        Elapsed run time (minutes): {elapsed}
        Timing information saved in {log_dir}/cpac_individual_timing_{pipeline}.csv
        System time of start:      {run_start}
        System time of completion: {run_finish}
        {output_check}
"""  # noqa: E501

    logger.info('%s', information.format(
        run_command=' '.join(['run', *sys.argv[1:]]),
        cpac_version=CPAC.__version__,
        dependency_versions=version_report().replace('\n', '\n    '),
        cores=c.pipeline_setup['system_config']['max_cores_per_participant'],
        participants=c.pipeline_setup['system_config'][
            'num_participants_at_once'],
        omp_threads=c.pipeline_setup['system_config']['num_OMP_threads'],
        ants_threads=c.pipeline_setup['system_config']['num_ants_threads'],
        max_cores=max_core_usage,
        random_seed=(
            '    Random seed: %s' %
            c.pipeline_setup['system_config']['random_seed']) if
        c.pipeline_setup['system_config']['random_seed'] is not None else '',
        license_notice=CPAC.license_notice.replace('\n', '\n    ')))
    subject_info = {}
    subject_info['subject_id'] = subject_id
    subject_info['start_time'] = pipeline_start_time

    check_centrality_degree = c.network_centrality['run'] and \
                              (len(c.network_centrality['degree_centrality'][
                                       'weight_options']) != 0 or \
                               len(c.network_centrality[
                                       'eigenvector_centrality'][
                                       'weight_options']) != 0)

    check_centrality_lfcd = c.network_centrality['run'] and \
                            len(c.network_centrality[
                                    'local_functional_connectivity_density'][
                                    'weight_options']) != 0

    if not test_config:
        # Check system dependencies
        check_ica_aroma = c.nuisance_corrections['1-ICA-AROMA']['run']
        if isinstance(check_ica_aroma, list):
            check_ica_aroma = True in check_ica_aroma
        check_system_deps(check_ants='ANTS' in c.registration_workflows[
            'anatomical_registration']['registration']['using'],
                        check_ica_aroma=check_ica_aroma,
                        check_centrality_degree=check_centrality_degree,
                        check_centrality_lfcd=check_centrality_lfcd)

    # absolute paths of the dirs
    c.pipeline_setup['working_directory']['path'] = os.path.join(
        os.path.abspath(c.pipeline_setup['working_directory']['path']),
        p_name)
    if 's3://' not in c.pipeline_setup['output_directory']['path']:
        c.pipeline_setup['output_directory']['path'] = os.path.abspath(
            c.pipeline_setup['output_directory']['path'])

    if c.pipeline_setup['system_config']['random_seed'] is not None:
        set_up_random_state_logger(log_dir)

    try:
        workflow = build_workflow(
            subject_id, sub_dict, c, p_name, num_ants_cores
        )
    except Exception as exception:
        logger.exception('Building workflow failed')
        raise exception

    wf_graph = c['pipeline_setup', 'log_directory', 'graphviz',
                 'entire_workflow']
    if wf_graph.get('generate'):
        for graph2use in wf_graph.get('graph2use'):
            dotfilename = os.path.join(log_dir, f'{p_name}_{graph2use}.dot')
            for graph_format in wf_graph.get('format'):
                try:
                    workflow.write_graph(dotfilename=dotfilename,
                                         graph2use=graph2use,
                                         format=graph_format,
                                         simple_form=wf_graph.get(
                                             'simple_form', True))
                except Exception as exception:
                    raise RuntimeError(f'Failed to visualize {p_name} ('
                                       f'{graph2use}, {graph_format})'
                                       ) from exception

    workflow_meta = WorkflowJSONMeta(pipeline_name=p_name, stage='pre')
    save_workflow_json(
        filename=os.path.join(log_dir, workflow_meta.filename()),
        workflow=workflow,
        meta=workflow_meta,
        custom_serializer=cpac_flowdump_serializer
    )

    if test_config:
        logger.info('This has been a test of the pipeline configuration '
                    'file, the pipeline was built successfully, but was '
                    'not run')
    else:
        working_dir = os.path.join(
            c.pipeline_setup['working_directory']['path'], workflow.name)

        # if c.write_debugging_outputs:
        #    with open(os.path.join(working_dir, 'resource_pool.pkl'), 'wb') as f:
        #        pickle.dump(strat_list, f)

        # if c.pipeline_setup['working_directory']['regenerate_outputs'] is True:

        #     erasable = list(find_files(working_dir, '*sink*')) + \
        #         list(find_files(working_dir, '*link*')) + \
        #         list(find_files(working_dir, '*log*'))

        #     for f in erasable:
        #         if os.path.isfile(f):
        #             os.remove(f)
        #         else:
        #             shutil.rmtree(f)

        if hasattr(c, 'trim') and c.trim:
            logger.warning("""
Trimming is an experimental feature, and if used wrongly, it can
lead to unreproducible results.
It is useful for performance optimization, but only if used correctly.
Please, make yourself aware of how it works and its assumptions:
    - The pipeline configuration has not changed;
    - The data configuration / BIDS directory has not changed;
    - The files from the output directory has not changed;
    - Your softwares versions has not changed;
    - Your C-PAC version has not changed;
    - You do not have access to the working directory.
""")

            workflow, _ = the_trimmer(
                workflow,
                output_dir=c.pipeline_setup['output_directory']['path'],
                s3_creds_path=input_creds_path,
            )

        pipeline_start_datetime = strftime("%Y-%m-%d %H:%M:%S")

        workflow_result = None
        try:
            subject_info['resource_pool'] = []

            # for strat_no, strat in enumerate(strat_list):
            #    strat_label = 'strat_%d' % strat_no
            #    subject_info[strat_label] = strat.get_name()
            #    subject_info['resource_pool'].append(strat.get_resource_pool())

            subject_info['status'] = 'Running'

            # Create callback logger
            cb_log_filename = os.path.join(log_dir, 'callback.log')

            try:
                if not os.path.exists(os.path.dirname(cb_log_filename)):
                    os.makedirs(os.path.dirname(cb_log_filename))
            except IOError:
                pass

            # Add handler to callback log file
            set_up_logger('callback', cb_log_filename, 'debug', log_dir,
                          mock=True)

            # Log initial information from all the nodes
            log_nodes_initial(workflow)

            # Add status callback function that writes in callback log
            nipype_version = REQUIREMENTS['nipype']
            if nipype.__version__ != nipype_version:
                logger.warning('This version of Nipype may not be compatible '
                               f'with CPAC v{CPAC.__version__}, please '
                               f'install Nipype version {nipype_version}\n')

            if plugin_args['n_procs'] == 1:
                plugin = 'Linear'
            if not plugin or plugin == 'LegacyMultiProc':
                plugin = LegacyMultiProcPlugin(plugin_args)
            elif plugin == 'MultiProc':
                plugin = MultiProcPlugin(plugin_args)

            try:
                # Actually run the pipeline now, for the current subject
                workflow_result = workflow.run(plugin=plugin,
                                               plugin_args=plugin_args)
            except UnicodeDecodeError:
                raise EnvironmentError(
                    "C-PAC migrated from Python 2 to Python 3 in v1.6.2 (see "
                    "release notes). Your working directory contains Python 2 "
                    "pickles, probably from an older version of C-PAC. If you "
                    "want to continue to use this working directory, run\n\n"
                    "docker run -i --rm --user $(id -u):$(id -g) "
                    "-v /path/to/working_dir:/working "
                    "fcpindi/c-pac:latest /bids_dir /outputs cli -- "
                    "utils repickle /working\n"
                    "\nor\n\n"
                    "singularity run "
                    "C-PAC_latest.sif /bids_dir /outputs cli -- "
                    "utils repickle /path/to/working_dir\n\n"
                    "before running C-PAC >=v1.6.2"
                )

            # PyPEER kick-off
            # if c.PyPEER['run']:
            #    from CPAC.pypeer.peer import prep_for_pypeer
            #    prep_for_pypeer(c.PyPEER['eye_scan_names'], c.PyPEER['data_scan_names'],
            #                    c.PyPEER['eye_mask_path'], c.pipeline_setup['output_directory']['path'], subject_id,
            #                    pipeline_ids, c.PyPEER['stimulus_path'], c.PyPEER['minimal_nuisance_correction']['peer_gsr'],
            #                    c.PyPEER['minimal_nuisance_correction']['peer_scrub'], c.PyPEER['minimal_nuisance_correction']['scrub_thresh'])

            # Dump subject info pickle file to subject log dir
            subject_info['status'] = 'Completed'

            subject_info_file = os.path.join(
                log_dir, 'subject_info_%s.pkl' % subject_id
            )
            with open(subject_info_file, 'wb') as info:
                pickle.dump(list(subject_info), info)

            # have this check in case the user runs cpac_runner from terminal and
            # the timing parameter list is not supplied as usual by the GUI
            if pipeline_timing_info is not None:

                # pipeline_timing_info list:
                #  [0] - unique pipeline ID
                #  [1] - pipeline start time stamp (first click of 'run' from GUI)
                #  [2] - number of subjects in subject list
                unique_pipeline_id = pipeline_timing_info[0]
                pipeline_start_stamp = pipeline_timing_info[1]
                num_subjects = pipeline_timing_info[2]

                # elapsed time data list:
                #  [0] - elapsed time in minutes
                elapsed_time_data = []

                elapsed_time_data.append(
                    int(((time.time() - pipeline_start_time) / 60)))

                # elapsedTimeBin list:
                #  [0] - cumulative elapsed time (minutes) across all subjects
                #  [1] - number of times the elapsed time has been appended
                #        (effectively a measure of how many subjects have run)

                # TODO
                # write more doc for all this
                # warning in .csv that some runs may be partial
                # code to delete .tmp file

                timing_temp_file_path = os.path.join(
                    c.pipeline_setup['log_directory']['path'],
                    '%s_pipeline_timing.tmp' % unique_pipeline_id)

                if not os.path.isfile(timing_temp_file_path):
                    elapsedTimeBin = []
                    elapsedTimeBin.append(0)
                    elapsedTimeBin.append(0)

                    with open(timing_temp_file_path, 'wb') as handle:
                        pickle.dump(elapsedTimeBin, handle)

                with open(timing_temp_file_path, 'rb') as handle:
                    elapsedTimeBin = pickle.loads(handle.read())

                elapsedTimeBin[0] = elapsedTimeBin[0] + elapsed_time_data[0]
                elapsedTimeBin[1] = elapsedTimeBin[1] + 1

                with open(timing_temp_file_path, 'wb') as handle:
                    pickle.dump(elapsedTimeBin, handle)

                # this happens once the last subject has finished running!
                if elapsedTimeBin[1] == num_subjects:

                    pipelineTimeDict = {}
                    pipelineTimeDict['Pipeline'] = c.pipeline_setup[
                        'pipeline_name']
                    pipelineTimeDict['Cores_Per_Subject'] = \
                    c.pipeline_setup['system_config'][
                        'max_cores_per_participant']
                    pipelineTimeDict['Simultaneous_Subjects'] = \
                    c.pipeline_setup['system_config'][
                        'num_participants_at_once']
                    pipelineTimeDict['Number_of_Subjects'] = num_subjects
                    pipelineTimeDict['Start_Time'] = pipeline_start_stamp
                    pipelineTimeDict['End_Time'] = strftime(
                        "%Y-%m-%d_%H:%M:%S")
                    pipelineTimeDict['Elapsed_Time_(minutes)'] = \
                    elapsedTimeBin[0]
                    pipelineTimeDict['Status'] = 'Complete'

                    gpaTimeFields = [
                        'Pipeline', 'Cores_Per_Subject',
                        'Simultaneous_Subjects',
                        'Number_of_Subjects', 'Start_Time',
                        'End_Time', 'Elapsed_Time_(minutes)',
                        'Status'
                    ]
                    timeHeader = dict(zip(gpaTimeFields, gpaTimeFields))

                    with open(os.path.join(
                            c.pipeline_setup['log_directory']['path'],
                                    'cpac_individual_timing_%s.csv' %
                                    c.pipeline_setup['pipeline_name']
                    ), 'a') as timeCSV, open(os.path.join(
                        c.pipeline_setup['log_directory']['path'],
                                'cpac_individual_timing_%s.csv' %
                                c.pipeline_setup['pipeline_name']
                    ), 'r') as readTimeCSV:

                        timeWriter = csv.DictWriter(timeCSV,
                                                    fieldnames=gpaTimeFields)
                        timeReader = csv.DictReader(readTimeCSV)

                        headerExists = False
                        for line in timeReader:
                            if 'Start_Time' in line:
                                headerExists = True

                        if headerExists is False:
                            timeWriter.writerow(timeHeader)

                        timeWriter.writerow(pipelineTimeDict)

                    # remove the temp timing file now that it is no longer needed
                    os.remove(timing_temp_file_path)

            # Upload logs to s3 if s3_str in output directory
            if c.pipeline_setup['output_directory'][
                'path'].lower().startswith('s3://'):

                try:
                    # Store logs in s3 output director/logs/...
                    s3_log_dir = os.path.join(
                        c.pipeline_setup['output_directory']['path'],
                        'logs',
                        os.path.basename(log_dir)
                    )
                    bucket_name = \
                    c.pipeline_setup['output_directory']['path'].split('/')[2]
                    bucket = fetch_creds.return_bucket(creds_path,
                                                       bucket_name)

                    # Collect local log files
                    local_log_files = []
                    for root, _, files in os.walk(log_dir):
                        local_log_files.extend([os.path.join(root, fil)
                                                for fil in files])
                    # Form destination keys
                    s3_log_files = [loc.replace(log_dir, s3_log_dir)
                                    for loc in local_log_files]
                    # Upload logs
                    aws_utils.s3_upload(bucket,
                                        (local_log_files, s3_log_files),
                                        encrypt=encrypt_data)
                    # Delete local log files
                    for log_f in local_log_files:
                        os.remove(log_f)

                except Exception as exc:
                    err_msg = (
                        'Unable to upload CPAC log files in: %s.\nError: %s')
                    logger.error(err_msg, log_dir, exc)

        except Exception:
            import traceback
            traceback.print_exc()
            execution_info = """

Error of subject workflow {workflow}

CPAC run error:

    Pipeline configuration: {pipeline}
    Subject workflow: {workflow}
    Elapsed run time (minutes): {elapsed}
    Timing information saved in {log_dir}/cpac_individual_timing_{pipeline}.csv
    System time of start:      {run_start}
    {output_check}
"""

        finally:

            if workflow:
                if os.path.exists(cb_log_filename):
                    resource_report(cb_log_filename,
                                    num_cores_per_sub, logger)

                logger.info('%s', execution_info.format(
                    workflow=workflow.name,
                    pipeline=c.pipeline_setup['pipeline_name'],
                    log_dir=c.pipeline_setup['log_directory']['path'],
                    elapsed=(time.time() - pipeline_start_time) / 60,
                    run_start=pipeline_start_datetime,
                    run_finish=strftime("%Y-%m-%d %H:%M:%S"),
                    output_check=check_outputs(
                                 c.pipeline_setup['output_directory']['path'],
                                 log_dir, c.pipeline_setup['pipeline_name'],
                                 c['subject_id'])
                ))

                if workflow_result is not None:
                    workflow_meta.stage = "post"
                    save_workflow_json(
                        filename=os.path.join(log_dir,
                                              workflow_meta.filename()),
                        workflow=workflow_result,
                        meta=workflow_meta,
                        custom_serializer=cpac_flowdump_serializer
                    )

                # Remove working directory when done
                if c.pipeline_setup['working_directory'][
                    'remove_working_dir']:
                    remove_workdir(working_dir)
                # Remove just .local from working directory
                else:
                    remove_workdir(os.path.join(os.environ["CPAC_WORKDIR"],
                                                '.local'))


def remove_workdir(wdpath: str) -> None:
    """Remove a given working directory if possible, warn if impossible

    Parameters
    ----------
    wdpath : str
        path to working directory to remove
    """
    try:
        if os.path.exists(wdpath):
            logger.info("Removing working dir: %s", wdpath)
            shutil.rmtree(wdpath)
    except (FileNotFoundError, PermissionError):
        logger.warning(
            'Could not remove working directory %s', wdpath)


def initialize_nipype_wf(cfg, sub_data_dct, name=""):

    if name:
        name = f'_{name}'

    workflow_name = f'cpac{name}_{sub_data_dct["subject_id"]}_{sub_data_dct["unique_id"]}'
    wf = pe.Workflow(name=workflow_name)
    wf.base_dir = cfg.pipeline_setup['working_directory']['path']
    wf.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(cfg.pipeline_setup['log_directory'][
                                         'path'])
    }

    return wf


def load_cpac_pipe_config(pipe_config):
    # Load in pipeline config file
    config_file = os.path.realpath(pipe_config)
    try:
        if not os.path.exists(config_file):
            raise IOError
        else:
            cfg = Configuration(yaml.safe_load(open(config_file, 'r')))
    except IOError:
        print("config file %s doesn't exist" % config_file)
        raise
    except yaml.parser.ParserError as e:
        error_detail = "\"%s\" at line %d" % (
            e.problem,
            e.problem_mark.line
        )
        raise Exception(
            "Error parsing config file: {0}\n\n"
            "Error details:\n"
            "    {1}"
            "\n\n".format(config_file, error_detail)
        )
    except Exception as e:
        raise Exception(
            "Error parsing config file: {0}\n\n"
            "Error details:\n"
            "    {1}"
            "\n\n".format(config_file, e)
        )
    return cfg


def build_anat_preproc_stack(rpool, cfg, pipeline_blocks=None):

    if not pipeline_blocks:
        pipeline_blocks = []

    # T1w Anatomical Preprocessing
    if not rpool.check_rpool('desc-reorient_T1w') and \
        not rpool.check_rpool('desc-preproc_T1w'):
        anat_init_blocks = [
            anatomical_init
        ]
        pipeline_blocks += anat_init_blocks
        
    using_brain_extraction = cfg.anatomical_preproc['brain_extraction']['using']

    if not rpool.check_rpool('freesurfer-subject-dir') and 'FreeSurfer-ABCD' not in using_brain_extraction:
        pipeline_blocks += [freesurfer_reconall]  # includes postproc

    if not rpool.check_rpool('desc-preproc_T1w'):

        # brain masking for ACPC alignment
        if cfg.anatomical_preproc['acpc_alignment']['acpc_target'] == 'brain':
                acpc_blocks = [
                    [brain_mask_acpc_afni,
                     brain_mask_acpc_fsl,
                     brain_mask_acpc_niworkflows_ants,
                     brain_mask_acpc_unet,
                     brain_mask_acpc_freesurfer_abcd,
                     brain_mask_acpc_freesurfer,
                     brain_mask_acpc_freesurfer_fsl_tight,
                     brain_mask_acpc_freesurfer_fsl_loose],
                    acpc_align_brain_with_mask,
                    brain_extraction_temp,
                    acpc_align_brain
                ]
        elif cfg.anatomical_preproc['acpc_alignment']['acpc_target'] == 'whole-head':
            if (rpool.check_rpool('space-T1w_desc-brain_mask') and \
                cfg.anatomical_preproc['acpc_alignment']['align_brain_mask']):
                acpc_blocks = [
                    acpc_align_head_with_mask
                    # outputs space-T1w_desc-brain_mask for later - keep the mask (the user provided)
                ]
            else:
                acpc_blocks = [
                    acpc_align_head  # does not output nor generate a mask
                ]
        

        anat_preproc_blocks = [
            (non_local_means, ('T1w', ['desc-preproc_T1w',
                                       'desc-reorient_T1w',
                                       'T1w'])),
            n4_bias_correction
        ]
        if cfg.anatomical_preproc['acpc_alignment']['run_before_preproc']:
            anat_blocks = acpc_blocks + anat_preproc_blocks
        else:
            anat_blocks = anat_preproc_blocks + acpc_blocks

        pipeline_blocks += anat_blocks

        pipeline_blocks += [freesurfer_abcd_preproc]

    if not rpool.check_rpool('freesurfer-subject-dir') and 'FreeSurfer-ABCD' in using_brain_extraction:
        pipeline_blocks += [freesurfer_reconall]  # includes postproc

    # Anatomical T1 brain masking

    anat_brain_mask_blocks = [
        [brain_mask_afni,
        brain_mask_fsl,
        brain_mask_niworkflows_ants,
        brain_mask_unet,
        brain_mask_freesurfer_abcd,
        brain_mask_freesurfer,
        brain_mask_freesurfer_fsl_tight,
        brain_mask_freesurfer_fsl_loose]
    ]
    pipeline_blocks += anat_brain_mask_blocks

    # T2w Anatomical Preprocessing
    if rpool.check_rpool('T2w'):
        if not rpool.check_rpool('desc-reorient_T2w'):
            anat_init_blocks_T2 = [
                anatomical_init_T2
            ]
            pipeline_blocks += anat_init_blocks_T2

        # TODO: T2 freesurfer_preproc?
        # pipeline_blocks += [freesurfer_preproc]

        if not rpool.check_rpool('desc-preproc_T2w'):

            # brain masking for ACPC alignment
            if cfg.anatomical_preproc['acpc_alignment']['acpc_target'] == 'brain':
                if rpool.check_rpool('space-T2w_desc-brain_mask'):
                    acpc_blocks_T2 = [
                        brain_extraction_temp_T2,
                        acpc_align_brain_with_mask_T2
                        # outputs space-T2w_desc-brain_mask for later - keep the mask (the user provided)
                    ]
                else:
                    acpc_blocks_T2 = [
                        brain_mask_acpc_T2,
                        # we don't want these masks to be used later, only used in brain_extraction_temp_T2
                        brain_extraction_temp_T2,
                        acpc_align_brain_T2
                    ]
            elif cfg.anatomical_preproc['acpc_alignment'][
                'acpc_target'] == 'whole-head':
                if rpool.check_rpool('space-T2w_desc-brain_mask'):
                    acpc_blocks_T2 = [
                        acpc_align_head_with_mask_T2
                        # outputs space-T2w_desc-brain_mask for later - keep the mask (the user provided)
                    ]
                else:
                    acpc_blocks_T2 = [
                        acpc_align_head_T2  # does not output nor generate a mask
                    ]

            anat_preproc_blocks_T2 = [
                registration_T2w_to_T1w,
                non_local_means_T2,
                n4_bias_correction_T2,
                t1t2_bias_correction
            ]
            if cfg.anatomical_preproc['acpc_alignment']['run_before_preproc']:
                anat_blocks_T2 = acpc_blocks_T2 + anat_preproc_blocks_T2
            else:
                anat_blocks_T2 = anat_preproc_blocks_T2 + acpc_blocks_T2

            pipeline_blocks += anat_blocks_T2

    # Anatomical T1 brain extraction
    if not rpool.check_rpool('desc-brain_T1w'):
        anat_brain_blocks = [
            brain_extraction
        ]
        pipeline_blocks += anat_brain_blocks

    # T2 brain masking
    if not rpool.check_rpool('space-T2w_desc-brain_mask'):
        anat_brain_mask_blocks_T2 = [
            brain_mask_T2
        ]
        pipeline_blocks += anat_brain_mask_blocks_T2

    if not rpool.check_rpool('desc-brain_T2w'):
        anat_brain_blocks_T2 = [
            brain_extraction_T2
        ]
        pipeline_blocks += anat_brain_blocks_T2

    return pipeline_blocks


def build_T1w_registration_stack(rpool, cfg, pipeline_blocks=None):

    if not pipeline_blocks:
        pipeline_blocks = []

    reg_blocks = []
    if not rpool.check_rpool('from-T1w_to-template_mode-image_xfm'):
        reg_blocks = [
            [register_ANTs_anat_to_template, register_FSL_anat_to_template],
            overwrite_transform_anat_to_template,
            warp_wholeheadT1_to_template,
            warp_T1mask_to_template
        ]

    if not rpool.check_rpool('desc-restore-brain_T1w'):
        reg_blocks.append(correct_restore_brain_intensity_abcd)

    if cfg.voxel_mirrored_homotopic_connectivity['run']:
        if not rpool.check_rpool('from-T1w_to-symtemplate_mode-image_xfm'):
            reg_blocks.append([register_symmetric_ANTs_anat_to_template,
                               register_symmetric_FSL_anat_to_template])
    pipeline_blocks += reg_blocks

    return pipeline_blocks


def build_segmentation_stack(rpool, cfg, pipeline_blocks=None):

    if not pipeline_blocks:
        pipeline_blocks = []

    if not rpool.check_rpool('label-CSF_mask') or \
            not rpool.check_rpool('label-WM_mask'):
        seg_blocks = [
            [tissue_seg_fsl_fast,
             tissue_seg_ants_prior,
	     tissue_seg_freesurfer]
        ]
        if 'T1_Template' in cfg.segmentation['tissue_segmentation'][
            'Template_Based']['template_for_segmentation']:
            seg_blocks.append(tissue_seg_T1_template_based)
        if 'EPI_Template' in cfg.segmentation['tissue_segmentation'][
            'Template_Based']['template_for_segmentation']:
            seg_blocks.append(tissue_seg_EPI_template_based)

        pipeline_blocks += seg_blocks

    if cfg.registration_workflows['anatomical_registration']['run'] and 'T1_Template' in cfg.segmentation[
    'tissue_segmentation']['Template_Based']['template_for_segmentation']:
        pipeline_blocks.append(warp_tissuemask_to_T1template)


    return pipeline_blocks


def list_blocks(pipeline_blocks, indent=None):
    """Function to list node blocks line by line

    Parameters
    ----------
    pipeline_blocks : list or tuple

    indent : int or None
       number of spaces after a tab indent

    Returns
    -------
    str
    """
    blockstring = yaml.dump([
        getattr(block, '__name__', getattr(block, 'name', yaml.safe_load(
            list_blocks(list(block))) if
            isinstance(block, (tuple, list, set)) else str(block))
        ) for block in pipeline_blocks])
    if isinstance(indent, int):
        blockstring = '\n'.join([
            '\t' + ' ' * indent + line.replace('- - ', '- ') for
            line in blockstring.split('\n')])
    return blockstring


def connect_pipeline(wf, cfg, rpool, pipeline_blocks):
    logger.info('\n'.join([
        'Connecting pipeline blocks:',
        list_blocks(pipeline_blocks, indent=1)]))

    previous_nb = None
    for block in pipeline_blocks:
        try:
            nb = NodeBlock(block, debug=cfg['pipeline_setup', 'Debugging',
                                            'verbose'])
            wf = nb.connect_block(wf, cfg, rpool)
        except LookupError as e:
            if nb.name == 'freesurfer_postproc':
                logger.warning(WARNING_FREESURFER_OFF_WITH_DATA)
                LOGTAIL['warnings'].append(WARNING_FREESURFER_OFF_WITH_DATA)
                continue
            previous_nb_str = (
                f"after node block '{previous_nb.get_name()}':"
            ) if previous_nb else 'at beginning:'
            # Alert user to block that raises error
            if isinstance(block, list):
                node_block_names = str([NodeBlock(b).get_name() for b in block])
                e.args = (
                    f'When trying to connect one of the node blocks '
                    f"{node_block_names} "
                    f"to workflow '{wf}' {previous_nb_str} {e.args[0]}",
                )
            else:
                node_block_names = NodeBlock(block).get_name()
                e.args = (
                    f'When trying to connect node block '
                    f"'{node_block_names}' "
                    f"to workflow '{wf}' {previous_nb_str} {e.args[0]}",
                )
            if cfg.pipeline_setup['Debugging']['verbose']:
                verbose_logger = getLogger('engine')
                verbose_logger.debug(e.args[0])
                verbose_logger.debug(rpool)
            raise
        previous_nb = nb

    return wf


def build_workflow(subject_id, sub_dict, cfg, pipeline_name=None,
                   num_ants_cores=1):
    from CPAC.utils.datasource import gather_extraction_maps

    # Workflow setup
    wf = initialize_nipype_wf(cfg, sub_dict)

    # Extract credentials path if it exists
    try:
        creds_path = sub_dict['creds_path']
        if creds_path and 'none' not in creds_path.lower():
            if os.path.exists(creds_path):
                input_creds_path = os.path.abspath(creds_path)
            else:
                err_msg = 'Credentials path: "%s" for subject "%s" was not ' \
                          'found. Check this path and try again.' % (
                              creds_path, subject_id)
                raise Exception(err_msg)
        else:
            input_creds_path = None
    except KeyError:
        input_creds_path = None

    cfg.pipeline_setup['input_creds_path'] = input_creds_path

    """""""""""""""""""""""""""""""""""""""""""""""""""
     PREPROCESSING
    """""""""""""""""""""""""""""""""""""""""""""""""""

    wf, rpool = initiate_rpool(wf, cfg, sub_dict)

    pipeline_blocks = build_anat_preproc_stack(rpool, cfg)

    # Anatomical to T1 template registration
    pipeline_blocks = build_T1w_registration_stack(rpool, cfg,
                                                   pipeline_blocks)

    # Anatomical tissue segmentation
    pipeline_blocks = build_segmentation_stack(rpool, cfg, pipeline_blocks)

    # Functional Preprocessing, including motion correction and BOLD masking
    if cfg.functional_preproc['run']:
        func_init_blocks = [
            func_reorient,
            func_scaling,
            func_truncate
        ]
        func_preproc_blocks = [
            func_despike,
            func_slice_time
        ]

        if not rpool.check_rpool('desc-mean_bold'):
            func_preproc_blocks.append(func_mean)

        func_mask_blocks = []
        if not rpool.check_rpool('space-bold_desc-brain_mask'):
            func_mask_blocks = [
                [bold_mask_afni, bold_mask_fsl, bold_mask_fsl_afni,
                 bold_mask_anatomical_refined, bold_mask_anatomical_based,
                 bold_mask_anatomical_resampled, bold_mask_ccs],
                bold_masking]

        func_prep_blocks = [
            calc_motion_stats,
            func_normalize,
            [coregistration_prep_vol,
             coregistration_prep_mean,
             coregistration_prep_fmriprep]
        ]

        # Distortion/Susceptibility Correction
        distcor_blocks = []
        if 'fmap' in sub_dict:
            fmap_keys = sub_dict['fmap']
            if 'phasediff' in fmap_keys or 'phase1' in fmap_keys:
                if 'magnitude' in fmap_keys or 'magnitude1' in fmap_keys:
                    distcor_blocks.append(distcor_phasediff_fsl_fugue)
            if len(fmap_keys) == 2:
                for key in fmap_keys:
                    if 'epi_' not in key:
                        break
                else:
                    distcor_blocks.append(distcor_blip_afni_qwarp) 
                    distcor_blocks.append(distcor_blip_fsl_topup)

        if distcor_blocks:
            if len(distcor_blocks) > 1:
                distcor_blocks = [distcor_blocks]
            func_prep_blocks += distcor_blocks

        func_motion_blocks = []
        if not rpool.check_rpool('desc-movementParameters_motion'):
            if cfg['functional_preproc']['motion_estimates_and_correction'][
                'motion_estimates']['calculate_motion_first']:
                func_motion_blocks = [
                    get_motion_ref,
                    func_motion_estimates,
                    motion_estimate_filter
                ]
                func_blocks = func_init_blocks + func_motion_blocks + \
                              func_preproc_blocks + [func_motion_correct_only] + \
                              func_mask_blocks + func_prep_blocks
            else:
                func_motion_blocks = [
                    get_motion_ref,
                    func_motion_correct,
                    motion_estimate_filter
                ]
                func_blocks = func_init_blocks + func_preproc_blocks + \
                              func_motion_blocks + func_mask_blocks + \
                              func_prep_blocks
        else:
            func_blocks = func_init_blocks + func_preproc_blocks + \
                          func_motion_blocks + func_mask_blocks + \
                          func_prep_blocks

        pipeline_blocks += func_blocks

    # BOLD to T1 coregistration
    if cfg.registration_workflows['functional_registration'][
        'coregistration']['run'] and \
            (not rpool.check_rpool('space-T1w_sbref') or
             not rpool.check_rpool('from-bold_to-T1w_mode-image_desc-linear_xfm')):
        coreg_blocks = [
            coregistration
        ]
        pipeline_blocks += coreg_blocks

    # BOLD to EPI-template registration (no T1w involved)
    if not rpool.check_rpool('space-EPItemplate_desc-brain_bold'):
        if coregistration not in pipeline_blocks:
            pipeline_blocks += [coregistration_prep_vol, coregistration_prep_mean]
        EPI_reg_blocks = [
            [register_ANTs_EPI_to_template, register_FSL_EPI_to_template]
        ]
        pipeline_blocks += EPI_reg_blocks

    if (cfg['registration_workflows', 'functional_registration',
            'EPI_registration', 'run'] and
        'EPI_Template' in cfg['segmentation', 'tissue_segmentation',
                              'Template_Based', 'template_for_segmentation']):
        pipeline_blocks.append(warp_tissuemask_to_EPItemplate)

    # Generate the composite transform for BOLD-to-template for the T1
    # anatomical template (the BOLD-to- EPI template is already created above)
    if cfg.registration_workflows['functional_registration'][
        'coregistration']['run'
    ] and 'T1_template' in cfg.registration_workflows[
        'functional_registration']['func_registration_to_template'][
            'target_template']['using']:
        pipeline_blocks += [create_func_to_T1template_xfm]

        if cfg.voxel_mirrored_homotopic_connectivity['run']:
            pipeline_blocks += [create_func_to_T1template_symmetric_xfm]

    # Nuisance Correction
    generate_only = True not in cfg['nuisance_corrections',
                                    '2-nuisance_regression', 'run']
    if not rpool.check_rpool('desc-cleaned_bold'):
        nuisance = [ICA_AROMA_ANTsreg, ICA_AROMA_FSLreg,
                    ICA_AROMA_ANTsEPIreg, ICA_AROMA_FSLEPIreg]

        nuisance_masks = [erode_mask_T1w,
                          erode_mask_CSF,
                          erode_mask_GM,
                          erode_mask_WM,
                          erode_mask_bold,
                          erode_mask_boldCSF,
                          erode_mask_boldGM,
                          erode_mask_boldWM]
        nuisance += nuisance_masks + choose_nuisance_blocks(cfg, rpool, \
                                    generate_only)

        pipeline_blocks += nuisance

    pipeline_blocks.append(ingress_regressors)

    apply_func_warp = {}
    _r_w_f_r = cfg.registration_workflows['functional_registration']
    # Warp the functional time series to template space
    apply_func_warp['T1'] = (
        _r_w_f_r['coregistration']['run'] and
        _r_w_f_r['func_registration_to_template']['run'])
    template_funcs = [
        'space-template_desc-preproc_bold',
        'space-template_bold'
    ]
    for func in template_funcs:
        if rpool.check_rpool(func):
            apply_func_warp['T1'] = False

    target_space_nuis = cfg.nuisance_corrections['2-nuisance_regression'][
        'space']
    target_space_alff = cfg.amplitude_low_frequency_fluctuation['target_space']
    target_space_reho = cfg.regional_homogeneity['target_space']

    if apply_func_warp['T1']:

        ts_to_T1template_block = [apply_phasediff_to_timeseries_separately,
                                  apply_blip_to_timeseries_separately,
                                  warp_timeseries_to_T1template,
                                  warp_timeseries_to_T1template_dcan_nhp]

        if 'Template' in target_space_alff or 'Template' in target_space_reho:
            ts_to_T1template_block += [warp_timeseries_to_T1template_deriv]

        if cfg.nuisance_corrections['2-nuisance_regression']['create_regressors']:
            ts_to_T1template_block += [(warp_timeseries_to_T1template_abcd, ('desc-preproc_bold', 'bold'))]
            ts_to_T1template_block.append(single_step_resample_timeseries_to_T1template)
        else:
            ts_to_T1template_block.append(warp_timeseries_to_T1template_abcd)
            ts_to_T1template_block.append(single_step_resample_timeseries_to_T1template)

        pipeline_blocks += [ts_to_T1template_block,
                            warp_sbref_to_T1template]

    if not rpool.check_rpool('space-template_desc-bold_mask'):
        pipeline_blocks += [warp_bold_mask_to_T1template,
                            warp_deriv_mask_to_T1template]

    pipeline_blocks += [func_despike_template]

    if 'Template' in target_space_alff and target_space_nuis == 'native':
        pipeline_blocks += [warp_denoiseNofilt_to_T1template]

    template = cfg.registration_workflows['functional_registration'][
        'func_registration_to_template']['target_template']['using']

    if 'T1_template' in template:
        apply_func_warp['EPI'] = (_r_w_f_r['coregistration']['run'] and
                                  _r_w_f_r['func_registration_to_template'
                                           ]['run_EPI'])
    else:
        apply_func_warp['EPI'] = (_r_w_f_r['func_registration_to_template'
                                           ]['run_EPI'])

    del _r_w_f_r

    template_funcs = [
        'space-EPItemplate_desc-cleaned_bold',
        'space-EPItemplate_desc-brain_bold',
        'space-EPItemplate_desc-motion_bold',
        'space-EPItemplate_desc-preproc_bold',
        'space-EPItemplate_bold'
    ]
    for func in template_funcs:
        if rpool.check_rpool(func):
            apply_func_warp['EPI'] = False

    if apply_func_warp['EPI']:
        pipeline_blocks += [warp_timeseries_to_EPItemplate,
                            warp_bold_mean_to_EPItemplate]

    if not rpool.check_rpool('space-EPItemplate_desc-bold_mask'):
        pipeline_blocks += [warp_bold_mask_to_EPItemplate,
                            warp_deriv_mask_to_EPItemplate]

    # Template-space nuisance regression
    nuisance_template = (cfg['nuisance_corrections', '2-nuisance_regression',
                             'space'] == 'template') and (not generate_only)
    if nuisance_template:
        pipeline_blocks += [nuisance_regression_template]
        # pipeline_blocks += [(nuisance_regression_template,
        #                     ("desc-preproc_bold", "desc-stc_bold"))]

    # PostFreeSurfer and fMRISurface
    if not rpool.check_rpool('space-fsLR_den-32k_bold.dtseries'):
        
        pipeline_blocks += [surface_postproc]

    if not rpool.check_rpool('surf_falff'):
        pipeline_blocks += [surface_falff]

    if not rpool.check_rpool('surf_alff'):
        pipeline_blocks += [surface_alff]

    if not rpool.check_rpool('surf-L_reho') or not rpool.check_rpool('surf-R_reho') :
        pipeline_blocks += [surface_reho]

    if not rpool.check_rpool('space-fsLR_den-32k_bold_surf-correlation_matrix'):
        pipeline_blocks += [surface_connectivity_matrix]

    # Extractions and Derivatives
    tse_atlases, sca_atlases = gather_extraction_maps(cfg)
    cfg.timeseries_extraction['tse_atlases'] = tse_atlases
    cfg.seed_based_correlation_analysis['sca_atlases'] = sca_atlases

    if not rpool.check_rpool('space-template_desc-Mean_timeseries') and \
                    'Avg' in tse_atlases:
        pipeline_blocks += [timeseries_extraction_AVG]

    if not rpool.check_rpool('desc-Voxel_timeseries') and \
                    'Voxel' in tse_atlases:
        pipeline_blocks += [timeseries_extraction_Voxel]

    if not rpool.check_rpool('desc-SpatReg_timeseries') and \
                    'SpatialReg' in tse_atlases:
        pipeline_blocks += [spatial_regression]

    if not rpool.check_rpool('space-template_desc-MeanSCA_correlations') and \
                    'Avg' in sca_atlases:
        pipeline_blocks += [SCA_AVG]

    if not rpool.check_rpool('space-template_desc-DualReg_correlations') and \
                    'DualReg' in sca_atlases:
        pipeline_blocks += [dual_regression]

    if not rpool.check_rpool('space-template_desc-MultReg_correlations') and \
                    'MultReg' in sca_atlases:
        pipeline_blocks += [multiple_regression]

    if 'Native' in target_space_alff:
        if not rpool.check_rpool('alff'):
            pipeline_blocks += [alff_falff]

    if 'Template' in target_space_alff:
        if not rpool.check_rpool('space-template_alff'):
            pipeline_blocks += [alff_falff_space_template]

    if 'Native' in target_space_reho:
        if not rpool.check_rpool('reho'):
            pipeline_blocks += [reho]

    if 'Template' in target_space_reho:
        if not rpool.check_rpool('space-template_reho'):
            pipeline_blocks += [reho_space_template]

    if not rpool.check_rpool('vmhc'):
        pipeline_blocks += [smooth_func_vmhc,
                            warp_timeseries_to_sym_template,
                            vmhc]

    if not rpool.check_rpool('centrality') and \
            any(cfg.network_centrality[option]['weight_options'] for
                option in valid_options['centrality']['method_options']):
        pipeline_blocks += [network_centrality]

    if cfg.pipeline_setup['output_directory']['quality_control'][
        'generate_xcpqc_files'
    ]:
        pipeline_blocks += [qc_xcp]

    if cfg.pipeline_setup['output_directory']['quality_control'][
        'generate_quality_control_images'
    ]:
        qc_stack, qc_montage_id_a, qc_montage_id_s, qc_hist_id, qc_plot_id = \
            create_qc_workflow(cfg)
        pipeline_blocks += qc_stack

    # Connect the entire pipeline!
    try:
        
        wf = connect_pipeline(wf, cfg, rpool, pipeline_blocks)
        
    except LookupError as lookup_error:
        missing_key = None
        errorstrings = [arg for arg in lookup_error.args[0].split('\n') if
                        arg.strip()]
        if lookup_error.args[0].startswith('When trying to connect node b'):
            missing_key = lookup_error.args[0].split("': ")[-1]
        for errorstring in [
            '[!] C-PAC says: The listed resource is not in the resource pool:',
            '[!] C-PAC says: None of the listed resources are in the resource '
            'pool:',
            '[!] C-PAC says: None of the listed resources in the node block '
            'being connected exist in the resource pool.\n\nResources:'
        ]:
            if errorstring in lookup_error.args[0]:
                missing_key = errorstrings[errorstrings.index(errorstring) + 1]
        if missing_key and missing_key.endswith('_bold'
                                                ) and 'func' not in sub_dict:
            raise FileNotFoundError(
                'The provided pipeline configuration requires functional '
                'data but no functional data were found for ' +
                '/'.join([sub_dict[key] for key in ['site', 'subject_id',
                         'unique_id'] if key in sub_dict]) + '. Please check '
                'your data and pipeline configurations.') from lookup_error
        raise lookup_error

    # Write out the data
    # TODO enforce value with schema validation
    try:
        encrypt_data = bool(cfg.pipeline_setup['Amazon-AWS']['s3_encryption'])
    except:
        encrypt_data = False

    # TODO enforce value with schema validation
    # Extract credentials path for output if it exists
    try:
        # Get path to creds file
        creds_path = ''
        if cfg.pipeline_setup['Amazon-AWS']['aws_output_bucket_credentials']:
            creds_path = str(cfg.pipeline_setup['Amazon-AWS'][
                                 'aws_output_bucket_credentials'])
            creds_path = os.path.abspath(creds_path)

        if cfg.pipeline_setup['output_directory'][
            'path'].lower().startswith('s3://'):
            # Test for s3 write access
            s3_write_access = \
                aws_utils.test_bucket_access(creds_path,
                                             cfg.pipeline_setup[
                                                 'output_directory']['path'])

            if not s3_write_access:
                raise Exception('Not able to write to bucket!')

    except Exception as e:
        if cfg.pipeline_setup['output_directory'][
            'path'].lower().startswith('s3://'):
            err_msg = 'There was an error processing credentials or ' \
                      'accessing the S3 bucket. Check and try again.\n' \
                      'Error: %s' % e
            raise Exception(err_msg)

    # Collect all pipeline variants and write to output directory
    rpool.gather_pipes(wf, cfg)

    return wf
