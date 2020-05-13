import os
import time
import six
import re
import csv
import shutil
import pickle
import copy
import json

import pandas as pd
import pkg_resources as p
import networkx as nx
import logging as cb_logging
from time import strftime

import nipype
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.interfaces.afni as afni
from nipype.interfaces.afni import preprocess
import nipype.interfaces.ants as ants
import nipype.interfaces.c3 as c3
from nipype.interfaces.utility import Merge
from nipype.pipeline.engine.utils import format_dot
from nipype import config
from nipype import logging

from indi_aws import aws_utils, fetch_creds

import CPAC
from CPAC.network_centrality.pipeline import (
    create_network_centrality_workflow
)
from CPAC.anat_preproc.anat_preproc import create_anat_preproc
from CPAC.anat_preproc.lesion_preproc import create_lesion_preproc
from CPAC.distortion_correction.distortion_correction import (
    create_EPI_DistCorr,
    blip_distcor_wf
)
from CPAC.func_preproc.func_preproc import (
    create_func_preproc,
    connect_func_preproc,
    slice_timing_wf,
    create_wf_edit_func,
    create_scale_func_wf
)
from CPAC.seg_preproc.seg_preproc import (
    create_seg_preproc,
    create_seg_preproc_template_based
)

from CPAC.seg_preproc.utils import mask_erosion

from CPAC.image_utils import (
    spatial_smooth_outputs,
    z_score_standardize,
    fisher_z_score_standardize,
    calc_avg
)

from CPAC.registration import (
    create_fsl_flirt_linear_reg,
    create_fsl_fnirt_nonlinear_reg,
    create_register_func_to_anat,
    create_bbregister_func_to_anat,
    create_register_func_to_epi,
    create_wf_calculate_ants_warp,
    output_func_to_standard
)

from CPAC.nuisance import create_regressor_workflow, \
    create_nuisance_regression_workflow, \
    filtering_bold_and_regressors, \
    bandpass_voxels, \
    NuisanceRegressor
from CPAC.aroma import create_aroma
from CPAC.median_angle import create_median_angle_correction
from CPAC.generate_motion_statistics import motion_power_statistics
from CPAC.scrubbing import create_scrubbing_preproc
from CPAC.timeseries import (
    get_roi_timeseries,
    get_voxel_timeseries,
    get_vertices_timeseries,
    get_spatial_map_timeseries
)

from CPAC.vmhc.vmhc import create_vmhc
from CPAC.reho.reho import create_reho
from CPAC.alff.alff import create_alff
from CPAC.sca.sca import create_sca, create_temporal_reg

from CPAC.connectome.pipeline import create_connectome

from CPAC.utils.datasource import (
    create_func_datasource,
    create_anat_datasource,
    create_fmap_datasource,
    get_fmap_phasediff_metadata,
    calc_deltaTE_and_asym_ratio,
    create_roi_mask_dataflow,
    create_spatial_map_dataflow,
    create_check_for_s3_node,
    resolve_resolution,
    resample_func_roi,
    match_epi_fmaps
)
from CPAC.utils.trimmer import the_trimmer
from CPAC.utils import Configuration, Strategy, Outputs, find_files
from CPAC.utils.interfaces.function import Function

from CPAC.utils.interfaces.datasink import DataSink

from CPAC.qc.pipeline import create_qc_workflow
from CPAC.qc.utils import generate_qc_pages

from CPAC.utils.utils import (
    extract_one_d,
    get_scan_params,
    get_tr,
    extract_txt,
    extract_output_mean,
    create_output_mean_csv,
    get_zscore,
    get_fisher_zscore,
    pick_wm
)

from CPAC.utils.monitoring import log_nodes_initial, log_nodes_cb

logger = logging.getLogger('nipype.workflow')
# config.enable_debug_mode()

def run_workflow(sub_dict, c, run, pipeline_timing_info=None,
                 p_name=None, plugin='MultiProc', plugin_args=None, test_config=False):
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
    plugin : string (optional); defaule='MultiProc'
        nipype plugin to utilize when the workflow is ran
    plugin_args : dictionary (optional); default=None
        plugin-specific arguments for the workflow plugin

    Returns
    -------
    workflow : nipype workflow
        the prepared nipype workflow object containing the parameters
        specified in the config
    '''

    # Import packages
    from CPAC.utils.utils import check_config_resources, check_system_deps

    # Assure that changes on config will not affect other parts
    c = copy.copy(c)

    subject_id = sub_dict['subject_id']
    if sub_dict['unique_id']:
        subject_id += "_" + sub_dict['unique_id']

    log_dir = os.path.join(c.logDirectory, 'pipeline_%s' % c.pipelineName, subject_id)
    if not os.path.exists(log_dir):
        os.makedirs(os.path.join(log_dir))

    # TODO ASH Enforce c.run_logging to be boolean
    # TODO ASH Schema validation
    config.update_config({
        'logging': {
            'log_directory': log_dir,
            'log_to_file': bool(getattr(c, 'run_logging', True))
        }
    })
    config.enable_resource_monitor()

    logging.update_logging(config)

    # Start timing here
    pipeline_start_time = time.time()
    # at end of workflow, take timestamp again, take time elapsed and check
    # tempfile add time to time data structure inside tempfile, and increment
    # number of subjects

    # Check pipeline config resources
    sub_mem_gb, num_cores_per_sub, num_ants_cores = \
        check_config_resources(c)

    if not plugin:
        plugin = 'MultiProc'

    if plugin_args:
        plugin_args['memory_gb'] = sub_mem_gb
        plugin_args['n_procs'] = num_cores_per_sub
    else:
        plugin_args = {'memory_gb': sub_mem_gb, 'n_procs': num_cores_per_sub}

    # perhaps in future allow user to set threads maximum
    # this is for centrality mostly
    # import mkl
    numThreads = '1'
    os.environ['OMP_NUM_THREADS'] = '1'  # str(num_cores_per_sub)
    os.environ['MKL_NUM_THREADS'] = '1'  # str(num_cores_per_sub)
    os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = str(num_ants_cores)

    # TODO: TEMPORARY
    # TODO: solve the UNet model hanging issue during MultiProc
    if "unet" in c.skullstrip_option:
        c.maxCoresPerParticipant = 1
        logger.info("\n\n[!] LOCKING CPUs PER PARTICIPANT TO 1 FOR U-NET "
                    "MODEL.\n\nThis is a temporary measure due to a known "
                    "issue preventing Nipype's parallelization from running "
                    "U-Net properly.\n\n")

    # calculate maximum potential use of cores according to current pipeline
    # configuration
    max_core_usage = int(c.maxCoresPerParticipant) * \
        int(c.numParticipantsAtOnce)

    ndmg_out = False
    try:
        if "ndmg" in c.output_tree:
            ndmg_out = True
    except:
        pass

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
        encrypt_data = bool(c.s3Encryption[0])
    except:
        encrypt_data = False


    information = """

    C-PAC version: {cpac_version}

    Setting maximum number of cores per participant to {cores}
    Setting number of participants at once to {participants}
    Setting OMP_NUM_THREADS to {threads}
    Setting MKL_NUM_THREADS to {threads}
    Setting ANTS/ITK thread usage to {ants_threads}
    Maximum potential number of cores that might be used during this run: {max_cores}

"""

    execution_info = """

    End of subject workflow {workflow}

    CPAC run complete:

        Pipeline configuration: {pipeline}
        Subject workflow: {workflow}
        Elapsed run time (minutes): {elapsed}
        Timing information saved in {log_dir}/cpac_individual_timing_{pipeline}.csv
        System time of start:      {run_start}
        System time of completion: {run_finish}

"""

    logger.info(information.format(
        cpac_version=CPAC.__version__,
        cores=c.maxCoresPerParticipant,
        participants=c.numParticipantsAtOnce,
        threads=numThreads,
        ants_threads=c.num_ants_threads,
        max_cores=max_core_usage
    ))


    subject_info = {}
    subject_info['subject_id'] = subject_id
    subject_info['start_time'] = pipeline_start_time

    check_centrality_degree = 1 in c.runNetworkCentrality and \
                              (True in c.degWeightOptions or \
                               True in c.eigWeightOptions)

    check_centrality_lfcd = 1 in c.runNetworkCentrality and \
                            True in c.lfcdWeightOptions

    # Check system dependencies
    check_system_deps(check_ants='ANTS' in c.regOption,
                      check_ica_aroma='1' in str(c.runICA[0]),
                      check_centrality_degree=check_centrality_degree,
                      check_centrality_lfcd=check_centrality_lfcd)


    # absolute paths of the dirs
    c.workingDirectory = os.path.abspath(c.workingDirectory)
    if 's3://' not in c.outputDirectory:
        c.outputDirectory = os.path.abspath(c.outputDirectory)

    workflow, strat_list, pipeline_ids = build_workflow(
        subject_id, sub_dict, c, p_name, num_ants_cores
    )

    forks = "\n\nStrategy forks:\n" + \
            "\n".join(["- " + pipe for pipe in sorted(set(pipeline_ids))]) + \
            "\n\n"

    logger.info(forks)

    if test_config:
        logger.info('This has been a test of the pipeline configuration '
                    'file, the pipeline was built successfully, but was '
                    'not run')
    else:

        working_dir = os.path.join(c.workingDirectory, workflow.name)

        #if c.write_debugging_outputs:
        #    with open(os.path.join(working_dir, 'resource_pool.pkl'), 'wb') as f:
        #        pickle.dump(strat_list, f)

        if c.reGenerateOutputs is True:

            erasable = list(find_files(working_dir, '*sink*')) + \
                list(find_files(working_dir, '*link*')) + \
                list(find_files(working_dir, '*log*'))

            for f in erasable:
                if os.path.isfile(f):
                    os.remove(f)
                else:
                    shutil.rmtree(f)

        if hasattr(c, 'trim') and c.trim:

            logger.warn("""
Trimming is an experimental feature, and if used wrongly, it can lead to unreproducible results.
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
                output_dir=c.outputDirectory,
                s3_creds_path=input_creds_path,
            )

        pipeline_start_datetime = strftime("%Y-%m-%d %H:%M:%S")

        try:
            subject_info['resource_pool'] = []

            for strat_no, strat in enumerate(strat_list):
                strat_label = 'strat_%d' % strat_no
                subject_info[strat_label] = strat.get_name()
                subject_info['resource_pool'].append(strat.get_resource_pool())

            subject_info['status'] = 'Running'

            # Create callback logger
            cb_log_filename = os.path.join(log_dir,
                                        'callback.log')

            try:
                if not os.path.exists(os.path.dirname(cb_log_filename)):
                    os.makedirs(os.path.dirname(cb_log_filename))
            except IOError:
                pass

            # Add handler to callback log file
            cb_logger = cb_logging.getLogger('callback')
            cb_logger.setLevel(cb_logging.DEBUG)
            handler = cb_logging.FileHandler(cb_log_filename)
            cb_logger.addHandler(handler)

            # Log initial information from all the nodes
            log_nodes_initial(workflow)

            # Add status callback function that writes in callback log
            if nipype.__version__ not in ('1.1.2'):
                err_msg = "This version of Nipype may not be compatible with " \
                            "CPAC v%s, please install Nipype version 1.1.2\n" \
                            % (CPAC.__version__)
                logger.error(err_msg)
            else:
                plugin_args['status_callback'] = log_nodes_cb

            if plugin_args['n_procs'] == 1:
                plugin = 'Linear'

            try:
                # Actually run the pipeline now, for the current subject
                workflow.run(plugin=plugin, plugin_args=plugin_args)
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
            if 1 in c.run_pypeer:
                from CPAC.pypeer.peer import prep_for_pypeer
                prep_for_pypeer(c.peer_eye_scan_names, c.peer_data_scan_names,
                                c.eye_mask_path, c.outputDirectory, subject_id,
                                pipeline_ids, c.peer_stimulus_path, c.peer_gsr,
                                c.peer_scrub, c.peer_scrub_thresh)

            # Dump subject info pickle file to subject log dir
            subject_info['status'] = 'Completed'

            subject_info_file = os.path.join(
                log_dir, 'subject_info_%s.pkl' % subject_id
            )
            with open(subject_info_file, 'wb') as info:
                pickle.dump(list(subject_info), info)

            # have this check in case the user runs cpac_runner from terminal and
            # the timing parameter list is not supplied as usual by the GUI
            if pipeline_timing_info != None:

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

                timing_temp_file_path = os.path.join(c.logDirectory,
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
                    pipelineTimeDict['Pipeline'] = c.pipelineName
                    pipelineTimeDict['Cores_Per_Subject'] = c.maxCoresPerParticipant
                    pipelineTimeDict['Simultaneous_Subjects'] = c.numParticipantsAtOnce
                    pipelineTimeDict['Number_of_Subjects'] = num_subjects
                    pipelineTimeDict['Start_Time'] = pipeline_start_stamp
                    pipelineTimeDict['End_Time'] = strftime("%Y-%m-%d_%H:%M:%S")
                    pipelineTimeDict['Elapsed_Time_(minutes)'] = elapsedTimeBin[0]
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
                        c.logDirectory,
                        'cpac_individual_timing_%s.csv' % c.pipelineName
                    ), 'a') as timeCSV, open(os.path.join(
                        c.logDirectory,
                        'cpac_individual_timing_%s.csv' % c.pipelineName
                    ), 'rb') as readTimeCSV:

                        timeWriter = csv.DictWriter(timeCSV, fieldnames=gpaTimeFields)
                        timeReader = csv.DictReader(readTimeCSV)

                        headerExists = False
                        for line in timeReader:
                            if 'Start_Time' in line:
                                headerExists = True

                        if headerExists == False:
                            timeWriter.writerow(timeHeader)

                        timeWriter.writerow(pipelineTimeDict)

                    # remove the temp timing file now that it is no longer needed
                    os.remove(timing_temp_file_path)

            # Upload logs to s3 if s3_str in output directory
            if c.outputDirectory.lower().startswith('s3://'):

                try:
                    # Store logs in s3 output director/logs/...
                    s3_log_dir = os.path.join(
                        c.outputDirectory,
                        'logs',
                        os.path.basename(log_dir)
                    )
                    bucket_name = c.outputDirectory.split('/')[2]
                    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

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
                    err_msg = 'Unable to upload CPAC log files in: %s.\nError: %s'
                    logger.error(err_msg, log_dir, exc)

        except Exception as e:
            import traceback; traceback.print_exc()
            execution_info = """

Error of subject workflow {workflow}

CPAC run error:

    Pipeline configuration: {pipeline}
    Subject workflow: {workflow}
    Elapsed run time (minutes): {elapsed}
    Timing information saved in {log_dir}/cpac_individual_timing_{pipeline}.csv
    System time of start:      {run_start}

"""

        finally:

            if 1 in c.generateQualityControlImages and not ndmg_out:
                for pip_id in pipeline_ids:
                    pipeline_base = os.path.join(c.outputDirectory,
                                                 'pipeline_{0}'.format(pip_id))

                    sub_output_dir = os.path.join(pipeline_base, subject_id)
                    qc_dir = os.path.join(sub_output_dir, 'qc')
                    generate_qc_pages(qc_dir)

            if workflow:

                logger.info(execution_info.format(
                    workflow=workflow.name,
                    pipeline=c.pipelineName,
                    log_dir=c.logDirectory,
                    elapsed=(time.time() - pipeline_start_time) / 60,
                    run_start=pipeline_start_datetime,
                    run_finish=strftime("%Y-%m-%d %H:%M:%S")
                ))

                # Remove working directory when done
                if c.removeWorkingDir:
                    try:
                        if os.path.exists(working_dir):
                            logger.info("Removing working dir: %s", working_dir)
                            shutil.rmtree(working_dir)
                    except:
                        logger.warn('Could not remove working directory %s',
                                    working_dir)


def build_workflow(subject_id, sub_dict, c, pipeline_name=None, num_ants_cores=1):

    # TODO ASH temporary code, remove
    # TODO ASH maybe scheme validation/normalization
    already_skullstripped = c.already_skullstripped[0]
    if already_skullstripped == 2:
        already_skullstripped = 0
    elif already_skullstripped == 3:
        already_skullstripped = 1

    if 'ANTS' in c.regOption:

        # if someone doesn't have anatRegANTSinterpolation in their pipe config,
        # it will default to LanczosWindowedSinc
        if not hasattr(c, 'anatRegANTSinterpolation'):
            setattr(c, 'anatRegANTSinterpolation', 'LanczosWindowedSinc')

        if c.anatRegANTSinterpolation not in ['Linear', 'BSpline', 'LanczosWindowedSinc']:
            err_msg = 'The selected ANTS interpolation method may be in the list of values: "Linear", "BSpline", "LanczosWindowedSinc"'
            raise Exception(err_msg)

        # if someone doesn't have funcRegANTSinterpolation in their pipe config,
        # it will default to LanczosWindowedSinc
        if not hasattr(c, 'funcRegANTSinterpolation'):
               setattr(c, 'funcRegANTSinterpolation', 'LanczosWindowedSinc')

        if c.funcRegANTSinterpolation not in ['Linear', 'BSpline', 'LanczosWindowedSinc']:
            err_msg = 'The selected ANTS interpolation method may be in the list of values: "Linear", "BSpline", "LanczosWindowedSinc"'
            raise Exception(err_msg)

    if 'FSL' in c.regOption:

        # if someone doesn't have anatRegFSLinterpolation in their pipe config,
        # it will default to sinc
        if not hasattr(c, 'anatRegFSLinterpolation'):
            setattr(c, 'anatRegFSLinterpolation', 'sinc')

        if c.anatRegFSLinterpolation not in ["trilinear", "sinc", "spline"]:
            err_msg = 'The selected FSL interpolation method may be in the list of values: "trilinear", "sinc", "spline"'
            raise Exception(err_msg)


    # Workflow setup
    workflow_name = 'resting_preproc_' + str(subject_id)
    workflow = pe.Workflow(name=workflow_name)
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

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

    # check if lateral_ventricles_mask exist
    if str(c.lateral_ventricles_mask).lower() in ['none', 'false']:
        ventricle_mask_exist = False
    else:
        ventricle_mask_exist = True

    # TODO ASH normalize file paths with schema validator
    template_keys = [
        ("anat", "templateSpecificationFile"),
        ("anat", "lateral_ventricles_mask"),
        ("anat", "PRIORS_CSF"),
        ("anat", "PRIORS_GRAY"),
        ("anat", "PRIORS_WHITE"),
        ("other", "configFileTwomm"),
        ("anat", "template_based_segmentation_CSF"),
        ("anat", "template_based_segmentation_GRAY"),
        ("anat", "template_based_segmentation_WHITE"),
    ]

    for key_type, key in template_keys:

        node = create_check_for_s3_node(
            key,
            getattr(c, key), key_type,
            input_creds_path, c.workingDirectory
        )

        setattr(c, key, node)

    """""""""""""""""""""""""""""""""""""""""""""""""""
     PREPROCESSING
    """""""""""""""""""""""""""""""""""""""""""""""""""

    strat_initial = Strategy()

    # The list of strategies that will be shared all along the pipeline creation
    strat_list = []
    num_strat = 0

    anat_flow = create_anat_datasource('anat_gather_%d' % num_strat)
    anat_flow.inputs.inputnode.subject = subject_id
    anat_flow.inputs.inputnode.anat = sub_dict['anat']
    anat_flow.inputs.inputnode.creds_path = input_creds_path
    anat_flow.inputs.inputnode.dl_dir = c.workingDirectory

    strat_initial.update_resource_pool({
        'anatomical': (anat_flow, 'outputspec.anat')
    })

    if 'brain_mask' in sub_dict.keys():
        if sub_dict['brain_mask'] and sub_dict['brain_mask'].lower() != 'none':
            brain_flow = create_anat_datasource('brain_gather_%d' % num_strat)
            brain_flow.inputs.inputnode.subject = subject_id
            brain_flow.inputs.inputnode.anat = sub_dict['brain_mask']
            brain_flow.inputs.inputnode.creds_path = input_creds_path
            brain_flow.inputs.inputnode.dl_dir = c.workingDirectory

            strat_initial.update_resource_pool({
                'anatomical_brain_mask': (brain_flow, 'outputspec.anat')
            })

    if 'lesion_mask' in sub_dict.keys():
        lesion_datasource = create_anat_datasource('lesion_gather_%d' % num_strat)
        lesion_datasource.inputs.inputnode.subject = subject_id
        lesion_datasource.inputs.inputnode.anat = sub_dict['lesion_mask']
        lesion_datasource.inputs.inputnode.creds_path = input_creds_path
        lesion_datasource.inputs.inputnode.dl_dir = c.workingDirectory

        strat_initial.update_resource_pool({
            'lesion_mask': (lesion_datasource, 'outputspec.anat')
        })

    templates_for_resampling = [
        (c.resolution_for_anat, c.template_brain_only_for_anat, 'template_brain_for_anat', 'resolution_for_anat'),
        (c.resolution_for_anat, c.template_skull_for_anat, 'template_skull_for_anat', 'resolution_for_anat'),
        (c.resolution_for_anat, c.template_symmetric_brain_only, 'template_symmetric_brain', 'resolution_for_anat'),
        (c.resolution_for_anat, c.template_symmetric_skull, 'template_symmetric_skull', 'resolution_for_anat'),
        (c.resolution_for_anat, c.dilated_symmetric_brain_mask, 'template_dilated_symmetric_brain_mask', 'resolution_for_anat'),
        (c.resolution_for_anat, c.ref_mask, 'template_ref_mask', 'resolution_for_anat'),
        (c.resolution_for_func_preproc, c.template_brain_only_for_func, 'template_brain_for_func_preproc', 'resolution_for_func_preproc'),
        (c.resolution_for_func_preproc, c.template_skull_for_func, 'template_skull_for_func_preproc', 'resolution_for_func_preproc'),
        (c.resolution_for_func_preproc, c.template_epi, 'template_epi', 'resolution_for_func_preproc'),  # no difference of skull and only brain
        (c.resolution_for_func_derivative, c.template_epi, 'template_epi_derivative', 'resolution_for_func_derivative'),  # no difference of skull and only brain
        (c.resolution_for_func_derivative, c.template_brain_only_for_func, 'template_brain_for_func_derivative', 'resolution_for_func_preproc'),
        (c.resolution_for_func_derivative, c.template_skull_for_func, 'template_skull_for_func_derivative', 'resolution_for_func_preproc'),
    ]

    if 1 in c.run_pypeer:
        templates_for_resampling.append((c.resolution_for_func_preproc, c.eye_mask_path, 'template_eye_mask', 'resolution_for_func_preproc'))
        Outputs.any.append("template_eye_mask")

    # update resampled template to resource pool
    for resolution, template, template_name, tag in templates_for_resampling:
        resampled_template = pe.Node(Function(input_names=['resolution', 'template', 'template_name', 'tag'],
                                              output_names=['resampled_template'],
                                              function=resolve_resolution,
                                              as_module=True),
                                        name='resampled_' + template_name)

        resampled_template.inputs.resolution = resolution
        resampled_template.inputs.template = template
        resampled_template.inputs.template_name = template_name
        resampled_template.inputs.tag = tag

        strat_initial.update_resource_pool({
            template_name: (resampled_template, 'resampled_template')
        })

    strat_list += [strat_initial]

    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        if 'anatomical_brain_mask' in strat:

            anat_preproc = create_anat_preproc(method='mask',
                                               config=c,
                                               wf_name='anat_preproc_mask_%d' % num_strat)

            new_strat = strat.fork()
            node, out_file = new_strat['anatomical']
            workflow.connect(node, out_file,
                             anat_preproc, 'inputspec.anat')

            node, out_file = strat['anatomical_brain_mask']
            workflow.connect(node, out_file,
                             anat_preproc, 'inputspec.brain_mask')

            new_strat.append_name(anat_preproc.name)
            new_strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
            new_strat.update_resource_pool({
                'anatomical_brain': (anat_preproc, 'outputspec.brain'),
                'anatomical_reorient': (anat_preproc, 'outputspec.reorient'),
            })
            new_strat.update_resource_pool({
                'anatomical_brain_mask': (anat_preproc, 'outputspec.brain_mask')
            }, override=True)

            new_strat_list += [new_strat]

            continue

        if already_skullstripped:

            anat_preproc = create_anat_preproc(method=None,
                                               already_skullstripped=True,
                                               config=c,
                                               wf_name='anat_preproc_already_%d' % num_strat)

            new_strat = strat.fork()
            node, out_file = new_strat['anatomical']
            workflow.connect(node, out_file,
                            anat_preproc, 'inputspec.anat')
            new_strat.append_name(anat_preproc.name)
            new_strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
            new_strat.update_resource_pool({
                'anatomical_brain': (anat_preproc, 'outputspec.brain'),
                'anatomical_reorient': (anat_preproc, 'outputspec.reorient'),
                'anatomical_brain_mask': (anat_preproc, 'outputspec.brain_mask'),
            })

            new_strat_list += [new_strat]

        else:
            if not any(o in c.skullstrip_option for o in ["AFNI", "FSL", "niworkflows-ants", "unet"]):
                err = '\n\n[!] C-PAC says: Your skull-stripping method options ' \
                    'setting does not include either \'AFNI\' or \'FSL\' or \'niworkflows-ants\'.\n\n' \
                    'Options you provided:\nskullstrip_option: {0}' \
                    '\n\n'.format(str(c.skullstrip_option))
                raise Exception(err)

            if "AFNI" in c.skullstrip_option:

                anat_preproc = create_anat_preproc(method='afni',
                                                   config=c,
                                                   wf_name='anat_preproc_afni_%d' % num_strat)

                anat_preproc.inputs.AFNI_options.set(
                    mask_vol=c.skullstrip_mask_vol,
                    shrink_factor=c.skullstrip_shrink_factor,
                    var_shrink_fac=c.skullstrip_var_shrink_fac,
                    shrink_fac_bot_lim=c.skullstrip_shrink_factor_bot_lim,
                    avoid_vent=c.skullstrip_avoid_vent,
                    niter=c.skullstrip_n_iterations,
                    pushout=c.skullstrip_pushout,
                    touchup=c.skullstrip_touchup,
                    fill_hole=c.skullstrip_fill_hole,
                    avoid_eyes=c.skullstrip_avoid_eyes,
                    use_edge=c.skullstrip_use_edge,
                    exp_frac=c.skullstrip_exp_frac,
                    smooth_final=c.skullstrip_smooth_final,
                    push_to_edge=c.skullstrip_push_to_edge,
                    use_skull=c.skullstrip_use_skull,
                    perc_int=c.skullstrip_perc_int,
                    max_inter_iter=c.skullstrip_max_inter_iter,
                    blur_fwhm=c.skullstrip_blur_fwhm,
                    fac=c.skullstrip_fac,
                    monkey=c.skullstrip_monkey,
                )

                new_strat = strat.fork()
                node, out_file = new_strat['anatomical']
                workflow.connect(node, out_file,
                                 anat_preproc, 'inputspec.anat')
                new_strat.append_name(anat_preproc.name)
                new_strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
                new_strat.update_resource_pool({
                    'anatomical_brain': (anat_preproc, 'outputspec.brain'),
                    'anatomical_reorient': (anat_preproc, 'outputspec.reorient'),
                    'anatomical_brain_mask': (anat_preproc, 'outputspec.brain_mask'),
                })

                new_strat_list += [new_strat]

            if "FSL" in c.skullstrip_option:
                anat_preproc = create_anat_preproc(method='fsl',
                                                   config=c,
                                                   wf_name='anat_preproc_bet_%d' % num_strat)

                anat_preproc.inputs.BET_options.set(
                    frac=c.bet_frac,
                    mask_boolean=c.bet_mask_boolean,
                    mesh_boolean=c.bet_mesh_boolean,
                    outline=c.bet_outline,
                    padding=c.bet_padding,
                    radius=c.bet_radius,
                    reduce_bias=c.bet_reduce_bias,
                    remove_eyes=c.bet_remove_eyes,
                    robust=c.bet_robust,
                    skull=c.bet_skull,
                    surfaces=c.bet_surfaces,
                    threshold=c.bet_threshold,
                    vertical_gradient=c.bet_vertical_gradient,
                )

                new_strat = strat.fork()
                node, out_file = new_strat['anatomical']
                workflow.connect(node, out_file,
                                anat_preproc, 'inputspec.anat')
                new_strat.append_name(anat_preproc.name)
                new_strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
                new_strat.update_resource_pool({
                    'anatomical_brain': (anat_preproc, 'outputspec.brain'),
                    'anatomical_reorient': (anat_preproc, 'outputspec.reorient'),
                    'anatomical_brain_mask': (anat_preproc, 'outputspec.brain_mask'),
                })

                new_strat_list += [new_strat]

            if "niworkflows-ants" in c.skullstrip_option:
                anat_preproc = create_anat_preproc(method='niworkflows-ants',
                                                   config=c,
                                                   wf_name='anat_preproc_niworkflows_ants_%d' % num_strat)

                new_strat = strat.fork()
                node, out_file = new_strat['anatomical']
                workflow.connect(node, out_file,
                                anat_preproc, 'inputspec.anat')
                new_strat.append_name(anat_preproc.name)
                new_strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
                new_strat.update_resource_pool({
                    'anatomical_brain': (anat_preproc, 'outputspec.brain'),
                    'anatomical_reorient': (anat_preproc, 'outputspec.reorient'),
                    'anatomical_brain_mask': (anat_preproc, 'outputspec.brain_mask'),
                })

                new_strat_list += [new_strat]

            if "unet" in c.skullstrip_option:
                anat_preproc = create_anat_preproc(method='unet',
                                                   config=c,
                                                   wf_name='anat_preproc_unet_%d' % num_strat)

                new_strat = strat.fork()
                node, out_file = new_strat['anatomical']
                workflow.connect(node, out_file,
                                anat_preproc, 'inputspec.anat')
                new_strat.append_name(anat_preproc.name)
                new_strat.set_leaf_properties(anat_preproc, 'outputspec.brain')
                new_strat.update_resource_pool({
                    'anatomical_brain': (anat_preproc, 'outputspec.brain'),
                    'anatomical_reorient': (anat_preproc, 'outputspec.reorient'),
                    'anatomical_brain_mask': (anat_preproc, 'outputspec.brain_mask'),
                })

                new_strat_list += [new_strat]


    strat_list = new_strat_list

    new_strat_list = []

    # either run FSL anatomical-to-MNI registration, or...
    if 'FSL' in c.regOption:
        for num_strat, strat in enumerate(strat_list):

            # this is to prevent the user from running FNIRT if they are
            # providing already-skullstripped inputs. this is because
            # FNIRT requires an input with the skull still on
            if already_skullstripped == 1:
                err_msg = '\n\n[!] CPAC says: FNIRT (for anatomical ' \
                          'registration) will not work properly if you ' \
                          'are providing inputs that have already been ' \
                          'skull-stripped.\n\nEither switch to using ' \
                          'ANTS for registration or provide input ' \
                          'images that have not been already ' \
                          'skull-stripped.\n\n'

                logger.info(err_msg)
                raise Exception

            flirt_reg_anat_mni = create_fsl_flirt_linear_reg(
                'anat_mni_flirt_register_%d' % num_strat
            )

            # Input registration parameters
            flirt_reg_anat_mni.inputs.inputspec.interp = c.anatRegFSLinterpolation

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file,
                             flirt_reg_anat_mni, 'inputspec.input_brain')

            # pass the reference files
            node, out_file = strat['template_brain_for_anat']
            workflow.connect(node, out_file,
                flirt_reg_anat_mni, 'inputspec.reference_brain')

            if 'ANTS' in c.regOption:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(flirt_reg_anat_mni.name)
            strat.set_leaf_properties(flirt_reg_anat_mni,
                                      'outputspec.output_brain')

            strat.update_resource_pool({
                'anatomical_to_mni_linear_xfm': (flirt_reg_anat_mni, 'outputspec.linear_xfm'),
                'mni_to_anatomical_linear_xfm': (flirt_reg_anat_mni, 'outputspec.invlinear_xfm'),
                'anatomical_to_standard': (flirt_reg_anat_mni, 'outputspec.output_brain')
            })

    strat_list += new_strat_list

    new_strat_list = []

    try:
        fsl_linear_reg_only = c.fsl_linear_reg_only
    except AttributeError:
        fsl_linear_reg_only = [0]

    if 'FSL' in c.regOption and 0 in fsl_linear_reg_only:

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            if 'anat_mni_flirt_register' in nodes:

                fnirt_reg_anat_mni = create_fsl_fnirt_nonlinear_reg(
                    'anat_mni_fnirt_register_%d' % num_strat
                )

                node, out_file = strat['anatomical_brain']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.input_brain')

                # pass the reference files
                node, out_file = strat['template_brain_for_anat']
                workflow.connect(node, out_file,
                    fnirt_reg_anat_mni, 'inputspec.reference_brain')

                node, out_file = strat['anatomical_reorient']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.input_skull')

                node, out_file = strat['anatomical_to_mni_linear_xfm']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.linear_aff')

                node, out_file = strat['template_skull_for_anat']
                workflow.connect(node, out_file,
                    fnirt_reg_anat_mni, 'inputspec.reference_skull')

                node, out_file = strat['template_ref_mask']
                workflow.connect(node, out_file,
                    fnirt_reg_anat_mni, 'inputspec.ref_mask')

                # assign the FSL FNIRT config file specified in pipeline
                # config.yml
                fnirt_reg_anat_mni.inputs.inputspec.fnirt_config = c.fnirtConfig

                if 1 in fsl_linear_reg_only:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(fnirt_reg_anat_mni.name)
                strat.set_leaf_properties(fnirt_reg_anat_mni,
                                          'outputspec.output_brain')

                strat.update_resource_pool({
                    'anatomical_to_mni_nonlinear_xfm': (fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm'),
                    'anatomical_to_standard': (fnirt_reg_anat_mni, 'outputspec.output_brain')
                }, override=True)

    strat_list += new_strat_list

    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        nodes = strat.get_nodes_names()

        # or run ANTS anatomical-to-MNI registration instead
        if 'ANTS' in c.regOption and \
            'anat_mni_flirt_register' not in nodes and \
                'anat_mni_fnirt_register' not in nodes:

            ants_reg_anat_mni = \
                create_wf_calculate_ants_warp(
                    'anat_mni_ants_register_%d' % num_strat,
                    num_threads=num_ants_cores,
                    reg_ants_skull = c.regWithSkull
                )

            # calculating the transform with the skullstripped is
            # reported to be better, but it requires very high
            # quality skullstripping. If skullstripping is imprecise
            # registration with skull is preferred

            # TODO ASH assess with schema validator
            if 1 in c.regWithSkull:

                if already_skullstripped == 1:
                    err_msg = '\n\n[!] CPAC says: You selected ' \
                        'to run anatomical registration with ' \
                        'the skull, but you also selected to ' \
                        'use already-skullstripped images as ' \
                        'your inputs. This can be changed ' \
                        'in your pipeline configuration ' \
                        'editor.\n\n'

                    logger.info(err_msg)
                    raise Exception

            # Input ANTs registration parameters
            if c.ANTs_para_T1_registration is None:
                err_msg = '\n\n[!] C-PAC says: \n'\
                    'You have selected \'regOption: [ANTS]\' as your anatomical registration method. \n'\
                            'However, ANTs parameters specified: {0}, is not supported. ' \
                            'Please specify ANTs parameters properly and try again'.format(str(c.ANTs_para_T1_registration))
                raise Exception(err_msg)
            else:
                ants_reg_anat_mni.inputs.inputspec.ants_para = c.ANTs_para_T1_registration

            ants_reg_anat_mni.inputs.inputspec.interp = c.anatRegANTSinterpolation

            # get the skull-stripped anatomical from resource pool
            node, out_file = strat['anatomical_brain']

            # pass the anatomical to the workflow
            workflow.connect(node, out_file,
                                ants_reg_anat_mni,
                                'inputspec.moving_brain')

            # pass the reference file
            node, out_file = strat['template_brain_for_anat']
            workflow.connect(node, out_file,
                ants_reg_anat_mni, 'inputspec.reference_brain')

            # get the reorient skull-on anatomical from resource pool
            node, out_file = strat['anatomical_reorient']

            # pass the anatomical to the workflow
            workflow.connect(node, out_file,
                                ants_reg_anat_mni,
                                'inputspec.moving_skull')

            # pass the reference file
            node, out_file = strat['template_skull_for_anat']
            workflow.connect(
                node, out_file,
                ants_reg_anat_mni, 'inputspec.reference_skull'
                )

            # Test if a lesion mask is found for the anatomical image
            if 'lesion_mask' in sub_dict and c.use_lesion_mask \
                    and 'lesion_preproc' not in nodes:
                # Create lesion preproc node to apply afni Refit and Resample
                lesion_preproc = create_lesion_preproc(
                    wf_name='lesion_preproc_%d' % num_strat
                )
                # Add the name of the node in the strat object
                strat.append_name(lesion_preproc.name)
                # I think I don't need to set this node as leaf but not sure
                # strat.set_leaf_properties(lesion_preproc, 'inputspec.lesion')

                # Add the lesion preprocessed to the resource pool
                strat.update_resource_pool({
                    'lesion_reorient': (lesion_preproc, 'outputspec.reorient')
                })
                # The Refit lesion is not added to the resource pool because
                # it is not used afterward

                # Retieve the lesion mask from the resource pool
                node, out_file = strat['lesion_mask']
                # Set the lesion mask as input of lesion_preproc
                workflow.connect(
                    node, out_file,
                    lesion_preproc, 'inputspec.lesion'
                )

                # Set the output of lesion preproc as parameter of ANTs
                # fixed_image_mask option
                workflow.connect(
                    lesion_preproc, 'outputspec.reorient',
                    ants_reg_anat_mni, 'inputspec.fixed_image_mask'
                )
            else:
                ants_reg_anat_mni.inputs.inputspec.fixed_image_mask = None

            strat.append_name(ants_reg_anat_mni.name)

            strat.set_leaf_properties(ants_reg_anat_mni,
                                      'outputspec.normalized_output_brain')

            strat.update_resource_pool({
                'ants_initial_xfm': (ants_reg_anat_mni, 'outputspec.ants_initial_xfm'),
                'ants_rigid_xfm': (ants_reg_anat_mni, 'outputspec.ants_rigid_xfm'),
                'ants_affine_xfm': (ants_reg_anat_mni, 'outputspec.ants_affine_xfm'),
                'anatomical_to_mni_nonlinear_xfm': (ants_reg_anat_mni, 'outputspec.warp_field'),
                'mni_to_anatomical_nonlinear_xfm': (ants_reg_anat_mni, 'outputspec.inverse_warp_field'),
                'anat_to_mni_ants_composite_xfm': (ants_reg_anat_mni, 'outputspec.composite_transform'),
                'anatomical_to_standard': (ants_reg_anat_mni, 'outputspec.normalized_output_brain')
            })

    strat_list += new_strat_list

    # [SYMMETRIC] T1 -> Symmetric Template, Non-linear registration (FNIRT/ANTS)

    new_strat_list = []

    if 1 in c.runVMHC and 1 in getattr(c, 'runFunctional', [1]):

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            if 'FSL' in c.regOption and \
               'anat_mni_ants_register' not in nodes:

                # this is to prevent the user from running FNIRT if they are
                # providing already-skullstripped inputs. this is because
                # FNIRT requires an input with the skull still on
                # TODO ASH normalize w schema validation to bool
                if already_skullstripped == 1:
                    err_msg = '\n\n[!] CPAC says: FNIRT (for anatomical ' \
                              'registration) will not work properly if you ' \
                              'are providing inputs that have already been ' \
                              'skull-stripped.\n\nEither switch to using ' \
                              'ANTS for registration or provide input ' \
                              'images that have not been already ' \
                              'skull-stripped.\n\n'

                    logger.info(err_msg)
                    raise Exception

                flirt_reg_anat_symm_mni = create_fsl_flirt_linear_reg(
                    'anat_symmetric_mni_flirt_register_%d' % num_strat
                )


                # Input registration parameters
                flirt_reg_anat_symm_mni.inputs.inputspec.interp = c.anatRegFSLinterpolation

                node, out_file = strat['anatomical_brain']
                workflow.connect(node, out_file,
                                 flirt_reg_anat_symm_mni,
                                 'inputspec.input_brain')

                # pass the reference files
                node, out_file = strat['template_symmetric_brain']
                workflow.connect(node, out_file,
                    flirt_reg_anat_symm_mni, 'inputspec.reference_brain')

                # if 'ANTS' in c.regOption:
                #    strat = strat.fork()
                #    new_strat_list.append(strat)

                strat.append_name(flirt_reg_anat_symm_mni.name)
                strat.set_leaf_properties(flirt_reg_anat_symm_mni,
                                          'outputspec.output_brain')

                strat.update_resource_pool({
                    'anatomical_to_symmetric_mni_linear_xfm': (
                    flirt_reg_anat_symm_mni, 'outputspec.linear_xfm'),
                    'symmetric_mni_to_anatomical_linear_xfm': (
                    flirt_reg_anat_symm_mni, 'outputspec.invlinear_xfm'),
                    'symmetric_anatomical_to_standard': (
                    flirt_reg_anat_symm_mni, 'outputspec.output_brain')
                })

        strat_list += new_strat_list

        new_strat_list = []

        try:
            fsl_linear_reg_only = c.fsl_linear_reg_only
        except AttributeError:
            fsl_linear_reg_only = [0]

        if 'FSL' in c.regOption and 0 in fsl_linear_reg_only:

            for num_strat, strat in enumerate(strat_list):

                nodes = strat.get_nodes_names()

                if 'anat_mni_flirt_register' in nodes:
                    fnirt_reg_anat_symm_mni = create_fsl_fnirt_nonlinear_reg(
                        'anat_symmetric_mni_fnirt_register_%d' % num_strat
                    )

                    node, out_file = strat['anatomical_brain']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.input_brain')

                    # pass the reference files
                    node, out_file = strat['template_brain_for_anat']
                    workflow.connect(node, out_file,
                        fnirt_reg_anat_symm_mni, 'inputspec.reference_brain')

                    node, out_file = strat['anatomical_reorient']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.input_skull')

                    node, out_file = strat['anatomical_to_mni_linear_xfm']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.linear_aff')

                    node, out_file = strat['template_symmetric_skull']
                    workflow.connect(node, out_file,
                        fnirt_reg_anat_symm_mni, 'inputspec.reference_skull')

                    node, out_file = strat['template_dilated_symmetric_brain_mask']
                    workflow.connect(node, out_file,
                        fnirt_reg_anat_symm_mni, 'inputspec.ref_mask')

                    strat.append_name(fnirt_reg_anat_symm_mni.name)
                    strat.set_leaf_properties(fnirt_reg_anat_symm_mni,
                                              'outputspec.output_brain')

                    strat.update_resource_pool({
                        'anatomical_to_symmetric_mni_nonlinear_xfm': (
                        fnirt_reg_anat_symm_mni, 'outputspec.nonlinear_xfm'),
                        'symmetric_anatomical_to_standard': (
                        fnirt_reg_anat_symm_mni, 'outputspec.output_brain')
                    }, override=True)

        strat_list += new_strat_list

        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            # or run ANTS anatomical-to-MNI registration instead
            if 'ANTS' in c.regOption and \
               'anat_mni_flirt_register' not in nodes and \
               'anat_mni_fnirt_register' not in nodes and \
               'anat_symmetric_mni_flirt_register' not in nodes and \
               'anat_symmetric_mni_fnirt_register' not in nodes:

                ants_reg_anat_symm_mni = \
                    create_wf_calculate_ants_warp(
                        'anat_symmetric_mni_ants_register_%d' % num_strat,
                        num_threads=num_ants_cores,
                        reg_ants_skull = c.regWithSkull
                    )

                # Input registration parameters
                ants_reg_anat_symm_mni.inputs.inputspec.ants_para = c.ANTs_para_T1_registration
                ants_reg_anat_symm_mni.inputs.inputspec.interp = c.anatRegANTSinterpolation

                # calculating the transform with the skullstripped is
                # reported to be better, but it requires very high
                # quality skullstripping. If skullstripping is imprecise
                # registration with skull is preferred
                if 1 in c.regWithSkull:

                    if already_skullstripped == 1:
                        err_msg = '\n\n[!] CPAC says: You selected ' \
                            'to run anatomical registration with ' \
                            'the skull, but you also selected to ' \
                            'use already-skullstripped images as ' \
                            'your inputs. This can be changed ' \
                            'in your pipeline configuration ' \
                            'editor.\n\n'

                        logger.info(err_msg)
                        raise Exception

                # get the skullstripped anatomical from resource pool
                node, out_file = strat['anatomical_brain']

                # pass the anatomical to the workflow
                workflow.connect(node, out_file,
                                    ants_reg_anat_symm_mni,
                                    'inputspec.moving_brain')

                # pass the reference file
                node, out_file = strat['template_symmetric_brain']
                workflow.connect(node, out_file,
                                ants_reg_anat_symm_mni, 'inputspec.reference_brain')

                # get the reorient skull-on anatomical from resource
                # pool
                node, out_file = strat['anatomical_reorient']

                # pass the anatomical to the workflow
                workflow.connect(node, out_file,
                                    ants_reg_anat_symm_mni,
                                    'inputspec.moving_skull')

                # pass the reference file
                node, out_file = strat['template_symmetric_skull']
                workflow.connect(node, out_file,
                                    ants_reg_anat_symm_mni, 'inputspec.reference_skull')


                if 'lesion_mask' in sub_dict and c.use_lesion_mask\
                        and 'lesion_preproc' not in nodes:
                    # Create lesion preproc node to apply afni Refit & Resample
                    lesion_preproc = create_lesion_preproc(
                        wf_name='lesion_preproc_%d' % num_strat
                    )
                    # Add the name of the node in the strat object
                    strat.append_name(lesion_preproc.name)

                    # I think I don't need to set this node as leaf but not sure
                    # strat.set_leaf_properties(lesion_preproc,
                    # 'inputspec.lesion')

                    # Add the lesion preprocessed to the resource pool
                    strat.update_resource_pool({
                        'lesion_reorient': (
                            lesion_preproc, 'outputspec.reorient')
                    })
                    # The Refit lesion is not added to the resource pool because
                    # it is not used afterward

                    # Retieve the lesion mask from the resource pool
                    node, out_file = strat['lesion_mask']
                    # Set the lesion mask as input of lesion_preproc
                    workflow.connect(
                        node, out_file,
                        lesion_preproc, 'inputspec.lesion'
                    )

                    # Set the output of lesion preproc as parameter of ANTs
                    # fixed_image_mask option
                    workflow.connect(
                        lesion_preproc, 'outputspec.reorient',
                        ants_reg_anat_symm_mni, 'inputspec.fixed_image_mask'
                    )
                else:
                    ants_reg_anat_symm_mni.inputs.inputspec.fixed_image_mask = \
                        None

                strat.append_name(ants_reg_anat_symm_mni.name)
                strat.set_leaf_properties(ants_reg_anat_symm_mni,
                                          'outputspec.normalized_output_brain')

                strat.update_resource_pool({
                    'ants_symmetric_initial_xfm': (ants_reg_anat_symm_mni, 'outputspec.ants_initial_xfm'),
                    'ants_symmetric_rigid_xfm': (ants_reg_anat_symm_mni, 'outputspec.ants_rigid_xfm'),
                    'ants_symmetric_affine_xfm': (ants_reg_anat_symm_mni, 'outputspec.ants_affine_xfm'),
                    'anatomical_to_symmetric_mni_nonlinear_xfm': (ants_reg_anat_symm_mni, 'outputspec.warp_field'),
                    'symmetric_mni_to_anatomical_nonlinear_xfm': (ants_reg_anat_symm_mni, 'outputspec.inverse_warp_field'),
                    'anat_to_symmetric_mni_ants_composite_xfm': (ants_reg_anat_symm_mni, 'outputspec.composite_transform'),
                    'symmetric_anatomical_to_standard': (ants_reg_anat_symm_mni, 'outputspec.normalized_output_brain')
                })

        strat_list += new_strat_list

    # Inserting Segmentation Preprocessing Workflow

    new_strat_list = []

    if 1 in c.runSegmentationPreprocessing:

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            seg_preproc = None

            if not any(o in c.seg_use_threshold for o in ["FSL-FAST Thresholding", "Customized Thresholding"]):
                err = '\n\n[!] C-PAC says: Your segmentation thresholding options ' \
                    'setting does not include either \'FSL-FAST Thresholding\' or \'Customized Thresholding\'.\n\n' \
                    'Options you provided:\nseg_use_threshold: {0}' \
                    '\n\n'.format(str(c.seg_use_threshold))
                raise Exception(err)

            # TODO ASH based on config, instead of nodes?
            if 'anat_mni_fnirt_register' in nodes or 'anat_mni_flirt_register' in nodes:
                use_ants = False
            elif 'anat_mni_ants_register' in nodes:
                use_ants = True

            seg_preproc = create_seg_preproc(use_ants=use_ants,
                                             use_priors=c.seg_use_priors,
                                             use_threshold=c.seg_use_threshold,
                                             csf_use_erosion=c.seg_csf_use_erosion,
                                             wm_use_erosion=c.seg_wm_use_erosion,
                                             gm_use_erosion=c.seg_gm_use_erosion,
                                             wf_name='seg_preproc_{0}'.format(num_strat))

            seg_preproc.inputs.csf_threshold.csf_threshold=c.seg_CSF_threshold_value
            seg_preproc.inputs.wm_threshold.wm_threshold=c.seg_WM_threshold_value
            seg_preproc.inputs.gm_threshold.gm_threshold=c.seg_GM_threshold_value

            seg_preproc.inputs.csf_erosion_prop.csf_erosion_prop=c.csf_erosion_prop
            seg_preproc.inputs.wm_erosion_prop.wm_erosion_prop=c.wm_erosion_prop
            seg_preproc.inputs.gm_erosion_prop.gm_erosion_prop=c.gm_erosion_prop

            seg_preproc.inputs.csf_mask_erosion_mm.csf_mask_erosion_mm=c.csf_mask_erosion_mm
            seg_preproc.inputs.wm_mask_erosion_mm.wm_mask_erosion_mm=c.wm_mask_erosion_mm
            seg_preproc.inputs.gm_mask_erosion_mm.gm_mask_erosion_mm=c.gm_mask_erosion_mm

            seg_preproc.inputs.csf_erosion_mm.csf_erosion_mm=c.csf_erosion_mm
            seg_preproc.inputs.wm_erosion_mm.wm_erosion_mm=c.wm_erosion_mm
            seg_preproc.inputs.gm_erosion_mm.gm_erosion_mm=c.gm_erosion_mm

            workflow.connect(anat_preproc, 'outputspec.brain_mask', seg_preproc, 'inputspec.brain_mask')

            # TODO ASH review
            if seg_preproc is None:
                continue

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file,
                             seg_preproc, 'inputspec.brain')

            if 'anat_mni_fnirt_register' in nodes or 'anat_mni_flirt_register' in nodes:
                node, out_file = strat['mni_to_anatomical_linear_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc,
                                 'inputspec.standard2highres_mat')

            elif 'anat_mni_ants_register' in nodes:
                node, out_file = strat['ants_initial_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc,
                                 'inputspec.standard2highres_init')

                node, out_file = strat['ants_rigid_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc,
                                 'inputspec.standard2highres_rig')

                node, out_file = strat['ants_affine_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc,
                                 'inputspec.standard2highres_mat')


            workflow.connect(c.PRIORS_CSF, 'local_path',
                                seg_preproc, 'inputspec.PRIOR_CSF')

            workflow.connect(c.PRIORS_GRAY, 'local_path',
                                seg_preproc, 'inputspec.PRIOR_GRAY')

            workflow.connect(c.PRIORS_WHITE, 'local_path',
                                seg_preproc, 'inputspec.PRIOR_WHITE')

            # TODO ASH review with forking function
            if 0 in c.runSegmentationPreprocessing:
                strat = strat.fork()
                new_strat_list.append(strat)

            if c.brain_use_erosion:
                ero_imports = ['import scipy.ndimage as nd' , 'import numpy as np', 'import nibabel as nb', 'import os']
                eroded_mask = pe.Node(util.Function(input_names = ['roi_mask', 'skullstrip_mask', 'mask_erosion_mm', 'mask_erosion_prop'],
                                                    output_names = ['output_roi_mask', 'eroded_skullstrip_mask'],
                                                    function = mask_erosion,
                                                    imports = ero_imports),
                                                    name='erode_skullstrip_brain_mask')
                eroded_mask.inputs.mask_erosion_mm = c.brain_mask_erosion_mm
                workflow.connect(anat_preproc, 'outputspec.brain_mask', eroded_mask, 'skullstrip_mask')
                workflow.connect(seg_preproc, 'outputspec.csf_probability_map', eroded_mask, 'roi_mask')

                strat.update_resource_pool({'anatomical_eroded_brain_mask': (eroded_mask, 'eroded_skullstrip_mask')})

            strat.append_name(seg_preproc.name)
            strat.update_resource_pool({
                'anatomical_gm_mask': (seg_preproc, 'outputspec.gm_mask'),
                'anatomical_csf_mask': (seg_preproc, 'outputspec.csf_mask'),
                'anatomical_wm_mask': (seg_preproc, 'outputspec.wm_mask'),
                'seg_probability_maps': (seg_preproc, 'outputspec.probability_maps'),
                'seg_mixeltype': (seg_preproc, 'outputspec.mixeltype'),
                'seg_partial_volume_map': (seg_preproc, 'outputspec.partial_volume_map'),
                'seg_partial_volume_files': (seg_preproc, 'outputspec.partial_volume_files')
            })

    strat_list += new_strat_list

    if 'T1_template' in c.template_based_segmentation:

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            if not any(o in c.template_based_segmentation for o in ['EPI_template', 'T1_template', 'None']):
                err = '\n\n[!] C-PAC says: Your template based segmentation ' \
                    'setting does not include either \'EPI_template\' or \'T1_template\'.\n\n' \
                    'Options you provided:\ntemplate_based_segmentation: {0}' \
                    '\n\n'.format(str(c.template_based_segmentation))
                raise Exception(err)

            # TODO ASH based on config, instead of nodes?
            if 'anat_mni_fnirt_register' in nodes or 'anat_mni_flirt_register' in nodes:
                use_ants = False
            elif 'anat_mni_ants_register' in nodes:
                use_ants = True

            seg_preproc_template_based = create_seg_preproc_template_based(use_ants=use_ants,
                                                             wf_name='seg_preproc_t1_template_{0}'.format(num_strat))

            # TODO ASH review
            if seg_preproc_template_based is None:
                continue

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file,
                             seg_preproc_template_based, 'inputspec.brain')

            if 'anat_mni_fnirt_register' in nodes or 'anat_mni_flirt_register' in nodes:
                node, out_file = strat['mni_to_anatomical_linear_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc_template_based,
                                 'inputspec.standard2highres_mat')

            elif 'anat_mni_ants_register' in nodes:
                node, out_file = strat['ants_initial_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc_template_based,
                                 'inputspec.standard2highres_init')

                node, out_file = strat['ants_rigid_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc_template_based,
                                 'inputspec.standard2highres_rig')

                node, out_file = strat['ants_affine_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc_template_based,
                                 'inputspec.standard2highres_mat')

            workflow.connect(c.template_based_segmentation_CSF, 'local_path',
                                seg_preproc_template_based, 'inputspec.CSF_template')

            workflow.connect(c.template_based_segmentation_GRAY, 'local_path',
                                seg_preproc_template_based, 'inputspec.GRAY_template')

            workflow.connect(c.template_based_segmentation_WHITE, 'local_path',
                                seg_preproc_template_based, 'inputspec.WHITE_template')

            # TODO ASH review with forking function
            if 'None' in c.template_based_segmentation:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(seg_preproc_template_based.name)
            strat.update_resource_pool({
                'anatomical_gm_mask': (seg_preproc_template_based, 'outputspec.gm_mask'),
                'anatomical_csf_mask': (seg_preproc_template_based, 'outputspec.csf_mask'),
                'anatomical_wm_mask': (seg_preproc_template_based, 'outputspec.wm_mask')
            })

    strat_list += new_strat_list

    # Inserting Functional Data workflow
    if ('func' in sub_dict or 'rest' in sub_dict) and \
            1 in getattr(c, 'runFunctional', [1]):
        #  pipeline needs to have explicit [0] to disable functional workflow

        for num_strat, strat in enumerate(strat_list):

            if 'func' in sub_dict:
                func_paths_dict = sub_dict['func']
            else:
                func_paths_dict = sub_dict['rest']

            func_wf = create_func_datasource(func_paths_dict,
                                             'func_gather_%d' % num_strat)
            func_wf.inputs.inputnode.set(
                subject=subject_id,
                creds_path=input_creds_path,
                dl_dir=c.workingDirectory
            )
            func_wf.get_node('inputnode').iterables = \
                ("scan", list(func_paths_dict.keys()))

            # Grab field maps
            diff = False
            blip = False
            if "fmap" in sub_dict:
                fmap_rp_list = []
                fmap_TE_list = []
                for key in sub_dict["fmap"]:
                    gather_fmap = create_fmap_datasource(sub_dict["fmap"],
                                                         "fmap_gather_"
                                                         "{0}".format(key))
                    gather_fmap.inputs.inputnode.set(
                        subject=subject_id,
                        creds_path=input_creds_path,
                        dl_dir=c.workingDirectory
                    )
                    gather_fmap.inputs.inputnode.scan = key
                    strat.update_resource_pool({
                        key: (gather_fmap, 'outputspec.rest'),
                        "{0}_scan_params".format(key): (gather_fmap,
                                                        'outputspec.scan_params')
                    })

                    fmap_rp_list.append(key)

                    if key == "diff_phase" or key == "diff_mag_one" or \
                            key == "diff_mag_two":
                        diff = True

                        get_fmap_metadata_imports = ['import json']
                        get_fmap_metadata = pe.Node(Function(
                            input_names=['data_config_scan_params'],
                            output_names=['echo_time',
                                          'dwell_time',
                                          'pe_direction'],
                            function=get_fmap_phasediff_metadata,
                            imports=get_fmap_metadata_imports),
                            name='{0}_get_metadata_{1}'.format(key,
                                                               num_strat))

                        node, out_file = strat["{}_scan_params".format(key)]
                        workflow.connect(node, out_file, get_fmap_metadata,
                                         'data_config_scan_params')

                        strat.update_resource_pool({
                            "{}_TE".format(key): (get_fmap_metadata,
                                                  'echo_time'),
                            "{}_dwell".format(key): (get_fmap_metadata,
                                                     'dwell_time'),
                            "{}_pedir".format(key): (get_fmap_metadata,
                                                     'pe_direction')
                        })
                        fmap_TE_list.append("{}_TE".format(key))

                    if key == "epi_AP" or key == "epi_PA":
                        blip = True

                if diff:
                    calc_delta_ratio = pe.Node(Function(
                        input_names=['dwell_time',
                                     'echo_time_one',
                                     'echo_time_two',
                                     'echo_time_three'],
                        output_names=['deltaTE',
                                      'dwell_asym_ratio'],
                        function=calc_deltaTE_and_asym_ratio),
                        name='diff_distcor_calc_delta_{}'.format(num_strat))

                    node, out_file = strat['diff_phase_dwell']
                    workflow.connect(node, out_file, calc_delta_ratio,
                                     'dwell_time')

                    node, out_file = strat[fmap_TE_list[0]]
                    workflow.connect(node, out_file, calc_delta_ratio,
                                     'echo_time_one')

                    node, out_file = strat[fmap_TE_list[1]]
                    workflow.connect(node, out_file, calc_delta_ratio,
                                     'echo_time_two')

                    if len(fmap_TE_list) > 2:
                        node, out_file = strat[fmap_TE_list[2]]
                        workflow.connect(node, out_file, calc_delta_ratio,
                                         'echo_time_three')

                    strat.update_resource_pool({
                        'deltaTE': (calc_delta_ratio, 'deltaTE'),
                        'dwell_asym_ratio': (calc_delta_ratio,
                                             'dwell_asym_ratio')
                    })

            # Add in nodes to get parameters from configuration file
            # a node which checks if scan_parameters are present for each scan
            scan_params = \
                pe.Node(Function(
                    input_names=['data_config_scan_params',
                                 'subject_id',
                                 'scan',
                                 'pipeconfig_tr',
                                 'pipeconfig_tpattern',
                                 'pipeconfig_start_indx',
                                 'pipeconfig_stop_indx'],
                    output_names=['tr',
                                  'tpattern',
                                  'ref_slice',
                                  'start_indx',
                                  'stop_indx',
                                  'pe_direction'],
                    function=get_scan_params,
                    as_module=True
                ), name='scan_params_%d' % num_strat)

            if "Selected Functional Volume" in c.func_reg_input:

                get_func_volume = pe.Node(interface=afni.Calc(),
                                          name='get_func_volume_{0}'.format(num_strat))

                get_func_volume.inputs.set(
                    expr='a',
                    single_idx=c.func_reg_input_volume,
                    outputtype='NIFTI_GZ'
                )
                workflow.connect(func_wf, 'outputspec.rest',
                                 get_func_volume, 'in_file_a')

            # wire in the scan parameter workflow
            workflow.connect(func_wf, 'outputspec.scan_params',
                             scan_params, 'data_config_scan_params')

            workflow.connect(func_wf, 'outputspec.subject',
                             scan_params, 'subject_id')

            workflow.connect(func_wf, 'outputspec.scan',
                             scan_params, 'scan')

            # connect in constants
            scan_params.inputs.set(
                pipeconfig_start_indx=c.startIdx,
                pipeconfig_stop_indx=c.stopIdx
            )

            strat.update_resource_pool({
                'raw_functional': (func_wf, 'outputspec.rest'),
                'scan_id': (func_wf, 'outputspec.scan'),
                'tr': (scan_params, 'tr'),
                'tpattern': (scan_params, 'tpattern'),
                'start_idx': (scan_params, 'start_indx'),
                'stop_idx': (scan_params, 'stop_indx'),
                'pe_direction': (scan_params, 'pe_direction'),
            })

            strat.set_leaf_properties(func_wf, 'outputspec.rest')

            if "Selected Functional Volume" in c.func_reg_input:
                strat.update_resource_pool({
                    'selected_func_volume': (get_func_volume, 'out_file')
                })


        # scale func data based on configuration information
        for num_strat, strat in enumerate(strat_list):

            scale_func_wf = create_scale_func_wf(
                runScaling=c.runScaling,
                scaling_factor=c.scaling_factor,
                wf_name="scale_func_%d" % (num_strat)
            )

            # connect the functional data from the leaf node into the wf
            node, out_file = strat.get_leaf_properties()
            workflow.connect(node, out_file,
                             scale_func_wf, 'inputspec.func')

            # replace the leaf node with the output from the recently added
            # workflow
            strat.set_leaf_properties(scale_func_wf, 'outputspec.scaled_func')


        # Truncate scan length based on configuration information
        for num_strat, strat in enumerate(strat_list):

            trunc_wf = create_wf_edit_func(
                wf_name="edit_func_%d" % (num_strat)
            )

            # connect the functional data from the leaf node into the wf
            node, out_file = strat.get_leaf_properties()
            workflow.connect(node, out_file,
                             trunc_wf, 'inputspec.func')

            # connect the other input parameters
            node, node_out = strat['start_idx']
            workflow.connect(node, node_out,
                             trunc_wf, 'inputspec.start_idx')

            node, node_out = strat['stop_idx']
            workflow.connect(node, node_out,
                             trunc_wf, 'inputspec.stop_idx')

            # replace the leaf node with the output from the recently added
            # workflow
            strat.set_leaf_properties(trunc_wf, 'outputspec.edited_func')
            strat.update_resource_pool({
                            'raw_functional_trunc': (trunc_wf, 'outputspec.edited_func'),
                        })

        # Motion Statistics Workflow
        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            if 0 in c.runMotionStatisticsFirst:

                new_strat_list += [strat.fork()]

            if 1 in c.runMotionStatisticsFirst:

                for skullstrip_tool in c.functionalMasking:

                    skullstrip_tool = skullstrip_tool.lower()

                    for motion_correct_ref in c.motion_correction_reference:

                        motion_correct_ref = motion_correct_ref.lower()

                        if " " in motion_correct_ref:
                            motion_correct_ref = motion_correct_ref.replace(" ", "_")

                        for motion_correct_tool in c.motion_correction:

                            motion_correct_tool = motion_correct_tool.lower()

                            new_strat = strat.fork()

                            func_preproc = create_func_preproc(
                                skullstrip_tool=skullstrip_tool,
                                motion_correct_tool=motion_correct_tool,
                                motion_correct_ref=motion_correct_ref,
                                config=c,
                                wf_name='func_preproc_before_stc_{0}_{1}_{2}_{3}'.format(skullstrip_tool,
                                                                                         motion_correct_ref,
                                                                                         motion_correct_tool,
                                                                                         num_strat)
                            )

                            node, out_file = new_strat['raw_functional_trunc']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.raw_func')

                            node, out_file = new_strat.get_leaf_properties()
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.func')

                            node, out_file = new_strat['anatomical_brain']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.anat_brain')

                            node, out_file = new_strat['anatomical_brain_mask']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.anatomical_brain_mask')

                            func_preproc.inputs.inputspec.twopass = \
                                getattr(c, 'functional_volreg_twopass', True)

                            new_strat.update_resource_pool({
                                    'movement_parameters': (func_preproc, 'outputspec.movement_parameters'),
                                    'max_displacement': (func_preproc, 'outputspec.max_displacement'),
                                    'functional_brain_mask_before_stc': (func_preproc, 'outputspec.mask'),
                                    'motion_correct_before_stc': (func_preproc, 'outputspec.motion_correct'),
                                    'coordinate_transformation': (func_preproc, 'outputspec.transform_matrices'),
                            })

                            gen_motion_stats = motion_power_statistics(
                                name = 'gen_motion_stats_before_stc_{0}_{1}_{2}_{3}'.format(skullstrip_tool,
                                                                                            motion_correct_ref,
                                                                                            motion_correct_tool,
                                                                                            num_strat),
                                motion_correct_tool = motion_correct_tool)

                            # Special case where the workflow is not getting outputs from
                            # resource pool but is connected to functional datasource
                            workflow.connect(func_wf, 'outputspec.subject',
                                            gen_motion_stats, 'inputspec.subject_id')

                            workflow.connect(func_wf, 'outputspec.scan',
                                            gen_motion_stats, 'inputspec.scan_id')

                            node, out_file = new_strat['motion_correct_before_stc']
                            workflow.connect(node, out_file,
                                            gen_motion_stats, 'inputspec.motion_correct')

                            node, out_file = new_strat['movement_parameters']
                            workflow.connect(node, out_file,
                                            gen_motion_stats,
                                            'inputspec.movement_parameters')

                            node, out_file = new_strat['max_displacement']
                            workflow.connect(node, out_file,
                                            gen_motion_stats, 'inputspec.max_displacement')

                            node, out_file = new_strat['functional_brain_mask_before_stc']
                            workflow.connect(node, out_file,
                                            gen_motion_stats, 'inputspec.mask')

                            node, out_file = new_strat['coordinate_transformation']
                            workflow.connect(node, out_file,
                                            gen_motion_stats, 'inputspec.transformations')

                            new_strat.append_name(gen_motion_stats.name)

                            new_strat.update_resource_pool({
                                'frame_wise_displacement_power': (gen_motion_stats, 'outputspec.FDP_1D'),
                                'frame_wise_displacement_jenkinson': (gen_motion_stats, 'outputspec.FDJ_1D'),
                                'dvars': (gen_motion_stats, 'outputspec.DVARS_1D'),
                                'power_params': (gen_motion_stats, 'outputspec.power_params'),
                                'motion_params': (gen_motion_stats, 'outputspec.motion_params')
                            })

                            new_strat_list.append(new_strat)

        strat_list = new_strat_list


        # Despike Workflow
        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            if 0 in c.runDespike:

                new_strat_list += [strat.fork()]

            if 1 in c.runDespike:

                new_strat = strat.fork()

                despike = pe.Node(interface=preprocess.Despike(),
                                name='func_despiked')
                despike.inputs.outputtype = 'NIFTI_GZ'

                node, out_file = new_strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                despike, 'in_file')

                new_strat.set_leaf_properties(despike, 'out_file')

                new_strat.update_resource_pool({
                    'despiked': (despike, 'out_file')
                })

                new_strat_list.append(new_strat)

        strat_list = new_strat_list


        # Slice Timing Correction Workflow
        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            if 0 in c.slice_timing_correction:

                new_strat_list += [strat.fork()]

            if 1 in c.slice_timing_correction:

                new_strat = strat.fork()

                slice_time = slice_timing_wf(name='func_slice_timing_correction_{0}'.format(num_strat))

                node, out_file = new_strat.get_leaf_properties()
                workflow.connect(node, out_file, slice_time,
                                'inputspec.func_ts')

                node, node_out = new_strat['tr']
                workflow.connect(node, node_out,
                                    slice_time, 'inputspec.tr')

                node, node_out = new_strat['tpattern']
                workflow.connect(node, node_out,
                                    slice_time, 'inputspec.tpattern')

                # add the name of the node to the strat name
                new_strat.append_name(slice_time.name)

                # set the leaf node
                new_strat.set_leaf_properties(slice_time, 'outputspec.slice_time_corrected')

                # add the outputs to the resource pool
                new_strat.update_resource_pool({
                    'slice_time_corrected': (slice_time, 'outputspec.slice_time_corrected')
                })

                new_strat_list.append(new_strat)

        # add new strats (if forked)
        strat_list = new_strat_list

        # Functional Image Preprocessing Workflow
        workflow, strat_list = connect_func_preproc(workflow, strat_list, c)

        # Distortion Correction
        new_strat_list = []

        # Distortion Correction - Field Map Phase-difference
        if "PhaseDiff" in c.distortion_correction and diff:
            for num_strat, strat in enumerate(strat_list):
                if 'BET' in c.fmap_distcorr_skullstrip:
                    epi_distcorr = create_EPI_DistCorr(
                        use_BET=True,
                        wf_name='diff_distcor_%d' % (num_strat)
                    )
                    epi_distcorr.inputs.bet_frac_input.bet_frac = c.fmap_distcorr_frac
                    epi_distcorr.get_node('bet_frac_input').iterables = \
                        ('bet_frac', c.fmap_distcorr_frac)
                else:
                    epi_distcorr = create_EPI_DistCorr(
                        use_BET=False,
                        wf_name='diff_distcor_%d' % (num_strat)
                    )
                    epi_distcorr.inputs.afni_threshold_input.afni_threshold = \
                        c.fmap_distcorr_threshold

                node, out_file = strat['anatomical_reorient']
                workflow.connect(node, out_file, epi_distcorr,
                                'inputspec.anat_file')

                node, out_file = strat['diff_phase']
                workflow.connect(node, out_file, epi_distcorr,
                                'inputspec.fmap_pha')

                node, out_file = strat['diff_mag_one']
                workflow.connect(node, out_file, epi_distcorr,
                                'inputspec.fmap_mag')

                node, out_file = strat['deltaTE']
                workflow.connect(node, out_file, epi_distcorr,
                                 'deltaTE_input.deltaTE')

                node, out_file = strat['diff_phase_dwell']
                workflow.connect(node, out_file, epi_distcorr,
                                 'dwellT_input.dwellT')

                node, out_file = strat['dwell_asym_ratio']
                workflow.connect(node, out_file, epi_distcorr,
                                 'dwell_asym_ratio_input.dwell_asym_ratio')

                # TODO ASH review forking
                if "None" in c.distortion_correction:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(epi_distcorr.name)

                strat.update_resource_pool({
                    'despiked_fieldmap': (epi_distcorr, 'outputspec.fmap_despiked'),
                    'fieldmap_mask': (epi_distcorr, 'outputspec.fieldmapmask'),
                })

        strat_list += new_strat_list

        diff_complete = False
        if "despiked_fieldmap" in strat and "fieldmap_mask" in strat:
            diff_complete = True

        # Distortion Correction - "Blip-Up / Blip-Down"
        if "Blip" in c.distortion_correction and blip:
            for num_strat, strat in enumerate(strat_list):
                match_epi_imports = ['import json']
                match_epi_fmaps_node = \
                    pe.Node(Function(input_names=['bold_pedir',
                                                  'epi_fmap_one',
                                                  'epi_fmap_params_one',
                                                  'epi_fmap_two',
                                                  'epi_fmap_params_two'],
                                     output_names=['opposite_pe_epi',
                                                   'same_pe_epi'],
                                     function=match_epi_fmaps,
                                     imports=match_epi_imports,
                                     as_module=True),
                            name='match_epi_fmaps_{0}'.format(num_strat))

                if fmap_rp_list:
                    epi_rp_key = fmap_rp_list[0]
                    epi_param_rp_key = "{0}_scan_params".format(epi_rp_key)
                    node, node_out = strat[epi_rp_key]
                    workflow.connect(node, node_out,
                                     match_epi_fmaps_node, 'epi_fmap_one')
                    node, node_out = strat[epi_param_rp_key]
                    workflow.connect(node, node_out,
                                     match_epi_fmaps_node, 'epi_fmap_params_one')
                    if len(epi_rp_key) > 1:
                        epi_rp_key = fmap_rp_list[1]
                        epi_param_rp_key = "{0}_scan_params".format(epi_rp_key)
                        node, node_out = strat[epi_rp_key]
                        workflow.connect(node, node_out,
                                         match_epi_fmaps_node, 'epi_fmap_two')
                        node, node_out = strat[epi_param_rp_key]
                        workflow.connect(node, node_out,
                                         match_epi_fmaps_node,
                                         'epi_fmap_params_two')

                node, node_out = strat['pe_direction']
                workflow.connect(node, node_out,
                                 match_epi_fmaps_node, 'bold_pedir')

                blip_correct = blip_distcor_wf(wf_name='blip_correct_{0}'.format(num_strat))

                node, out_file = strat["mean_functional"]
                workflow.connect(node, out_file,
                                 blip_correct, 'inputspec.func_mean')

                workflow.connect(match_epi_fmaps_node, 'opposite_pe_epi',
                                 blip_correct, 'inputspec.opposite_pe_epi')

                workflow.connect(match_epi_fmaps_node, 'same_pe_epi',
                                 blip_correct, 'inputspec.same_pe_epi')

            if "None" in c.distortion_correction:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(blip_correct.name)

            strat.update_resource_pool({
                'blip_warp': (blip_correct, 'outputspec.blip_warp'),
                'blip_warp_inverse': (blip_correct, 'outputspec.blip_warp_inverse'),
                'mean_functional': (blip_correct, 'outputspec.new_func_mean'),
                'functional_brain_mask': (blip_correct, 'outputspec.new_func_mask')
            }, override=True)

        strat_list += new_strat_list


        for num_strat, strat in enumerate(strat_list):

            # Resample brain mask with derivative resolution
            node, out_file = strat['functional_brain_mask']
            resampled_template = pe.Node(
                Function(
                    input_names=['resolution', 'template', 'template_name'],
                    output_names=['resampled_template'],
                    function=resolve_resolution,
                    as_module=True
                ),
                name='functional_brain_mask_derivative_%d' % (num_strat)
            )

            resampled_template.inputs.resolution = c.resolution_for_func_derivative
            resampled_template.inputs.template_name = 'functional_brain_mask_derivative'
            workflow.connect(node, out_file, resampled_template, 'template')
            strat.update_resource_pool({
                'functional_brain_mask_derivative': (
                    resampled_template,
                    'resampled_template'
                )
            })

        # Func -> T1 Registration (Initial Linear reg)

        # Depending on configuration, either passes output matrix to
        # Func -> Template ApplyWarp, or feeds into linear reg of BBReg operation
        # (if BBReg is enabled)

        new_strat_list = []

        if 1 in c.runRegisterFuncToAnat:

            for num_strat, strat in enumerate(strat_list):

                nodes = strat.get_nodes_names()

                # if field map-based distortion correction is on, but BBR is off,
                # send in the distortion correction files here
                # TODO: is this robust to the possibility of forking both
                # TODO: distortion correction and BBR at the same time?
                # TODO: (note if you are forking with BBR on/off, at this point
                # TODO:  there is still only one strat, so you would have to fork
                # TODO:  here instead to have a func->anat with fieldmap and
                # TODO:  without, and send the without-fieldmap to the BBR fork)

                dist_corr = False
                if 'diff_distcor' in nodes and 1 not in c.runBBReg:
                    dist_corr = True
                    # TODO: for now, disabling dist corr when BBR is disabled
                    err = "\n\n[!] Field map distortion correction is enabled, " \
                        "but Boundary-Based Registration is off- BBR is " \
                        "required for distortion correction.\n\n"
                    raise Exception(err)

                func_to_anat = create_register_func_to_anat(diff_complete,
                                                            'func_to_anat_FLIRT'
                                                            '_%d' % num_strat)

                # Input registration parameters
                func_to_anat.inputs.inputspec.interp = 'trilinear'


                # TODO ASH normalize strings with enums?
                if 'Mean Functional' in c.func_reg_input:
                    # Input functional image (mean functional)
                    node, out_file = strat['mean_functional']
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.func')

                elif 'Selected Functional Volume' in c.func_reg_input:
                    # Input functional image (specific volume)
                    node, out_file = strat['selected_func_volume']
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.func')

                # Input skull-stripped anatomical (anat.nii.gz)
                node, out_file = strat['anatomical_brain']
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.anat')

                if diff_complete:
                    # apply field map distortion correction outputs to
                    # the func->anat registration
                    node, out_file = strat['diff_phase_dwell']
                    workflow.connect(node, out_file,
                                     func_to_anat,
                                     'echospacing_input.echospacing')

                    node, out_file = strat['diff_phase_pedir']
                    workflow.connect(node, out_file,
                                     func_to_anat, 'pedir_input.pedir')

                    node, out_file = strat["despiked_fieldmap"]
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.fieldmap')

                    node, out_file = strat["fieldmap_mask"]
                    workflow.connect(node, out_file,
                                     func_to_anat, 'inputspec.fieldmapmask')

                # TODO ASH review forking
                if 0 in c.runRegisterFuncToAnat:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(func_to_anat.name)

                strat.update_resource_pool({
                    'mean_functional_in_anat': (func_to_anat, 'outputspec.anat_func_nobbreg'),
                    'functional_to_anat_linear_xfm': (func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg')
                })

        strat_list += new_strat_list

        # Func -> T1 Registration (BBREG)

        # Outputs 'functional_to_anat_linear_xfm', a matrix file of the
        # functional-to-anatomical registration warp to be applied LATER in
        # func_mni_warp, which accepts it as input 'premat'

        new_strat_list = []

        if 1 in c.runRegisterFuncToAnat and 1 in c.runBBReg:

            for num_strat, strat in enumerate(strat_list):

                nodes = strat.get_nodes_names()

                # this is needed here in case tissue segmentation is set on/off
                # and you have bbreg enabled- this will ensure bbreg will run for
                # the strat that has segmentation but will not run (thus avoiding
                # a crash) on the strat without segmentation
                if 'seg_preproc' in nodes or 'seg_preproc_t1_template' in nodes:

                    func_to_anat_bbreg = create_bbregister_func_to_anat(
                        diff_complete,
                        'func_to_anat_bbreg_%d' % num_strat
                    )

                    # Input registration parameters
                    func_to_anat_bbreg.inputs.inputspec.bbr_schedule = \
                        c.boundaryBasedRegistrationSchedule

                    # TODO ASH normalize strings with enums?
                    if 'Mean Functional' in c.func_reg_input:
                        # Input functional image (mean functional)
                        node, out_file = strat['mean_functional']
                        workflow.connect(node, out_file,
                                        func_to_anat_bbreg, 'inputspec.func')

                    elif 'Selected Functional Volume' in c.func_reg_input:
                        # Input functional image (specific volume)
                        node, out_file = strat['selected_func_volume']
                        workflow.connect(node, out_file,
                                        func_to_anat_bbreg, 'inputspec.func')

                    # Input anatomical whole-head image (reoriented)
                    node, out_file = strat['anatomical_reorient']
                    workflow.connect(node, out_file,
                                    func_to_anat_bbreg,
                                    'inputspec.anat_skull')

                    node, out_file = strat['functional_to_anat_linear_xfm']
                    workflow.connect(node, out_file,
                                    func_to_anat_bbreg,
                                    'inputspec.linear_reg_matrix')

                    if 'T1_template' in c.template_based_segmentation or 'EPI_template' in c.template_based_segmentation:
                        # Input segmentation mask,
                        # since template_based_segmentation cannot generate probability maps
                        node, out_file = strat['anatomical_wm_mask']
                        workflow.connect(node, out_file,
                                        func_to_anat_bbreg,
                                        'inputspec.anat_wm_segmentation')

                    else:
                        # Input segmentation probability maps for white matter
                        # segmentation
                        node, out_file = strat['seg_probability_maps']
                        workflow.connect(node, (out_file, pick_wm),
                                         func_to_anat_bbreg,
                                         'inputspec.anat_wm_segmentation')

                    #if dist_corr and phase_diff:
                    # apply field map distortion correction outputs to
                    # the func->anat registration
                    if diff_complete:
                        node, out_file = strat['diff_phase_dwell']
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg,
                                         'echospacing_input.echospacing')

                        node, out_file = strat['diff_phase_pedir']
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg,
                                         'pedir_input.pedir')

                        node, out_file = strat["despiked_fieldmap"]
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg,
                                         'inputspec.fieldmap')

                        node, out_file = strat["fieldmap_mask"]
                        workflow.connect(node, out_file,
                                         func_to_anat_bbreg,
                                         'inputspec.fieldmapmask')


                    # TODO ASH review forking
                    if 0 in c.runBBReg:
                        strat = strat.fork()
                        new_strat_list.append(strat)

                    strat.append_name(func_to_anat_bbreg.name)

                    strat.update_resource_pool({
                        'mean_functional_in_anat': (func_to_anat_bbreg, 'outputspec.anat_func'),
                        'functional_to_anat_linear_xfm': (func_to_anat_bbreg, 'outputspec.func_to_anat_linear_xfm')
                    }, override=True)

                else:
                    # TODO ASH review
                    # anatomical segmentation is not being run in this particular
                    # strategy/fork - we don't want this to stop workflow building
                    # unless there is only one strategy
                    if len(strat_list) > 1:
                        pass
                    else:
                        err = "\n\n[!] Boundary-based registration (BBR) for " \
                            "functional-to-anatomical registration is " \
                            "enabled, but anatomical segmentation is not. " \
                            "BBR requires the outputs of segmentation. " \
                            "Please modify your pipeline configuration and " \
                            "run again.\n\n"
                        raise Exception(err)


        strat_list += new_strat_list

        # CC This is the first opportunity to write some of the outputs of basic
        # func preproc, such as the brain mask and mean EPI. Doing it any later
        # might result in multiple versions of these files being needlessly generated
        # do to strategies created by denoising, which do not impact the mean or brainmask
        # preproc Func -> T1/EPI Template

        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            if 'EPI_template' in c.runRegisterFuncToTemplate:

                for reg in c.regOption:

                    if 'T1_template' in c.runRegisterFuncToTemplate:
                        strat = strat.fork()

                    func_to_epi = \
                        create_register_func_to_epi(
                            name='func_to_epi_{0}_{1}'.format(reg.lower(), num_strat),
                            reg_option=reg,
                            reg_ants_skull=c.regWithSkull
                        )

                    # Input registration parameters
                    if c.ANTs_para_EPI_registration is None:
                        err_msg = '\n\n[!] C-PAC says: \n'\
                            "You have selected \'regOption: [{0}]\' and \'runRegisterFuncToTemplate :  ['{1}']\'. \n"\
                                 'However, no EPI-to-template ANTs parameters were specified. ' \
                                    'Please specify ANTs parameters properly and try again'.format(str(c.regOption),
                                                                                                   str(c.runRegisterFuncToTemplate))
                        raise Exception(err_msg)
                    else:
                        func_to_epi.inputs.inputspec.ants_para = c.ANTs_para_EPI_registration

                    func_to_epi.inputs.inputspec.interp = c.funcRegANTSinterpolation

                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file, func_to_epi, 'inputspec.func_4d')

                    if 'Mean Functional' in c.func_reg_input:
                        node, out_file = strat['mean_functional']
                        workflow.connect(node, out_file, func_to_epi, 'inputspec.func_3d')

                    elif 'Selected Functional Volume' in c.func_reg_input:
                        node, out_file = strat['selected_func_volume']
                        workflow.connect(node, out_file, func_to_epi, 'inputspec.func_3d')

                    node, out_file = strat['template_epi']
                    workflow.connect(node, out_file, func_to_epi, 'inputspec.epi')

                    node, out_file = strat['functional_brain_mask']
                    workflow.connect(node, out_file, func_to_epi, 'inputspec.func_3d_mask')

                    # update resource pool
                    strat.update_resource_pool({
                        'functional_to_epi-standard': (func_to_epi, 'outputspec.func_in_epi'),
                    })

                    if reg == 'FSL':
                        strat.update_resource_pool({
                            'func_to_epi_linear_xfm': (func_to_epi, 'outputspec.fsl_flirt_xfm'),  
                            'func_to_epi_nonlinear_xfm': (func_to_epi, 'outputspec.fsl_fnirt_xfm'),
                            'epi_to_func_linear_xfm': (func_to_epi, 'outputspec.invlinear_xfm'),
                        })

                    elif reg == 'ANTS':
                        strat.update_resource_pool({
                            'func_to_epi_ants_initial_xfm': (func_to_epi, 'outputspec.ants_initial_xfm'),
                            'func_to_epi_ants_rigid_xfm': (func_to_epi, 'outputspec.ants_rigid_xfm'),
                            'func_to_epi_ants_affine_xfm': (func_to_epi, 'outputspec.ants_affine_xfm'),
                            'func_to_epi_nonlinear_xfm': (func_to_epi, 'outputspec.warp_field'),
                            'epi_to_func_nonlinear_xfm': (func_to_epi, 'outputspec.inverse_warp_field'), # rename
                        })

                    strat.append_name(func_to_epi.name)

                    for output_name, func_key, ref_key, image_type in [ \
                            ('functional_brain_mask_to_standard', 'functional_brain_mask', 'template_skull_for_func_preproc', 'func_mask'),
                            ('functional_brain_mask_to_standard_derivative', 'functional_brain_mask', 'template_skull_for_func_derivative', 'func_mask'),
                            ('mean_functional_to_standard', 'mean_functional', 'template_brain_for_func_preproc', 'func_derivative'),
                            ('mean_functional_to_standard_derivative', 'mean_functional', 'template_brain_for_func_derivative', 'func_derivative'),
                            ('motion_correct_to_standard', 'motion_correct', 'template_brain_for_func_preproc', 'func_4d'),
                    ]:
                        output_func_to_standard(workflow, func_key, ref_key,
                                                output_name, strat,
                                                num_strat, c,
                                                input_image_type=image_type,
                                                registration_template='epi',
                                                func_type='non-ica-aroma')

                    if 'T1_template' in c.runRegisterFuncToTemplate:
                        new_strat_list.append(strat)

        strat_list += new_strat_list


        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            if 'T1_template' in c.runRegisterFuncToTemplate and \
                    'func_to_epi_ants' not in nodes and \
                    'func_to_epi_fsl' not in nodes:

                for output_name, func_key, ref_key, image_type in [ \
                        ('functional_brain_mask_to_standard', 'functional_brain_mask', 'template_skull_for_func_preproc', 'func_mask'),
                        ('functional_brain_mask_to_standard_derivative', 'functional_brain_mask', 'template_skull_for_func_derivative', 'func_mask'),
                        ('mean_functional_to_standard', 'mean_functional', 'template_brain_for_func_preproc', 'func_derivative'),
                        ('mean_functional_to_standard_derivative', 'mean_functional', 'template_brain_for_func_derivative', 'func_derivative'),
                        ('motion_correct_to_standard', 'motion_correct', 'template_brain_for_func_preproc', 'func_4d'),
                ]:
                    output_func_to_standard(workflow, func_key, ref_key,
                                            output_name, strat, num_strat, c,
                                            input_image_type=image_type,
                                            registration_template='t1',
                                            func_type='non-ica-aroma')


        # Inserting epi-template-based-segmentation Workflow
        new_strat_list = []

        if 'EPI_template' in c.template_based_segmentation:

            for num_strat, strat in enumerate(strat_list):

                nodes = strat.get_nodes_names()

                if 'func_to_epi_fsl' not in nodes and 'func_to_epi_ants' not in nodes:
                    continue

                if not any(o in c.template_based_segmentation for o in ['EPI_template', 'T1_template', 'None']):
                    err = '\n\n[!] C-PAC says: Your template based segmentation ' \
                        'setting does not include either \'EPI_template\' or \'T1_template\'.\n\n' \
                        'Options you provided:\ntemplate_based_segmentation: {0}' \
                        '\n\n'.format(str(c.template_based_segmentation))
                    raise Exception(err)

                # TODO ASH based on config, instead of nodes?
                if 'func_to_epi_fsl' in nodes:
                    use_ants = False
                elif 'func_to_epi_ants' in nodes:
                    use_ants = True

                seg_preproc_template_based = create_seg_preproc_template_based(use_ants=use_ants,
                                                                wf_name='seg_preproc_epi_template_{0}'.format(num_strat))

                # TODO ASH review
                if seg_preproc_template_based is None:
                    continue

                if 'func_to_epi_fsl' in nodes:

                    node, out_file = strat['mean_functional']
                    workflow.connect(node, out_file,
                                    seg_preproc_template_based, 'inputspec.brain')

                    node, out_file = strat['func_to_epi_invlinear_xfm']
                    workflow.connect(node, out_file,
                                    seg_preproc_template_based,
                                    'inputspec.standard2highres_mat')

                elif 'func_to_epi_ants' in nodes:

                    node, out_file = strat['mean_functional']
                    workflow.connect(node, out_file,
                                    seg_preproc_template_based, 'inputspec.brain')

                    node, out_file = strat['func_to_epi_ants_initial_xfm']
                    # node, out_file = strat['ants_initial_xfm']
                    workflow.connect(node, out_file,
                                    seg_preproc_template_based,
                                    'inputspec.standard2highres_init')

                    node, out_file = strat['func_to_epi_ants_rigid_xfm']
                    # node, out_file = strat['ants_rigid_xfm']
                    workflow.connect(node, out_file,
                                    seg_preproc_template_based,
                                    'inputspec.standard2highres_rig')

                    node, out_file = strat['func_to_epi_ants_affine_xfm']
                    # node, out_file = strat['ants_affine_xfm']
                    workflow.connect(node, out_file,
                                    seg_preproc_template_based,
                                    'inputspec.standard2highres_mat')

                workflow.connect(c.template_based_segmentation_CSF, 'local_path',
                                    seg_preproc_template_based, 'inputspec.CSF_template')

                workflow.connect(c.template_based_segmentation_GRAY, 'local_path',
                                    seg_preproc_template_based, 'inputspec.GRAY_template')

                workflow.connect(c.template_based_segmentation_WHITE, 'local_path',
                                    seg_preproc_template_based, 'inputspec.WHITE_template')

                # TODO ASH review with forking function
                if 'None' in c.template_based_segmentation:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(seg_preproc_template_based.name)
                strat.update_resource_pool({
                    'epi_gm_mask': (seg_preproc_template_based, 'outputspec.gm_mask'),
                    'epi_csf_mask': (seg_preproc_template_based, 'outputspec.csf_mask'),
                    'epi_wm_mask': (seg_preproc_template_based, 'outputspec.wm_mask')
                })

        strat_list += new_strat_list


        # Inserting Generate Motion Statistics Workflow
        new_strat_list = []
        for num_strat, strat in enumerate(strat_list):
            nodes = strat.get_nodes_names()
            if not any("gen_motion_stats_before_stc" in node for node in nodes):

                motion_stats_node = [x for x in nodes if "func_preproc" in x][0]

                # skullstripping tool
                if "fsl_afni" in motion_stats_node:
                    skullstrip_tool = 'fsl_afni'
                elif "fsl" in motion_stats_node:
                    skullstrip_tool = 'fsl'
                elif "afni" in motion_stats_node:
                    skullstrip_tool = 'afni'
                elif "anatomical_refined" in motion_stats_node:
                    skullstrip_tool = 'anatomical_refined'

                # motion correction reference
                if "mean" in motion_stats_node:
                    motion_correct_ref = "mean"
                elif "median" in motion_stats_node:
                    motion_correct_ref = "median"
                elif "selected_volume" in motion_stats_node:
                    motion_correct_ref = "selected_volume"

                # motion correction tool
                if "3dvolreg" in motion_stats_node:
                    motion_correct_tool = "3dvolreg"
                elif "mcflirt" in motion_stats_node:
                    motion_correct_tool = "mcflirt"

                gen_motion_stats = motion_power_statistics(
                                name='gen_motion_stats_{0}_{1}_{2}_{3}'.format(skullstrip_tool,
                                                                               motion_correct_ref,
                                                                               motion_correct_tool,
                                                                               num_strat),
                                motion_correct_tool=motion_correct_tool)

                # Special case where the workflow is not getting outputs from
                # resource pool but is connected to functional datasource
                workflow.connect(func_wf, 'outputspec.subject',
                                 gen_motion_stats, 'inputspec.subject_id')

                workflow.connect(func_wf, 'outputspec.scan',
                                 gen_motion_stats, 'inputspec.scan_id')

                node, out_file = strat['motion_correct']
                workflow.connect(node, out_file,
                                 gen_motion_stats, 'inputspec.motion_correct')

                node, out_file = strat['movement_parameters']
                workflow.connect(node, out_file,
                                 gen_motion_stats,
                                 'inputspec.movement_parameters')

                node, out_file = strat['max_displacement']
                workflow.connect(node, out_file,
                                 gen_motion_stats,
                                 'inputspec.max_displacement')

                node, out_file = strat['functional_brain_mask']
                workflow.connect(node, out_file,
                                 gen_motion_stats, 'inputspec.mask')

                node, out_file = strat['coordinate_transformation']
                workflow.connect(node, out_file,
                                 gen_motion_stats,
                                 'inputspec.transformations')

                strat.update_resource_pool({
                    'frame_wise_displacement_power': (gen_motion_stats, 'outputspec.FDP_1D'),
                    'frame_wise_displacement_jenkinson': (gen_motion_stats, 'outputspec.FDJ_1D'),
                    'dvars': (gen_motion_stats, 'outputspec.DVARS_1D'),
                    'power_params': (gen_motion_stats, 'outputspec.power_params'),
                    'motion_params': (gen_motion_stats, 'outputspec.motion_params')
                })


        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            if 1 in c.runICA:

                if 0 in c.runICA:
                    new_strat_list += [strat.fork()]

                nodes = strat.get_nodes_names()

                if 'none' in str(c.TR).lower():
                    TR = None
                else:
                    TR = float(c.TR)

                # FNIRT ONLY! ANTS further below!
                if 'FSL' in c.regOption and \
                        'anat_symmetric_mni_ants_register' not in nodes and \
                            'anat_mni_ants_register' not in nodes:

                    aroma_preproc = create_aroma(tr=TR,
                                                wf_name='create_aroma_%d' % num_strat)

                    aroma_preproc.inputs.params.denoise_type = c.aroma_denoise_type

                    node, out_file = strat.get_leaf_properties()
                    workflow.connect(node, out_file, aroma_preproc,
                                    'inputspec.denoise_file')

                    node, out_file = strat['functional_to_anat_linear_xfm']
                    workflow.connect(node, out_file, aroma_preproc,
                                    'inputspec.mat_file')

                    node, out_file = strat['anatomical_to_mni_nonlinear_xfm']
                    workflow.connect(node, out_file, aroma_preproc,
                                    'inputspec.fnirt_warp_file')

                    if c.aroma_denoise_type == 'nonaggr':

                        strat.set_leaf_properties(aroma_preproc,
                                                  'outputspec.nonaggr_denoised_file')

                        strat.update_resource_pool({
                            'ica_aroma_denoised_functional': (
                                aroma_preproc, 'outputspec.nonaggr_denoised_file')
                            }
                        )

                    elif c.aroma_denoise_type == 'aggr':
                        strat.set_leaf_properties(aroma_preproc,
                                                  'outputspec.aggr_denoised_file')

                        strat.update_resource_pool({
                            'ica_aroma_denoised_functional': (
                                aroma_preproc, 'outputspec.aggr_denoised_file')
                            }
                        )

                    strat.append_name(aroma_preproc.name)

                elif 'ANTS' in c.regOption and \
                    'anat_symmetric_mni_flirt_register' not in nodes and \
                    'anat_symmetric_mni_fnirt_register' not in nodes and \
                    'anat_mni_flirt_register' not in nodes and \
                    'anat_mni_fnirt_register' not in nodes:

                    # we don't have the FNIRT warp file, so we need to calculate
                    # ICA-AROMA de-noising in template space

                    if 'func_to_epi_ants' in nodes:

                        for output_name, func_key, ref_key, image_type in [ \
                                ('ica_aroma_functional_to_standard', 'leaf', 'template_brain_for_func_preproc', 'func_4d'),
                        ]:
                            output_func_to_standard(workflow, func_key, ref_key, output_name, strat, num_strat, c, input_image_type=image_type, registration_template='epi', func_type='ica-aroma')

                    elif 'T1_template' in c.runRegisterFuncToTemplate:

                        for output_name, func_key, ref_key, image_type in [ \
                                ('ica_aroma_functional_to_standard', 'leaf', 'template_brain_for_func_preproc', 'func_4d'),
                        ]:
                            output_func_to_standard(workflow, func_key, ref_key, output_name, strat, num_strat, c, input_image_type=image_type, registration_template='t1',func_type='ica-aroma')

                    aroma_preproc = create_aroma(tr=TR, wf_name='create_aroma_{0}'.format(num_strat))
                    aroma_preproc.inputs.params.denoise_type = c.aroma_denoise_type

                    node, out_file = strat['ica_aroma_functional_to_standard']
                    workflow.connect(node, out_file, aroma_preproc,
                                    'inputspec.denoise_file')

                    # warp back
                    if c.aroma_denoise_type == 'nonaggr':
                        node, out_file = (
                            aroma_preproc, 'outputspec.nonaggr_denoised_file'
                        )

                    elif c.aroma_denoise_type == 'aggr':
                        node, out_file = (
                            aroma_preproc, 'outputspec.aggr_denoised_file'
                        )

                    strat.update_resource_pool({
                        'ica_aroma_denoised_functional_to_standard': (node, out_file)
                        }
                    )

                    if 'func_to_epi_ants' in nodes:

                        for output_name, func_key, ref_key, image_type in [ \
                                ('ica_aroma_denoised_functional', 'ica_aroma_denoised_functional_to_standard', 'mean_functional', 'func_4d'),
                        ]:
                            output_func_to_standard(workflow, func_key, ref_key, output_name, strat, num_strat, c, input_image_type=image_type, inverse=True, registration_template='epi', func_type='ica-aroma')

                    else:

                        for output_name, func_key, ref_key, image_type in [ \
                                ('ica_aroma_denoised_functional', 'ica_aroma_denoised_functional_to_standard', 'mean_functional', 'func_4d'),
                        ]:
                            output_func_to_standard(workflow, func_key, ref_key, output_name, strat, num_strat, c, input_image_type=image_type, inverse=True, registration_template='t1', func_type='ica-aroma')

                    node, out_file = strat["ica_aroma_denoised_functional"]
                    strat.set_leaf_properties(node, out_file)

                    strat.update_resource_pool({
                        'ica_aroma_denoised_functional': (node, out_file)
                        }, override=True
                    )

        strat_list += new_strat_list

        # Inserting Nuisance Regressor Workflow
        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            # for each strategy, create a new one without nuisance
            if 0 in c.runNuisance or 1 in c.run_pypeer:
                new_strat_list.append(strat.fork())

            nodes = strat.get_nodes_names()

            has_segmentation = 'seg_preproc' in nodes or 'seg_preproc_t1_template' in nodes or 'seg_preproc_epi_template' in nodes
            use_ants = 'anat_mni_fnirt_register' not in nodes and 'anat_mni_flirt_register' not in nodes and 'func_to_epi_fsl' not in nodes

            for regressors_selector_i, regressors_selector in enumerate(c.Regressors):

                new_strat = strat.fork()

                # Before start nuisance_wf, covert OrderedDict(regressors_selector) to dict
                from CPAC.utils.utils import ordereddict_to_dict
                regressors_selector = ordereddict_to_dict(regressors_selector)

                # to guarantee immutability
                regressors_selector = NuisanceRegressor(
                    copy.deepcopy(regressors_selector)
                )

                # remove tissue regressors when there is no segmentation
                # on the strategy
                if not has_segmentation:
                    for reg in ['aCompCor',
                                'WhiteMatter',
                                'GreyMatter',
                                'CerebrospinalFluid']:

                        if reg in regressors_selector:
                            del regressors_selector[reg]

                regressor_workflow = create_regressor_workflow(
                    regressors_selector,
                    use_ants=use_ants,
                    ventricle_mask_exist=ventricle_mask_exist,
                    name='nuisance_regressor_{0}_{1}'.format(regressors_selector_i, num_strat)
                )

                node, node_out = strat['tr']
                workflow.connect(node, node_out,
                                 regressor_workflow, 'inputspec.tr')

                node, out_file = new_strat['anatomical_brain']
                workflow.connect(
                    node, out_file,
                    regressor_workflow, 'inputspec.anatomical_file_path'
                )

                if has_segmentation:

                    workflow.connect(
                        c.lateral_ventricles_mask, 'local_path',
                        regressor_workflow, 'inputspec.lat_ventricles_mask_file_path'
                    )

                    if 'seg_preproc' in nodes or 'seg_preproc_t1_template' in nodes :

                        node, out_file = new_strat['anatomical_gm_mask']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow, 'inputspec.gm_mask_file_path'
                        )

                        node, out_file = new_strat['anatomical_wm_mask']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow, 'inputspec.wm_mask_file_path'
                        )

                        node, out_file = new_strat['anatomical_csf_mask']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow, 'inputspec.csf_mask_file_path'
                        )

                    if 'seg_preproc_epi_template' in nodes:

                        node, out_file = new_strat['epi_gm_mask']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow, 'inputspec.gm_mask_file_path'
                        )

                        node, out_file = new_strat['epi_wm_mask']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow, 'inputspec.wm_mask_file_path'
                        )

                        node, out_file = new_strat['epi_csf_mask']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow, 'inputspec.csf_mask_file_path'
                        )

                node, out_file = new_strat['movement_parameters']
                workflow.connect(
                    node, out_file,
                    regressor_workflow,
                    'inputspec.motion_parameters_file_path'
                )

                node, out_file= new_strat['functional_to_anat_linear_xfm']
                workflow.connect(
                    node, out_file,
                    regressor_workflow,
                    'inputspec.func_to_anat_linear_xfm_file_path'
                )

                node, out_file = new_strat.get_leaf_properties()
                workflow.connect(
                    node, out_file,
                    regressor_workflow,
                    'inputspec.functional_file_path'
                )

                new_strat.update_resource_pool({
                    'functional_freq_unfiltered': (
                        node, out_file
                    ),
                })

                node, out_file = new_strat['frame_wise_displacement_jenkinson']
                workflow.connect(
                    node, out_file,
                    regressor_workflow,
                    'inputspec.fd_j_file_path'
                )

                node, out_file = new_strat['frame_wise_displacement_power']
                workflow.connect(
                    node, out_file,
                    regressor_workflow,
                    'inputspec.fd_p_file_path'
                )

                node, out_file = new_strat['dvars']
                workflow.connect(
                    node, out_file,
                    regressor_workflow,
                    'inputspec.dvars_file_path'
                )

                node, out_file = new_strat['functional_brain_mask']
                workflow.connect(
                    node, out_file,
                    regressor_workflow,
                    'inputspec.functional_brain_mask_file_path'
                )

                if c.brain_use_erosion:
                    node, out_file = new_strat['anatomical_eroded_brain_mask']
                    workflow.connect(
                        node, out_file,
                        regressor_workflow, 'inputspec.anatomical_eroded_brain_mask_file_path'
                    )

                regressor_workflow.get_node('inputspec').iterables = ([
                    ('selector', [regressors_selector]),
                ])

                if use_ants:
                    if 'func_to_epi_ants' in nodes:
                        # pass the ants_affine_xfm to the input for the
                        # INVERSE transform, but ants_affine_xfm gets inverted
                        # within the workflow

                        node, out_file = new_strat['func_to_epi_ants_initial_xfm']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow,
                            'inputspec.anat_to_mni_initial_xfm_file_path'
                        )

                        node, out_file = new_strat['func_to_epi_ants_rigid_xfm']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow,
                            'inputspec.anat_to_mni_rigid_xfm_file_path'
                        )

                        node, out_file = new_strat['func_to_epi_ants_affine_xfm']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow,
                            'inputspec.anat_to_mni_affine_xfm_file_path'
                        )

                    elif 'T1_template' in c.runRegisterFuncToTemplate:
                        # pass the ants_affine_xfm to the input for the
                        # INVERSE transform, but ants_affine_xfm gets inverted
                        # within the workflow

                        node, out_file = new_strat['ants_initial_xfm']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow,
                            'inputspec.anat_to_mni_initial_xfm_file_path'
                        )

                        node, out_file = new_strat['ants_rigid_xfm']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow,
                            'inputspec.anat_to_mni_rigid_xfm_file_path'
                        )

                        node, out_file = new_strat['ants_affine_xfm']
                        workflow.connect(
                            node, out_file,
                            regressor_workflow,
                            'inputspec.anat_to_mni_affine_xfm_file_path'
                        )
                else:
                    node, out_file = new_strat['mni_to_anatomical_linear_xfm']
                    workflow.connect(
                        node, out_file,
                        regressor_workflow,
                        'inputspec.mni_to_anat_linear_xfm_file_path'
                    )

                new_strat.update_resource_pool({
                    'nuisance_regression_selector': regressors_selector,
                    'functional_nuisance_regressors': (
                        regressor_workflow,
                        'outputspec.regressors_file_path'
                    ),
                })

                # Inserting Nuisance REGRESSION Workflow
                if 1 in c.runNuisance:

                    if 'Bandpass' in regressors_selector:
                        nuis_name = 'nuisance_regression_before-filt_{0}_' \
                                    '{1}'.format(regressors_selector_i, num_strat)
                    else:
                        nuis_name = 'nuisance_regression_{0}_' \
                                    '{1}'.format(regressors_selector_i, num_strat)

                    nuisance_regression_before_workflow = create_nuisance_regression_workflow(
                        regressors_selector,
                        name=nuis_name)

                    if 'Bandpass' in regressors_selector:
                        filtering = filtering_bold_and_regressors(regressors_selector,
                                                                  name='frequency_filtering_'
                                                                       '{0}_{1}'.format(regressors_selector_i, num_strat))

                    node, out_file = new_strat.get_leaf_properties()

                    workflow.connect(
                        node, out_file,
                        nuisance_regression_before_workflow,
                        'inputspec.functional_file_path'
                    )

                    if 'Bandpass' in regressors_selector:
                        workflow.connect(
                            regressor_workflow,
                            'outputspec.regressors_file_path',
                            filtering,
                            'inputspec.regressors_file_path'
                        )

                    workflow.connect(
                        regressor_workflow,
                        'outputspec.regressors_file_path',
                        nuisance_regression_before_workflow,
                        'inputspec.regressor_file'
                    )

                    node, out_file = new_strat['functional_brain_mask']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_before_workflow,
                        'inputspec.functional_brain_mask_file_path'
                    )

                    node, out_file = new_strat['frame_wise_displacement_jenkinson']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_before_workflow,
                        'inputspec.fd_j_file_path'
                    )

                    node, out_file = new_strat['frame_wise_displacement_power']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_before_workflow,
                        'inputspec.fd_p_file_path'
                    )

                    node, out_file = new_strat['dvars']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_before_workflow,
                        'inputspec.dvars_file_path'
                    )

                    if 'Bandpass' in regressors_selector:
                        if 'Before' in c.filtering_order:
                            nuisance_regression_after_workflow = create_nuisance_regression_workflow(
                                regressors_selector,
                                name='nuisance_regression_after-filt_{0}_'
                                     '{1}'.format(regressors_selector_i, num_strat))

                            workflow.connect(
                                filtering,
                                'outputspec.residual_file_path',
                                nuisance_regression_after_workflow,
                                'inputspec.functional_file_path'
                            )

                            workflow.connect(
                                filtering,
                                'outputspec.residual_regressor',
                                nuisance_regression_after_workflow,
                                'inputspec.regressor_file'
                            )

                            node, out_file = new_strat['functional_brain_mask']
                            workflow.connect(
                                node, out_file,
                                nuisance_regression_after_workflow,
                                'inputspec.functional_brain_mask_file_path'
                            )

                            node, out_file = new_strat['frame_wise_displacement_jenkinson']
                            workflow.connect(
                                node, out_file,
                                nuisance_regression_after_workflow,
                                'inputspec.fd_j_file_path'
                            )

                            node, out_file = new_strat['frame_wise_displacement_power']
                            workflow.connect(
                                node, out_file,
                                nuisance_regression_after_workflow,
                                'inputspec.fd_p_file_path'
                            )

                            node, out_file = new_strat['dvars']
                            workflow.connect(
                                node, out_file,
                                nuisance_regression_after_workflow,
                                'inputspec.dvars_file_path'
                            )

                            node, out_file = new_strat.get_leaf_properties()
                            workflow.connect(
                                node, out_file,
                                filtering,
                                'inputspec.functional_file_path'
                            )

                            new_strat.set_leaf_properties(
                                nuisance_regression_after_workflow,
                                'outputspec.residual_file_path'
                            )

                            new_strat.update_resource_pool({
                                'functional_freq_filtered': (
                                    filtering,
                                    'outputspec.residual_file_path'
                                ),
                            })

                            new_strat.update_resource_pool({
                                 'functional_nuisance_residuals': (
                                    nuisance_regression_after_workflow,
                                    'outputspec.residual_file_path'
                                ),
                            })

                            new_strat.append_name(nuisance_regression_after_workflow.name)

                        elif 'After' in c.filtering_order:
                            workflow.connect(
                                nuisance_regression_before_workflow,
                                'outputspec.residual_file_path',
                                filtering,
                                'inputspec.functional_file_path'
                            )

                            new_strat.set_leaf_properties(
                                filtering,
                                'outputspec.residual_file_path'
                            )

                            new_strat.update_resource_pool({
                                'functional_nuisance_residuals': (
                                    nuisance_regression_before_workflow,
                                    'outputspec.residual_file_path'
                                ),
                            })

                            new_strat.update_resource_pool({
                                'functional_freq_filtered': (
                                    filtering,
                                    'outputspec.residual_file_path'
                                ),
                            })

                    else:
                        new_strat.set_leaf_properties(
                            nuisance_regression_before_workflow,
                            'outputspec.residual_file_path'
                        )

                        new_strat.update_resource_pool({
                            'functional_nuisance_residuals': (
                                nuisance_regression_before_workflow,
                                'outputspec.residual_file_path'
                            ),
                        })

                    new_strat.update_resource_pool({
                        'functional_freq_unfiltered': (
                            nuisance_regression_before_workflow,
                            'outputspec.residual_file_path'
                        ),
                    }, override=True)

                    new_strat.append_name(regressor_workflow.name)
                    new_strat.append_name(nuisance_regression_before_workflow.name)
                    new_strat.append_name(filtering.name)

                new_strat_list.append(new_strat)

        # Be aware that this line is supposed to override the current strat_list: it is not a typo/mistake!
        # Each regressor forks the strategy, instead of reusing it, to keep the code simple
        strat_list = new_strat_list


        # Inserting Median Angle Correction Workflow
        new_strat_list = []

        # TODO ASH normalize w schema val
        if 1 in c.runMedianAngleCorrection:

            for num_strat, strat in enumerate(strat_list):

                # for each strategy, create a new one without median angle
                if 0 in c.runMedianAngleCorrection:
                    new_strat_list.append(strat.fork())

                median_angle_corr = create_median_angle_correction(
                    'median_angle_corr_%d' % num_strat
                )

                median_angle_corr.get_node('median_angle_correct').iterables = \
                    ('target_angle_deg', c.targetAngleDeg)

                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                median_angle_corr, 'inputspec.subject')

                strat.append_name(median_angle_corr.name)

                strat.set_leaf_properties(median_angle_corr, 'outputspec.subject')

                strat.update_resource_pool({
                    'functional_median_angle_corrected': (median_angle_corr, 'outputspec.subject')
                })

        strat_list += new_strat_list


        # Denoised Func -> Template, uses antsApplyTransforms (ANTS) or ApplyWarp (FSL) to
        #  apply the warp; also includes mean functional warp
        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            if 'func_to_epi_fsl' in nodes or 'func_to_epi_ants' in nodes :
                for output_name, func_key, ref_key, image_type in [ \
                        ('functional_to_standard', 'leaf', 'template_brain_for_func_preproc', 'func_4d'),
                ]:
                    output_func_to_standard(workflow, func_key, ref_key, output_name, strat, num_strat, c, input_image_type=image_type, registration_template='epi', func_type='non-ica-aroma')

            elif 'T1_template' in c.runRegisterFuncToTemplate:
                for output_name, func_key, ref_key, image_type in [ \
                        ('functional_to_standard', 'leaf', 'template_brain_for_func_preproc', 'func_4d'),
                ]:
                    output_func_to_standard(workflow, func_key, ref_key, output_name, strat, num_strat, c, input_image_type=image_type, registration_template='t1', func_type='non-ica-aroma')

        strat_list += new_strat_list


        # Derivatives

        # Inserting ALFF/fALFF workflow
        #     NOTE: this is calculated using the functional time series from
        #           before frequency filtering and beyond
        new_strat_list = []

        if 1 in c.runALFF:
            for num_strat, strat in enumerate(strat_list):

                alff = create_alff('alff_falff_{0}'.format(num_strat))

                alff.inputs.hp_input.hp = c.highPassFreqALFF
                alff.inputs.lp_input.lp = c.lowPassFreqALFF
                alff.get_node('hp_input').iterables = ('hp',
                                                    c.highPassFreqALFF)
                alff.get_node('lp_input').iterables = ('lp',
                                                    c.lowPassFreqALFF)

                node, out_file = strat['functional_freq_unfiltered']
                workflow.connect(node, out_file,
                                 alff, 'inputspec.rest_res')
                node, out_file = strat['functional_brain_mask']
                workflow.connect(node, out_file,
                                 alff, 'inputspec.rest_mask')

                strat.append_name(alff.name)

                strat.update_resource_pool({
                    'alff': (alff, 'outputspec.alff_img'),
                    'falff': (alff, 'outputspec.falff_img')
                })

        strat_list += new_strat_list

        # Inserting VMHC Workflow

        new_strat_list = []

        if 1 in c.runVMHC:

            for num_strat, strat in enumerate(strat_list):

                create_vmhc(workflow, num_strat, strat, c,
                        output_name='vmhc_{0}'.format(num_strat))

        strat_list += new_strat_list

        # Inserting REHO Workflow

        if 1 in c.runReHo:

            for num_strat, strat in enumerate(strat_list):

                preproc = create_reho()
                cluster_size = c.clusterSize

                # TODO ASH schema validator
                # Check the cluster size is supported
                if cluster_size not in [7, 19, 27]:
                    err_msg = 'Cluster size specified: %d, is not supported. ' \
                            'Change to 7, 19, or 27 and try again' % cluster_size
                    raise Exception(err_msg)
                else:
                    preproc.inputs.inputspec.cluster_size = cluster_size
                    reho = preproc.clone('reho_%d' % num_strat)

                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                reho, 'inputspec.rest_res_filt')

                node, out_file = strat['functional_brain_mask']
                workflow.connect(node, out_file,
                                reho, 'inputspec.rest_mask')

                strat.update_resource_pool({
                    'reho': (reho, 'outputspec.raw_reho_map')
                })

        ts_analysis_dict = {}
        sca_analysis_dict = {}

        # TODO ASH normalize w schema val
        if c.tsa_roi_paths:

            tsa_roi_dict = c.tsa_roi_paths[0]

            # Timeseries and SCA config selections processing

            # flip the dictionary
            for roi_path in tsa_roi_dict.keys():
                ts_analysis_to_run = [
                    x.strip() for x in tsa_roi_dict[roi_path].split(",")
                ]

                if any(
                    corr in ts_analysis_to_run for corr in [
                        "PearsonCorr", "PartialCorr"
                    ]
                ) and "Avg" not in ts_analysis_to_run:
                    ts_analysis_to_run += ["Avg"]

                for analysis_type in ts_analysis_to_run:
                    if analysis_type not in ts_analysis_dict.keys():
                        ts_analysis_dict[analysis_type] = []
                    ts_analysis_dict[analysis_type].append(roi_path)

        # TODO ASH normalize w schema val
        if 1 in c.runROITimeseries:

            # TODO ASH normalize w schema val
            if not c.tsa_roi_paths:
                err = "\n\n[!] CPAC says: Time Series Extraction is " \
                    "set to run, but no ROI NIFTI file paths were provided!" \
                    "\n\n"
                raise Exception(err)

        # TODO ASH normalize w schema val
        if 1 in c.runSCA:

            # TODO ASH normalize w schema val
            if c.sca_roi_paths:
                sca_roi_dict = c.sca_roi_paths[0]
            else:
                err = "\n\n[!] CPAC says: Seed-based Correlation Analysis is " \
                    "set to run, but no ROI NIFTI file paths were provided!" \
                    "\n\n"
                raise Exception(err)

            # flip the dictionary
            for roi_path in sca_roi_dict.keys():
                # update analysis dict
                for analysis_type in sca_roi_dict[roi_path].split(","):
                    analysis_type = analysis_type.replace(" ", "")

                    if analysis_type not in sca_analysis_dict.keys():
                        sca_analysis_dict[analysis_type] = []

                    sca_analysis_dict[analysis_type].append(roi_path)

        # Section: Spatial Regression Based Time Series

        new_strat_list = []

        if "SpatialReg" in ts_analysis_dict.keys(
        ) or "DualReg" in sca_analysis_dict.keys():

            for num_strat, strat in enumerate(strat_list):

                if "SpatialReg" in ts_analysis_dict.keys():

                    resample_spatial_map_to_native_space = pe.Node(
                        interface=fsl.FLIRT(),
                        name='resample_spatial_map_to_native_space_%d' % num_strat
                    )
                    resample_spatial_map_to_native_space.inputs.set(
                        interp='nearestneighbour',
                        apply_xfm=True,
                        in_matrix_file=c.identityMatrix
                    )

                    spatial_map_dataflow = create_spatial_map_dataflow(
                        ts_analysis_dict["SpatialReg"],
                        'spatial_map_dataflow_%d' % num_strat
                    )

                    spatial_map_dataflow.inputs.inputspec.set(
                        creds_path=input_creds_path,
                        dl_dir=c.workingDirectory
                    )

                    spatial_map_timeseries = get_spatial_map_timeseries(
                        'spatial_map_timeseries_%d' % num_strat
                    )
                    spatial_map_timeseries.inputs.inputspec.demean = True  # c.spatialDemean

                    node, out_file = strat['functional_to_standard']
                    node2, out_file2 = strat['functional_brain_mask_to_standard']

                    # resample the input functional file and functional mask
                    # to spatial map
                    workflow.connect(node, out_file,
                                    resample_spatial_map_to_native_space,
                                    'reference')
                    workflow.connect(spatial_map_dataflow,
                                    'select_spatial_map.out_file',
                                    resample_spatial_map_to_native_space,
                                    'in_file')

                    # connect it to the spatial_map_timeseries
                    workflow.connect(resample_spatial_map_to_native_space,
                                    'out_file',
                                    spatial_map_timeseries,
                                    'inputspec.spatial_map')
                    workflow.connect(node2, out_file2,
                                    spatial_map_timeseries,
                                    'inputspec.subject_mask')
                    workflow.connect(node, out_file,
                                    spatial_map_timeseries,
                                    'inputspec.subject_rest')

                    strat.append_name(spatial_map_timeseries.name)

                    strat.update_resource_pool({
                        'spatial_map_timeseries': (spatial_map_timeseries, 'outputspec.subject_timeseries')
                    })

                if "DualReg" in sca_analysis_dict.keys():
                    resample_spatial_map_to_native_space_for_dr = pe.Node(
                        interface=fsl.FLIRT(),
                        name='resample_spatial_map_to_native_space_for_DR_%d' % num_strat
                    )
                    resample_spatial_map_to_native_space_for_dr.inputs.set(
                        interp='nearestneighbour',
                        apply_xfm=True,
                        in_matrix_file=c.identityMatrix
                    )

                    spatial_map_dataflow_for_dr = create_spatial_map_dataflow(
                        sca_analysis_dict["DualReg"],
                        'spatial_map_dataflow_for_DR_%d' % num_strat
                    )

                    spatial_map_dataflow_for_dr.inputs.inputspec.set(
                        creds_path=input_creds_path,
                        dl_dir=c.workingDirectory
                    )

                    spatial_map_timeseries_for_dr = get_spatial_map_timeseries(
                        'spatial_map_timeseries_for_DR_%d' % num_strat
                    )

                    spatial_map_timeseries_for_dr.inputs.inputspec.demean = True  # c.spatialDemean

                    node, out_file = strat['functional_to_standard']
                    node2, out_file2 = strat['functional_brain_mask_to_standard']

                    # resample the input functional file and functional mask
                    # to spatial map
                    workflow.connect(node, out_file,
                                    resample_spatial_map_to_native_space_for_dr,
                                    'reference')
                    workflow.connect(spatial_map_dataflow_for_dr,
                                    'select_spatial_map.out_file',
                                    resample_spatial_map_to_native_space_for_dr,
                                    'in_file')

                    # connect it to the spatial_map_timeseries
                    workflow.connect(
                        resample_spatial_map_to_native_space_for_dr,
                        'out_file',
                        spatial_map_timeseries_for_dr,
                        'inputspec.spatial_map'
                    )

                    workflow.connect(node, out_file,
                                     spatial_map_timeseries_for_dr,
                                     'inputspec.subject_rest')
                    strat.append_name(spatial_map_timeseries_for_dr.name)

                    strat.update_resource_pool({
                        'spatial_map_timeseries_for_DR': (
                        spatial_map_timeseries_for_dr,
                        'outputspec.subject_timeseries')
                    })

        strat_list += new_strat_list

        if 1 in c.runROITimeseries and (
            "Avg" in ts_analysis_dict.keys() or
            "Avg" in sca_analysis_dict.keys() or
            "MultReg" in sca_analysis_dict.keys()
        ):

            # ROI Based Time Series
            new_strat_list = []

            for num_strat, strat in enumerate(strat_list):

                if "Avg" in ts_analysis_dict.keys():
                    resample_functional_roi = pe.Node(Function(input_names = ['in_func',
                                                                              'in_roi',
                                                                              'realignment',
                                                                              'identity_matrix'],
                                              output_names = ['out_func',
                                                              'out_roi'],
                                              function = resample_func_roi,
                                              as_module = True),
                                        name = 'resample_functional_roi_{0}'.format(num_strat))

                    resample_functional_roi.inputs.realignment = c.realignment
                    resample_functional_roi.inputs.identity_matrix = c.identityMatrix

                    roi_dataflow = create_roi_mask_dataflow(
                        ts_analysis_dict["Avg"],
                        'roi_dataflow_%d' % num_strat
                    )

                    roi_dataflow.inputs.inputspec.set(
                        creds_path=input_creds_path,
                        dl_dir=c.workingDirectory
                    )

                    roi_timeseries = get_roi_timeseries(
                        'roi_timeseries_%d' % num_strat
                    )
                    roi_timeseries.inputs.inputspec.output_type = c.roiTSOutputs

                    node, out_file = strat['functional_to_standard']

                    # resample the input functional file to roi
                    workflow.connect(node, out_file,
                                     resample_functional_roi, 'in_func')
                    workflow.connect(roi_dataflow, 'outputspec.out_file',
                                     resample_functional_roi, 'in_roi')

                    # connect it to the roi_timeseries
                    # workflow.connect(roi_dataflow, 'outputspec.out_file',
                    #                  roi_timeseries, 'input_roi.roi')
                    workflow.connect(resample_functional_roi, 'out_roi',
                                     roi_timeseries, 'input_roi.roi')
                    workflow.connect(resample_functional_roi, 'out_func',
                                     roi_timeseries, 'inputspec.rest')

                    strat.append_name(roi_timeseries.name)
                    strat.update_resource_pool({
                        'roi_timeseries': (roi_timeseries, 'outputspec.roi_outputs'),
                        'functional_to_roi': (resample_functional_roi, 'out_func')
                    })

                    # create the graphs
                    from CPAC.utils.ndmg_utils import ndmg_create_graphs

                    ndmg_graph = pe.MapNode(Function(
                        input_names=['ts', 'labels'],
                        output_names=['out_file'],
                        function=ndmg_create_graphs,
                        as_module=True
                    ), name='ndmg_graphs_%d' % num_strat,
                        iterfield=['labels'])

                    workflow.connect(roi_timeseries, 'outputspec.roi_ts', ndmg_graph, 'ts')
                    workflow.connect(roi_dataflow,
                                     'outputspec.out_file',
                                     ndmg_graph, 'labels')

                    strat.update_resource_pool({
                        'ndmg_graph': (ndmg_graph, 'out_file')
                    })

                if "Avg" in sca_analysis_dict.keys():

                    # same workflow, except to run TSE and send it to the resource
                    # pool so that it will not get sent to SCA
                    resample_functional_roi_for_sca = pe.Node(Function(input_names = ['in_func',
                                                                                      'in_roi',
                                                                                      'realignment',
                                                                                      'identity_matrix'],
                                              output_names = ['out_func',
                                                              'out_roi'],
                                              function = resample_func_roi,
                                              as_module = True),
                                        name = 'resample_functional_roi_for_sca_{0}'.format(num_strat))

                    resample_functional_roi_for_sca.inputs.realignment = c.realignment
                    resample_functional_roi_for_sca.inputs.identity_matrix = c.identityMatrix

                    roi_dataflow_for_sca = create_roi_mask_dataflow(
                        sca_analysis_dict["Avg"],
                        'roi_dataflow_for_sca_%d' % num_strat
                    )

                    roi_dataflow_for_sca.inputs.inputspec.set(
                        creds_path=input_creds_path,
                        dl_dir=c.workingDirectory
                    )

                    roi_timeseries_for_sca = get_roi_timeseries(
                        'roi_timeseries_for_sca_%d' % num_strat
                    )

                    node, out_file = strat['functional_to_standard']

                    # resample the input functional file to roi
                    workflow.connect(node, out_file,
                                     resample_functional_roi_for_sca, 'in_func')
                    workflow.connect(roi_dataflow_for_sca, 'outputspec.out_file',
                                     resample_functional_roi_for_sca, 'in_roi')

                    # connect it to the roi_timeseries
                    workflow.connect(resample_functional_roi_for_sca, 'out_roi',
                                     roi_timeseries_for_sca, 'input_roi.roi')
                    workflow.connect(resample_functional_roi_for_sca, 'out_func',
                                     roi_timeseries_for_sca, 'inputspec.rest')

                    strat.append_name(roi_timeseries_for_sca.name)
                    strat.update_resource_pool({
                        'roi_timeseries_for_SCA': (roi_timeseries_for_sca, 'outputspec.roi_outputs'),
                        'functional_to_roi_for_SCA': (resample_functional_roi, 'out_func')
                    })

                if "MultReg" in sca_analysis_dict.keys():

                    # same workflow, except to run TSE and send it to the resource
                    # pool so that it will not get sent to SCA
                    resample_functional_roi_for_multreg = pe.Node(Function(input_names = ['in_func',
                                                                                          'in_roi',
                                                                                          'realignment',
                                                                                          'identity_matrix'],
                                              output_names = ['out_func',
                                                              'out_roi'],
                                              function = resample_func_roi,
                                              as_module = True),
                                        name = 'resample_functional_roi_for_multreg_{0}'.format(num_strat))

                    resample_functional_roi_for_multreg.inputs.realignment = c.realignment
                    resample_functional_roi_for_multreg.inputs.identity_matrix = c.identityMatrix

                    roi_dataflow_for_multreg = create_roi_mask_dataflow(
                        sca_analysis_dict["MultReg"],
                        'roi_dataflow_for_mult_reg_%d' % num_strat
                    )

                    roi_dataflow_for_multreg.inputs.inputspec.set(
                        creds_path=input_creds_path,
                        dl_dir=c.workingDirectory
                    )

                    roi_timeseries_for_multreg = get_roi_timeseries(
                        'roi_timeseries_for_mult_reg_%d' % num_strat
                    )

                    node, out_file = strat['functional_to_standard']

                    # resample the input functional file to roi
                    workflow.connect(node, out_file,
                                    resample_functional_roi_for_multreg,
                                    'in_func')
                    workflow.connect(roi_dataflow_for_multreg,
                                    'outputspec.out_file',
                                    resample_functional_roi_for_multreg,
                                    'in_roi')

                    # connect it to the roi_timeseries
                    workflow.connect(resample_functional_roi_for_multreg,
                                    'out_roi',
                                    roi_timeseries_for_multreg,
                                    'input_roi.roi')
                    workflow.connect(resample_functional_roi_for_multreg,
                                    'out_func',
                                    roi_timeseries_for_multreg,
                                    'inputspec.rest')

                    strat.append_name(roi_timeseries_for_multreg.name)
                    strat.update_resource_pool({
                        'roi_timeseries_for_SCA_multreg': (roi_timeseries_for_multreg, 'outputspec.roi_outputs')
                    })

        strat_list += new_strat_list


        # Connectome
        if "PearsonCorr" in ts_analysis_dict.keys() or "PartialCorr" in ts_analysis_dict.keys():

            for num_strat, strat in enumerate(strat_list):

                if "PearsonCorr" in ts_analysis_dict.keys():
                    connectome_wf = create_connectome('connectome_PearsonCorr_%d' % num_strat)
                    connectome_wf.inputs.inputspec.method = "PearsonCorr"

                    node, out_file = strat['roi_timeseries']

                    workflow.connect(node,
                                    out_file,
                                    connectome_wf,
                                    'inputspec.time_series')

                    strat.update_resource_pool({
                        'connectome_PearsonCorr': (connectome_wf, 'outputspec.connectome')
                    })

                if "PartialCorr" in ts_analysis_dict.keys():
                    connectome_wf = create_connectome('connectome_PartialCorr_%d' % num_strat)
                    connectome_wf.inputs.inputspec.method = "PartialCorr"

                    node, out_file = strat['roi_timeseries']

                    workflow.connect(node,
                                    out_file,
                                    connectome_wf,
                                    'inputspec.time_series')

                    strat.update_resource_pool({
                        'connectome_PartialCorr': (connectome_wf, 'outputspec.connectome')
                    })

        # Voxel Based Time Series
        new_strat_list = []

        if "Voxel" in ts_analysis_dict.keys():

            for num_strat, strat in enumerate(strat_list):

                resample_functional_to_mask = pe.Node(Function(input_names = ['in_func',
                                                                            'in_roi',
                                                                            'realignment',
                                                                            'identity_matrix'],
                                              output_names = ['out_func',
                                                              'out_roi'],
                                              function = resample_func_roi,
                                              as_module = True),
                                        name = 'resample_functional_to_mask_{0}'.format(num_strat))

                resample_functional_to_mask.inputs.realignment = c.realignment
                resample_functional_to_mask.inputs.identity_matrix = c.identityMatrix

                mask_dataflow = create_roi_mask_dataflow(ts_analysis_dict["Voxel"],
                                                        'mask_dataflow_%d' % num_strat)

                voxel_timeseries = get_voxel_timeseries(
                    'voxel_timeseries_%d' % num_strat)
                voxel_timeseries.inputs.inputspec.output_type = c.roiTSOutputs

                node, out_file = strat['functional_to_standard']

                # resample the input functional file to mask
                workflow.connect(node, out_file,
                                resample_functional_to_mask, 'in_func')
                workflow.connect(mask_dataflow, 'outputspec.out_file',
                                resample_functional_to_mask, 'in_roi')

                # connect it to the voxel_timeseries
                workflow.connect(resample_functional_to_mask, 'out_roi',
                                voxel_timeseries, 'input_mask.mask')
                workflow.connect(resample_functional_to_mask, 'out_func',
                                voxel_timeseries, 'inputspec.rest')

                strat.append_name(voxel_timeseries.name)
                strat.update_resource_pool({
                    'voxel_timeseries': (voxel_timeseries, 'outputspec.mask_outputs')
                })

        strat_list += new_strat_list

        # Inserting SCA workflow for ROI INPUT

        new_strat_list = []

        if "Avg" in sca_analysis_dict.keys():

            for num_strat, strat in enumerate(strat_list):
                sca_roi = create_sca('sca_roi_%d' % num_strat)

                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                sca_roi, 'inputspec.functional_file')

                node, out_file = strat['roi_timeseries_for_SCA']
                workflow.connect(node, (out_file, extract_one_d),
                                sca_roi, 'inputspec.timeseries_one_d')

                strat.update_resource_pool({
                    'sca_roi_files': (sca_roi, 'outputspec.correlation_files')
                })

                strat.append_name(sca_roi.name)

        strat_list += new_strat_list

        # (Dual Regression) Temporal Regression for Dual Regression

        new_strat_list = []

        if "DualReg" in sca_analysis_dict.keys():

            for num_strat, strat in enumerate(strat_list):

                dr_temp_reg = create_temporal_reg(
                    'temporal_dual_regression_%d' % num_strat
                )
                dr_temp_reg.inputs.inputspec.normalize = c.mrsNorm
                dr_temp_reg.inputs.inputspec.demean = True

                node, out_file = strat['spatial_map_timeseries_for_DR']

                node2, out_file2 = strat.get_leaf_properties()
                node3, out_file3 = strat['functional_brain_mask']

                workflow.connect(node2, out_file2,
                                dr_temp_reg, 'inputspec.subject_rest')

                workflow.connect(node, out_file,
                                dr_temp_reg, 'inputspec.subject_timeseries')

                workflow.connect(node3, out_file3,
                                dr_temp_reg, 'inputspec.subject_mask')

                strat.update_resource_pool({
                    'dr_tempreg_maps_files': (dr_temp_reg, 'outputspec.temp_reg_map_files'),
                    'dr_tempreg_maps_zstat_files': (dr_temp_reg, 'outputspec.temp_reg_map_z_files')
                })

                strat.append_name(dr_temp_reg.name)

        strat_list += new_strat_list

        # (Multiple Regression) Temporal Regression for SCA

        new_strat_list = []

        if "MultReg" in sca_analysis_dict.keys():

            for num_strat, strat in enumerate(strat_list):

                sc_temp_reg = create_temporal_reg(
                    'temporal_regression_sca_%d' % num_strat,
                    which='RT'
                )
                sc_temp_reg.inputs.inputspec.normalize = c.mrsNorm
                sc_temp_reg.inputs.inputspec.demean = True

                node, out_file = strat['functional_to_standard']
                node2, out_file2 = strat['roi_timeseries_for_SCA_multreg']
                node3, out_file3 = strat['functional_brain_mask_to_standard']

                workflow.connect(node, out_file,
                                sc_temp_reg, 'inputspec.subject_rest')

                workflow.connect(node2, (out_file2, extract_one_d),
                                sc_temp_reg, 'inputspec.subject_timeseries')

                workflow.connect(node3, out_file3,
                                sc_temp_reg, 'inputspec.subject_mask')

                strat.update_resource_pool({
                    'sca_tempreg_maps_files': (sc_temp_reg, 'outputspec.temp_reg_map_files'),
                    'sca_tempreg_maps_zstat_files': (sc_temp_reg, 'outputspec.temp_reg_map_z_files')
                })

                strat.append_name(sc_temp_reg.name)

        strat_list += new_strat_list

        # Section: Network centrality

        # TODO ASH handle as boolean on schema validator / normalizer
        if 1 in c.runNetworkCentrality:

            # TODO ASH move to schema validator
            # validate the mask file path
            # if not c.templateSpecificationFile.endswith(".nii") and \
            #         not c.templateSpecificationFile.endswith(".nii.gz"):
            #     err = "\n\n[!] CPAC says: The Network Centrality mask " \
            #           "specification file must be a NIFTI file (ending in .nii " \
            #           "or .nii.gz).\nFile path you provided: %s\n\n" \
            #           % c.templateSpecificationFile

            #     raise Exception(err)

            strat_list = create_network_centrality_workflow(
                workflow, c, strat_list)

        '''
        Loop through the resource pool and connect the nodes for:
            - applying warps to standard
            - z-score standardization
            - smoothing
            - calculating output averages
        '''

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            if 'func_to_epi_fsl' in nodes or 'func_to_epi_ants' in nodes :

                rp = strat.get_resource_pool()

                for key in sorted(rp.keys()):

                    if key in Outputs.native_nonsmooth:
                        image_type = 'func_derivative'
                    elif key in Outputs.native_nonsmooth_mult:
                        image_type = 'func_derivative_multi'
                    else:
                        continue

                    output_name = '{0}_to_standard'.format(key)
                    if output_name not in strat:
                        output_func_to_standard(workflow, key, 'template_epi_derivative',
                            '{0}_to_standard'.format(key), strat, num_strat, c, input_image_type=image_type, registration_template='epi', func_type='non-ica-aroma')

            elif 'T1_template' in c.runRegisterFuncToTemplate:

                rp = strat.get_resource_pool()

                for key in sorted(rp.keys()):

                    if key in Outputs.native_nonsmooth:
                        image_type = 'func_derivative'
                    elif key in Outputs.native_nonsmooth_mult:
                        image_type = 'func_derivative_multi'
                    else:
                        continue

                    output_name = '{0}_to_standard'.format(key)
                    if output_name not in strat:
                        output_func_to_standard(workflow, key, 'template_brain_for_func_derivative',
                            '{0}_to_standard'.format(key), strat, num_strat, c, input_image_type=image_type, registration_template='t1', func_type='non-ica-aroma')

            if "Before" in c.smoothing_order:

                # run smoothing before Z-scoring
                if 1 in c.run_smoothing:
                    rp = strat.get_resource_pool()
                    for key in sorted(rp.keys()):
                        if 'centrality' in key or key in Outputs.native_nonsmooth + Outputs.native_nonsmooth_mult + \
                                Outputs.template_nonsmooth + Outputs.template_nonsmooth_mult:
                            spatial_smooth_outputs(workflow, key, strat, num_strat, c)
                            # c.smoothing_mehod can be FSL or AFNI, FSL as default

                if 1 in c.runZScoring:
                    rp = strat.get_resource_pool()

                    for key in sorted(rp.keys()):
                        # connect nodes for z-score standardization
                        if "sca_roi_files_to_standard" in key:
                            # correlation files need the r-to-z
                            strat = fisher_z_score_standardize(workflow, key,
                                                            "roi_timeseries_for_SCA",
                                                            strat, num_strat,
                                                            map_node=True)
                        elif "centrality" in key:
                            # specific mask
                            strat = z_score_standardize(workflow, key,
                                                        c.templateSpecificationFile,
                                                        strat, num_strat,
                                                        map_node=True)
                        elif key in Outputs.template_raw:
                            # raw score, in template space
                            strat = z_score_standardize(workflow, key,
                                                        "functional_brain_mask_to_standard_derivative",
                                                        strat, num_strat)

                        elif key in Outputs.template_raw_mult:
                            # same as above but multiple files so mapnode required
                            strat = z_score_standardize(workflow, key,
                                                        "functional_brain_mask_to_standard_derivative",
                                                        strat, num_strat,
                                                        map_node=True)

            elif "After" in c.smoothing_order:
                # run smoothing after Z-scoring
                if 1 in c.runZScoring:
                    rp = strat.get_resource_pool()
                    for key in sorted(rp.keys()):
                        # connect nodes for z-score standardization
                        if "sca_roi_files_to_standard" in key:
                            # correlation files need the r-to-z
                            strat = fisher_z_score_standardize(workflow, key,
                                                            "roi_timeseries_for_SCA",
                                                            strat, num_strat,
                                                            map_node=True)
                        elif "centrality" in key:
                            # specific mask
                            strat = z_score_standardize(workflow, key,
                                                        c.templateSpecificationFile,
                                                        strat, num_strat,
                                                        map_node=True)
                        elif key in Outputs.template_raw:
                            # raw score, in template space
                            strat = z_score_standardize(workflow, key,
                                                        "functional_brain_mask_to_standard_derivative",
                                                        strat, num_strat)
                        elif key in Outputs.template_raw_mult:
                            # same as above but multiple files so mapnode required
                            strat = z_score_standardize(workflow, key,
                                                        "functional_brain_mask_to_standard_derivative",
                                                        strat, num_strat,
                                                        map_node=True)

                if 1 in c.run_smoothing:

                    rp = strat.get_resource_pool()

                    for key in sorted(rp.keys()):
                        if 'centrality' in key or key in Outputs.native_nonsmooth + Outputs.native_nonsmooth_mult + \
                                Outputs.template_nonsmooth + Outputs.template_nonsmooth_mult:
                            spatial_smooth_outputs(workflow, key, strat, num_strat, c)

            rp = strat.get_resource_pool()
            for key in sorted(rp.keys()):
                # connect nodes to calculate averages
                if key in Outputs.average:
                    # the outputs we need the averages for
                    strat = calc_avg(workflow, key, strat, num_strat)
                elif key in Outputs.average_mult:
                    # those outputs, but the ones with multiple files (map nodes)
                    strat = calc_avg(workflow, key, strat,
                                    num_strat, map_node=True)

    # Quality Control
    if 1 in c.generateQualityControlImages:
        create_qc_workflow(workflow, c, strat_list, Outputs.qc)


    ndmg_out = False
    try:
        if "ndmg" in c.output_tree:
            ndmg_out = True
    except:
        pass


    # TODO enforce value with schema validation
    try:
        encrypt_data = bool(c.s3Encryption[0])
    except:
        encrypt_data = False


    # TODO enforce value with schema validation
    # Extract credentials path for output if it exists
    try:
        # Get path to creds file
        creds_path = ''
        if c.awsOutputBucketCredentials:
            creds_path = str(c.awsOutputBucketCredentials)
            creds_path = os.path.abspath(creds_path)

        if c.outputDirectory.lower().startswith('s3://'):
            # Test for s3 write access
            s3_write_access = \
                aws_utils.test_bucket_access(creds_path,
                                                c.outputDirectory)

            if not s3_write_access:
                raise Exception('Not able to write to bucket!')

    except Exception as e:
        if c.outputDirectory.lower().startswith('s3://'):
            err_msg = 'There was an error processing credentials or ' \
                        'accessing the S3 bucket. Check and try again.\n' \
                        'Error: %s' % e
            raise Exception(err_msg)


    # this section creates names for the different branched strategies.
    # it identifies where the pipeline has forked and then appends the
    # name of the forked nodes to the branch name in the output directory

    fork_points_labels = Strategy.get_forking_labels(strat_list)

    # DataSinks
    pipeline_ids = []

    scan_ids = ['scan_anat']
    if 'func' in sub_dict:
        scan_ids += ['scan_' + str(scan_id)
                        for scan_id in sub_dict['func']]
    if 'rest' in sub_dict:
        scan_ids += ['scan_' + str(scan_id)
                        for scan_id in sub_dict['rest']]


    for num_strat, strat in enumerate(strat_list):

        if pipeline_name is None or pipeline_name == 'None':
            pipeline_id = c.pipelineName
        else:
            pipeline_id = pipeline_name

        if fork_points_labels[strat]:
            pipeline_id += '_' + fork_points_labels[strat]

        pipeline_ids.append(pipeline_id)

        rp = strat.get_resource_pool()

        for resource_i, resource in enumerate(sorted(rp.keys())):

            if not resource.startswith('qc___') and resource not in Outputs.any:
                continue

            if resource not in Outputs.override_optional and not ndmg_out:

                if 1 not in c.write_func_outputs:
                    if resource in Outputs.extra_functional:
                        continue

                if 1 not in c.write_debugging_outputs:
                    if resource in Outputs.debugging:
                        continue

                if 'Off' not in c.runRegisterFuncToTemplate:
                    if resource in Outputs.native_nonsmooth or \
                        resource in Outputs.native_nonsmooth_mult or \
                            resource in Outputs.native_smooth:
                        continue

                if 0 not in c.runZScoring:
                    # write out only the z-scored outputs
                    if resource in Outputs.template_raw or \
                            resource in Outputs.template_raw_mult:
                        continue

                if 0 not in c.run_smoothing:
                    # write out only the smoothed outputs
                    if resource in Outputs.native_nonsmooth or \
                        resource in Outputs.template_nonsmooth or \
                            resource in Outputs.native_nonsmooth_mult or \
                            resource in Outputs.template_nonsmooth_mult:
                        continue

            if ndmg_out:
                ds = pe.Node(DataSink(),
                                name='sinker_{}_{}'.format(num_strat,
                                                        resource_i))
                ds.inputs.base_directory = c.outputDirectory
                ds.inputs.creds_path = creds_path
                ds.inputs.encrypt_bucket_keys = encrypt_data
                ds.inputs.parameterization = True
                ds.inputs.regexp_substitutions = [
                    (r'_rename_(.)*/', ''),
                    (r'_scan_', 'scan-'),
                    (r'/_mask_', '/roi-'),
                    (r'file_s3(.)*/', ''),
                    (r'ndmg_atlases', ''),
                    (r'func_atlases', ''),
                    (r'label', ''),
                    (r'res-.+\/', ''),
                    (r'_mask_', 'roi-'),
                    (r'mask_sub-', 'sub-'),
                    (r'/_selector_', '_nuis-'),
                    (r'_selector_pc', ''),
                    (r'.linear', ''),
                    (r'.wm', ''),
                    (r'.global', ''),
                    (r'.motion', ''),
                    (r'.quadratic', ''),
                    (r'.gm', ''),
                    (r'.compcor', ''),
                    (r'.csf', ''),
                    (r'_sub-', '/sub-'),
                    (r'(\.\.)', '')
                ]

                container = 'pipeline_{0}'.format(pipeline_id)

                sub_ses_id = subject_id.split('_')

                if 'sub-' not in sub_ses_id[0]:
                    sub_tag = 'sub-{0}'.format(sub_ses_id[0])
                else:
                    sub_tag = sub_ses_id[0]

                ses_tag = 'ses-1'
                if len(sub_ses_id) > 1:
                    if 'ses-' not in sub_ses_id[1]:
                        ses_tag = 'ses-{0}'.format(sub_ses_id[1])
                    else:
                        ses_tag = sub_ses_id[1]

                id_tag = '_'.join([sub_tag, ses_tag])

                anat_template_tag = 'standard'
                func_template_tag = 'standard'

                try:
                    if 'FSL' in c.regOption and 'ANTS' not in c.regOption:
                        if 'MNI152' in c.fnirtConfig:
                            anat_template_tag = 'MNI152'
                            func_template_tag = 'MNI152'
                except:
                    pass

                anat_res_tag = c.resolution_for_anat.replace('mm','')
                func_res_tag = c.resolution_for_func_preproc.replace('mm','')

                ndmg_key_dct = {
                    'anatomical_brain': (
                        'anat',
                        'preproc',
                        '{0}_T1w_preproc_brain'.format(id_tag)
                    ),
                    'anatomical_to_standard': (
                        'anat',
                        'registered',
                        '{0}_T1w_space-{1}_res-{2}x{2}x{2}_registered'
                        .format(id_tag, anat_template_tag, anat_res_tag)
                    ),
                    'functional_preprocessed': (
                        'func',
                        'preproc',
                        '{0}_bold_preproc'
                        .format(id_tag)
                    ),
                    'functional_nuisance_residuals': (
                        'func',
                        'clean',
                        '{0}_bold_space-{1}_res-{2}x{2}x{2}_clean'
                        .format(id_tag, func_template_tag, func_res_tag)
                    ),
                    'functional_to_standard': (
                        'func',
                        'registered',
                        '{0}_bold_space-{1}_res-{2}x{2}x{2}_registered'
                        .format(id_tag, func_template_tag, func_res_tag)
                    ),
                    'functional_brain_mask_to_standard': (
                        'func',
                        'registered',
                        '{0}_bold_space-{1}_res-{2}x{2}x{2}_registered_mask'
                        .format(id_tag, func_template_tag, func_res_tag)
                    ),
                    'roi_timeseries': (
                        'func',
                        'roi-timeseries',
                        '{0}_bold_res-{1}x{1}x{1}_variant-mean_timeseries'
                        .format(id_tag, func_res_tag)
                    ),
                    'ndmg_graph': (
                        'func',
                        'roi-connectomes',
                        '{0}_bold_res-{1}x{1}x{1}_measure-correlation'
                        .format(id_tag, func_res_tag)
                    )
                }

                if resource not in ndmg_key_dct.keys():
                    continue

                ds.inputs.container = '{0}/{1}'.format(container,
                                                        ndmg_key_dct[resource][0])
                node, out_file = rp[resource]

                # rename the file
                if 'roi_' in resource or 'ndmg_graph' in resource:
                    rename_file = pe.MapNode(
                        interface=util.Rename(),
                        name='rename__{}_{}'.format(num_strat, resource_i),
                        iterfield=['in_file']
                    )
                else:
                    rename_file = pe.Node(
                        interface=util.Rename(),
                        name='rename_{}_{}'.format(num_strat, resource_i)
                    )
                rename_file.inputs.keep_ext = True
                rename_file.inputs.format_string = ndmg_key_dct[resource][2]

                workflow.connect(node, out_file,
                                    rename_file, 'in_file')
                workflow.connect(rename_file, 'out_file',
                                    ds, ndmg_key_dct[resource][1])

            else:
                output_sink_nodes = []

                # regular datasink
                ds = pe.Node(
                    DataSink(),
                    name='sinker_{}_{}'.format(num_strat, resource)
                )
                ds.inputs.base_directory = c.outputDirectory
                ds.inputs.creds_path = creds_path
                ds.inputs.encrypt_bucket_keys = encrypt_data
                ds.inputs.container = os.path.join(
                    'pipeline_{0}'.format(pipeline_id), subject_id
                )
                ds.inputs.regexp_substitutions = [
                    (r"/_sca_roi(.)*[/]", '/'),
                    (r"/_smooth_centrality_(\d)+[/]", '/'),
                    (r"/_z_score(\d)+[/]", "/"),
                    (r"/_dr_tempreg_maps_zstat_files_smooth_(\d)+[/]", "/"),
                    (r"/_sca_tempreg_maps_zstat_files_smooth_(\d)+[/]", "/"),
                    (r"/qc___", '/qc/')
                ]

                output_sink_nodes = []
                node, out_file = rp[resource]

                # exclue Nonetype transforms
                if resource == 'ants_initial_xfm' or resource == 'ants_rigid_xfm' or resource == 'ants_affine_xfm' \
                    or resource == 'ants_symmetric_initial_xfm' or resource == 'ants_symmetric_rigid_xfm' or resource == 'ants_symmetric_affine_xfm':

                    ants_para = c.ANTs_para_T1_registration
                    for para_index in range(len(ants_para)):
                        for para_type in ants_para[para_index]:
                            if para_type == 'initial-moving-transform':
                                if ants_para[para_index][para_type]['initializationFeature'] == 0 and resource == 'ants_initial_xfm':
                                    workflow.connect(node, out_file, ds, resource)
                            elif para_type == 'transforms':
                                for trans_index in range(len(ants_para[para_index][para_type])):
                                    for trans_type in ants_para[para_index][para_type][trans_index]:
                                        if trans_type == 'Rigid' and resource == 'ants_rigid_xfm':
                                            workflow.connect(node, out_file, ds, resource)
                                        if trans_type == 'Affine' and resource == 'ants_affine_xfm':
                                            workflow.connect(node, out_file, ds, resource)
                # exclue Nonetype transforms
                if resource == 'func_to_epi_ants_initial_xfm' or resource == 'func_to_epi_ants_rigid_xfm' or resource == 'func_to_epi_ants_affine_xfm':
                    ants_para = c.ANTs_para_EPI_registration
                    for para_index in range(len(ants_para)):
                        for para_type in ants_para[para_index]:
                            if para_type == 'initial-moving-transform':
                                if ants_para[para_index][para_type]['initializationFeature'] == 0 and resource == 'func_to_epi_ants_initial_xfm':
                                    workflow.connect(node, out_file, ds, resource)
                            elif para_type == 'transforms':
                                for trans_index in range(len(ants_para[para_index][para_type])):
                                    for trans_type in ants_para[para_index][para_type][trans_index]:
                                        if trans_type == 'Rigid' and resource == 'func_to_epi_ants_rigid_xfm':
                                            workflow.connect(node, out_file, ds, resource)
                                        if trans_type == 'Affine' and resource == 'func_to_epi_ants_affine_xfm':
                                            workflow.connect(node, out_file, ds, resource)
                if resource not in ['ants_initial_xfm', 'ants_rigid_xfm', 'ants_affine_xfm', 'func_to_epi_ants_initial_xfm', 'func_to_epi_ants_rigid_xfm', 'func_to_epi_ants_affine_xfm',\
                    'ants_symmetric_initial_xfm','ants_symmetric_rigid_xfm','ants_symmetric_affine_xfm']:
                    workflow.connect(node, out_file, ds, resource)

                output_sink_nodes += [(ds, 'out_file')]

    logger.info('\n\n' + 'Pipeline building completed.' + '\n\n')

    return workflow, strat_list, pipeline_ids
