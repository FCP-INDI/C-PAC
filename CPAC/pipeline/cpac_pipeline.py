import os
import time
import csv
import shutil
import pickle
import copy
import yaml

import logging as cb_logging
from time import strftime

import nipype
from CPAC.pipeline import nipype_pipeline_engine as pe
from nipype import config
from nipype import logging

from indi_aws import aws_utils, fetch_creds

import CPAC

from CPAC.pipeline.engine import NodeBlock, initiate_rpool
from CPAC.anat_preproc.anat_preproc import (
    freesurfer_preproc,
    anatomical_init,
    acpc_align_head,
    acpc_align_head_with_mask,
    acpc_align_brain,
    acpc_align_brain_with_mask,
    non_local_means,
    n4_bias_correction,
    brain_mask_afni,
    brain_mask_fsl,
    brain_mask_niworkflows_ants,
    brain_mask_unet,
    brain_mask_freesurfer,
    brain_mask_acpc_afni,
    brain_mask_acpc_fsl,
    brain_mask_acpc_niworkflows_ants,
    brain_mask_acpc_unet,
    brain_mask_acpc_freesurfer,
    brain_extraction
)

from CPAC.registration.registration import (
    register_ANTs_anat_to_template,
    register_FSL_anat_to_template,
    register_symmetric_ANTs_anat_to_template,
    register_symmetric_FSL_anat_to_template,
    register_ANTs_EPI_to_template,
    register_FSL_EPI_to_template,
    coregistration_prep_vol,
    coregistration_prep_mean,
    coregistration,
    create_func_to_T1template_xfm,
    create_func_to_T1template_symmetric_xfm,
    warp_timeseries_to_T1template,
    warp_bold_mean_to_T1template,
    warp_bold_mask_to_T1template,
    warp_deriv_mask_to_T1template,
    warp_timeseries_to_EPItemplate,
    warp_bold_mean_to_EPItemplate,
    warp_bold_mask_to_EPItemplate,
    warp_deriv_mask_to_EPItemplate
)

from CPAC.seg_preproc.seg_preproc import (
    tissue_seg_fsl_fast,
    tissue_seg_T1_template_based,
    tissue_seg_EPI_template_based,
    tissue_seg_ants_prior,
    tissue_seg_freesurfer
)

from CPAC.func_preproc.func_preproc import (
    func_scaling,
    func_truncate,
    func_despike,
    func_slice_time,
    func_reorient,
    bold_mask_afni,
    bold_mask_fsl,
    bold_mask_fsl_afni,
    bold_mask_anatomical_refined,
    bold_mask_anatomical_based,
    bold_masking,
    func_mean,
    func_normalize,
    get_motion_ref,
    func_motion_estimates,
    motion_estimate_filter,
    calc_motion_stats,
    func_motion_correct_only,
    func_motion_correct
)

from CPAC.distortion_correction.distortion_correction import (
    distcor_phasediff_fsl_fugue,
    distcor_blip_afni_qwarp
)

from CPAC.nuisance.nuisance import (
    ICA_AROMA_ANTsreg,
    ICA_AROMA_FSLreg,
    ICA_AROMA_EPIreg,
    nuisance_regression_complete,
    erode_mask_T1w,
    erode_mask_CSF,
    erode_mask_GM,
    erode_mask_WM,
    nuisance_regression_EPItemplate,
    erode_mask_bold,
    erode_mask_boldCSF,
    erode_mask_boldGM,
    erode_mask_boldWM
)

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

from CPAC.alff.alff import alff_falff
from CPAC.reho.reho import reho

from CPAC.vmhc.vmhc import (
    smooth_func_vmhc,
    warp_timeseries_to_sym_template,
    vmhc
)

from CPAC.network_centrality.pipeline import (
    network_centrality
)

from CPAC.utils.datasource import (
    gather_extraction_maps
)
from CPAC.pipeline.schema import valid_options
from CPAC.utils.trimmer import the_trimmer
from CPAC.utils import Configuration

from CPAC.qc.pipeline import create_qc_workflow
from CPAC.qc.utils import generate_qc_pages

from CPAC.utils.utils import (
    check_config_resources,
    check_system_deps,
)

from CPAC.utils.monitoring import log_nodes_cb, log_nodes_initial
from CPAC.utils.monitoring.draw_gantt_chart import resource_report

logger = logging.getLogger('nipype.workflow')

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

    # Assure that changes on config will not affect other parts
    c = copy.copy(c)

    subject_id = sub_dict['subject_id']
    if sub_dict['unique_id']:
        subject_id += "_" + sub_dict['unique_id']

    log_dir = os.path.join(c.pipeline_setup['log_directory']['path'],
                           f'pipeline_{c.pipeline_setup["pipeline_name"]}',
                           subject_id)
    if not os.path.exists(log_dir):
        os.makedirs(os.path.join(log_dir))

    # TODO ASH Enforce c.run_logging to be boolean
    # TODO ASH Schema validation
    config.update_config({
        'logging': {
            'log_directory': log_dir,
            'log_to_file': bool(getattr(c.pipeline_setup['log_directory'],
                                        'run_logging', True))
        },
        'execution': {
            'crashfile_format': 'txt'
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
    sub_mem_gb, num_cores_per_sub, num_ants_cores, num_omp_cores = check_config_resources(
        c)

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

    C-PAC version: {cpac_version}

    Setting maximum number of cores per participant to {cores}
    Setting number of participants at once to {participants}
    Setting OMP_NUM_THREADS to {omp_threads}
    Setting MKL_NUM_THREADS to 1
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
        cores=c.pipeline_setup['system_config']['max_cores_per_participant'],
        participants=c.pipeline_setup['system_config'][
            'num_participants_at_once'],
        omp_threads=c.pipeline_setup['system_config']['num_OMP_threads'],
        ants_threads=c.pipeline_setup['system_config']['num_ants_threads'],
        max_cores=max_core_usage
    ))

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
    c.pipeline_setup['working_directory']['path'] = os.path.abspath(
        c.pipeline_setup['working_directory']['path'])
    if 's3://' not in c.pipeline_setup['output_directory']['path']:
        c.pipeline_setup['output_directory']['path'] = os.path.abspath(
            c.pipeline_setup['output_directory']['path'])

    workflow = build_workflow(
        subject_id, sub_dict, c, p_name, num_ants_cores
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
                output_dir=c.pipeline_setup['output_directory']['path'],
                s3_creds_path=input_creds_path,
            )

        pipeline_start_datetime = strftime("%Y-%m-%d %H:%M:%S")

        try:
            subject_info['resource_pool'] = []

            # for strat_no, strat in enumerate(strat_list):
            #    strat_label = 'strat_%d' % strat_no
            #    subject_info[strat_label] = strat.get_name()
            #    subject_info['resource_pool'].append(strat.get_resource_pool())

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
            if nipype.__version__ not in ('1.5.1'):
                err_msg = "This version of Nipype may not be compatible with " \
                          "CPAC v%s, please install Nipype version 1.5.1\n" \
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

                        if headerExists == False:
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
                    err_msg = 'Unable to upload CPAC log files in: %s.\nError: %s'
                    logger.error(err_msg, log_dir, exc)

        except Exception as e:
            import traceback;
            traceback.print_exc()
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

            if workflow:

                resource_report(cb_log_filename,
                                num_cores_per_sub, logger)

                logger.info(execution_info.format(
                    workflow=workflow.name,
                    pipeline=c.pipeline_setup['pipeline_name'],
                    log_dir=c.pipeline_setup['log_directory']['path'],
                    elapsed=(time.time() - pipeline_start_time) / 60,
                    run_start=pipeline_start_datetime,
                    run_finish=strftime("%Y-%m-%d %H:%M:%S")
                ))

                # Remove working directory when done
                if c.pipeline_setup['working_directory'][
                    'remove_working_dir']:
                    try:
                        if os.path.exists(working_dir):
                            logger.info("Removing working dir: %s",
                                        working_dir)
                            shutil.rmtree(working_dir)
                    except (FileNotFoundError, PermissionError):
                        logger.warn('Could not remove working directory %s',
                                    working_dir)


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
    if not rpool.check_rpool('desc-reorient_T1w'):
        anat_init_blocks = [
            anatomical_init
        ]
        pipeline_blocks += anat_init_blocks

    pipeline_blocks += [freesurfer_preproc]

    if not rpool.check_rpool('desc-preproc_T1w'):

        # brain masking for ACPC alignment
        if cfg.anatomical_preproc['acpc_alignment']['acpc_target'] == 'brain':
            if rpool.check_rpool('space-T1w_desc-brain_mask') or \
                    cfg.surface_analysis['run_freesurfer']:
                acpc_blocks = [
                    brain_extraction,
                    acpc_align_brain_with_mask
                    # outputs space-T1w_desc-brain_mask for later - keep the mask (the user provided)
                ]
            else:
                acpc_blocks = [
                    [brain_mask_acpc_afni,
                     brain_mask_acpc_fsl,
                     brain_mask_acpc_niworkflows_ants,
                     brain_mask_acpc_unet],
                       #brain_mask_acpc_freesurfer
                    # we don't want these masks to be used later
                    brain_extraction,
                    acpc_align_brain
                ]
        elif cfg.anatomical_preproc['acpc_alignment'][
            'acpc_target'] == 'whole-head':
            if rpool.check_rpool('space-T1w_desc-brain_mask') or \
                    cfg.surface_analysis['run_freesurfer']:
                acpc_blocks = [
                    acpc_align_head_with_mask
                    # outputs space-T1w_desc-brain_mask for later - keep the mask (the user provided)
                ]
            else:
                acpc_blocks = [
                    acpc_align_head  # does not output nor generate a mask
                ]

        anat_preproc_blocks = [
            non_local_means,
            n4_bias_correction
        ]
        if cfg.anatomical_preproc['acpc_alignment']['run_before_preproc']:
            anat_blocks = acpc_blocks + anat_preproc_blocks
        else:
            anat_blocks = anat_preproc_blocks + acpc_blocks

        pipeline_blocks += anat_blocks

    # Anatomical brain masking
    if not rpool.check_rpool('space-T1w_desc-brain_mask') or \
            cfg.surface_analysis['run_freesurfer']:
        anat_brain_mask_blocks = [
            [brain_mask_afni,
             brain_mask_fsl,
             brain_mask_niworkflows_ants,
             brain_mask_unet]
               #brain_mask_freesurfer
        ]
        pipeline_blocks += anat_brain_mask_blocks

    if not rpool.check_rpool('desc-brain_T1w'):
        anat_brain_blocks = [
            brain_extraction
        ]
        pipeline_blocks += anat_brain_blocks

    return pipeline_blocks


def build_T1w_registration_stack(rpool, cfg, pipeline_blocks=None):

    if not pipeline_blocks:
        pipeline_blocks = []

    reg_blocks = []
    if not rpool.check_rpool('space-template_desc-brain_T1w'):
        reg_blocks = [
            [register_ANTs_anat_to_template, register_FSL_anat_to_template]
        ]
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
             tissue_seg_ants_prior]
             #tissue_seg_freesurfer
        ]
        if 'T1_Template' in cfg.segmentation['tissue_segmentation'][
            'Template_Based']['template_for_segmentation']:
            seg_blocks = [
                [tissue_seg_fsl_fast,
                 tissue_seg_ants_prior,
                 tissue_seg_T1_template_based]
                # tissue_seg_freesurfer
            ]
        if 'EPI_Template' in cfg.segmentation['tissue_segmentation'][
            'Template_Based']['template_for_segmentation']:
            seg_blocks = [
                [tissue_seg_fsl_fast,
                 tissue_seg_ants_prior,
                 tissue_seg_EPI_template_based]
                # tissue_seg_freesurfer
            ]
        pipeline_blocks += seg_blocks

    return pipeline_blocks


def connect_pipeline(wf, cfg, rpool, pipeline_blocks):

    for block in pipeline_blocks:
        try:
            nb = NodeBlock(block)
            wf = nb.connect_block(wf, cfg, rpool)
        except LookupError as e:
            previous_nb_str = (
                f"after node block '{previous_nb.get_name()}': "
            ) if previous_nb else 'at beginning:'
            # Alert user to block that raises error
            e.args = (
                'When trying to connect node block '
                f"'{NodeBlock(block).get_name()}' "
                f"to workflow '{wf}' " + previous_nb_str + e.args[0],
            )
            raise
        previous_nb = nb

    return wf


def build_workflow(subject_id, sub_dict, cfg, pipeline_name=None,
                   num_ants_cores=1):

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
    if cfg.functional_preproc['run'] and \
            (not rpool.check_rpool('desc-brain_bold') or
             not rpool.check_rpool('space-bold_desc-brain_mask') or
             not rpool.check_rpool('movement-parameters')):
        func_init_blocks = [
            func_scaling,
            func_truncate
        ]
        func_preproc_blocks = [
            func_despike,
            func_slice_time,
            func_reorient
        ]
        func_prep_blocks = [
            [bold_mask_afni, bold_mask_fsl, bold_mask_fsl_afni,
             bold_mask_anatomical_refined, bold_mask_anatomical_based],
            bold_masking,
            calc_motion_stats,
            func_mean,
            func_normalize
        ]

        # Distortion/Susceptibility Correction
        distcor_blocks = []
        if rpool.check_rpool('diffphase') and rpool.check_rpool('diffmag'):
            distcor_blocks.append(distcor_phasediff_fsl_fugue)

        if rpool.check_rpool('epi_1'):
            distcor_blocks.append(distcor_blip_afni_qwarp)

        if distcor_blocks:
            if len(distcor_blocks) > 1:
                distcor_blocks = [distcor_blocks]
            func_prep_blocks += distcor_blocks

        if cfg['functional_preproc']['motion_estimates_and_correction'][
            'calculate_motion_first']:
            func_motion_blocks = [
                get_motion_ref,
                func_motion_estimates,
                motion_estimate_filter
            ]
            func_blocks = func_init_blocks + func_motion_blocks + \
                          func_preproc_blocks + [func_motion_correct_only] + \
                          func_prep_blocks
        else:
            func_motion_blocks = [
                get_motion_ref,
                func_motion_correct,
                motion_estimate_filter
            ]
            func_blocks = func_init_blocks + func_preproc_blocks + \
                          func_motion_blocks + func_prep_blocks

        pipeline_blocks += func_blocks

    # BOLD to T1 coregistration
    if cfg.registration_workflows['functional_registration'][
        'coregistration']['run'] and \
            (not rpool.check_rpool('space-T1w_desc-mean_bold') or
             not rpool.check_rpool('from-bold_to-T1w_mode-image_desc-linear_xfm')):
        coreg_blocks = [
            [coregistration_prep_vol, coregistration_prep_mean],
            coregistration
        ]
        pipeline_blocks += coreg_blocks

    # BOLD to EPI-template registration (no T1w involved)
    if not rpool.check_rpool('space-template_desc-brain_bold'):
        if coregistration not in pipeline_blocks:
            pipeline_blocks += [coregistration_prep_vol, coregistration_prep_mean]
        EPI_reg_blocks = [
            [register_ANTs_EPI_to_template, register_FSL_EPI_to_template]
        ]
        pipeline_blocks += EPI_reg_blocks
        
    if 'EPI_Template' in cfg.segmentation['tissue_segmentation'][
        'Template_Based']['template_for_segmentation']:
        if not rpool.check_rpool('space-bold_label-CSF_mask') or \
                not rpool.check_rpool('space-bold_label-WM_mask'):
            pipeline_blocks += [tissue_seg_EPI_template_based]

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
    if not rpool.check_rpool('desc-cleaned_bold'):
        nuisance = []
        if 'T1_template' in \
                cfg.registration_workflows['functional_registration'][
                    'func_registration_to_template']['target_template'][
                    'using']:
            nuisance += [ICA_AROMA_ANTsreg, ICA_AROMA_FSLreg]
        if 'EPI_template' in \
                cfg.registration_workflows['functional_registration'][
                    'func_registration_to_template']['target_template'][
                    'using']:
            nuisance += [ICA_AROMA_EPIreg]

        if cfg.nuisance_corrections['2-nuisance_regression'][
                'Regressors']:
            if 'T1_template' in \
                cfg.registration_workflows['functional_registration'][
                    'func_registration_to_template']['target_template'][
                    'using']:
                nuisance_blocks = [
                    erode_mask_T1w,
                    erode_mask_CSF,
                    erode_mask_GM,
                    erode_mask_WM,
                    nuisance_regression_complete
                ]
                nuisance += nuisance_blocks
                
            if 'EPI_template' in \
                cfg.registration_workflows['functional_registration'][
                    'func_registration_to_template']['target_template'][
                    'using']:
                epi_nuisance_blocks = [
                    erode_mask_bold,
                    erode_mask_boldCSF,
                    erode_mask_boldGM,
                    erode_mask_boldWM,
                    nuisance_regression_EPItemplate
                ]
                nuisance += epi_nuisance_blocks

        pipeline_blocks += nuisance

    # Warp the functional time series to template space
    apply_func_warp = cfg.registration_workflows['functional_registration'][
        'coregistration']['run'] or 'EPI_template' in cfg.registration_workflows[
        'functional_registration']['func_registration_to_template'][
        'target_template']['using']
    template_funcs = [
        'space-template_desc-cleaned_bold',
        'space-template_desc-brain_bold',
        'space-template_desc-motion_bold',
        'space-template_desc-preproc_bold',
        'space-template_bold'
    ]
    for func in template_funcs:
        if rpool.check_rpool(func):
            apply_func_warp = False

    if apply_func_warp:
        pipeline_blocks += [[warp_timeseries_to_T1template,
                             warp_timeseries_to_EPItemplate],
                            warp_bold_mean_to_T1template,
                            warp_bold_mean_to_EPItemplate]
                            
    if not rpool.check_rpool('space-template_desc-bold_mask'):
        pipeline_blocks += [warp_bold_mask_to_T1template,
                            warp_deriv_mask_to_T1template,
                            warp_bold_mask_to_EPItemplate,
                            warp_deriv_mask_to_EPItemplate]

    # Extractions and Derivatives
    tse_atlases, sca_atlases = gather_extraction_maps(cfg)
    cfg.timeseries_extraction['tse_atlases'] = tse_atlases
    cfg.seed_based_correlation_analysis['sca_atlases'] = sca_atlases

    if not rpool.check_rpool('desc-Mean_timeseries') and \
                    'Avg' in tse_atlases:
        pipeline_blocks += [timeseries_extraction_AVG]

    if not rpool.check_rpool('desc-Voxel_timeseries') and \
                    'Voxel' in tse_atlases:
        pipeline_blocks += [timeseries_extraction_Voxel]

    if not rpool.check_rpool('desc-SpatReg_timeseries') and \
                    'SpatialReg' in tse_atlases:
        pipeline_blocks += [spatial_regression]

    if not rpool.check_rpool('desc-MeanSCA_correlations') and \
                    'Avg' in sca_atlases:
        pipeline_blocks += [SCA_AVG]

    if not rpool.check_rpool('desc-DualReg_correlations') and \
                    'DualReg' in sca_atlases:
        pipeline_blocks += [dual_regression]

    if not rpool.check_rpool('desc-MultReg_correlations') and \
                    'MultReg' in sca_atlases:
        pipeline_blocks += [multiple_regression]

    if not rpool.check_rpool('alff'):
        pipeline_blocks += [alff_falff]

    if not rpool.check_rpool('reho'):
        pipeline_blocks += [reho]

    if not rpool.check_rpool('vmhc'):
        pipeline_blocks += [smooth_func_vmhc,
                            warp_timeseries_to_sym_template,
                            vmhc]

    if not rpool.check_rpool('centrality') and \
            any([cfg.network_centrality[option]['weight_options'] for option in valid_options['centrality']['method_options']]):
        pipeline_blocks += [network_centrality]

    if cfg.pipeline_setup['output_directory'][
        'generate_quality_control_images']:
        qc_stack, qc_montage_id_a, qc_montage_id_s, qc_hist_id, qc_plot_id = \
            create_qc_workflow(cfg)
        pipeline_blocks += qc_stack

    # Connect the entire pipeline!
    wf = connect_pipeline(wf, cfg, rpool, pipeline_blocks)

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
