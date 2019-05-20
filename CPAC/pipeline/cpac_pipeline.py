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
from CPAC.EPI_DistCorr.EPI_DistCorr import create_EPI_DistCorr
from CPAC.func_preproc.func_preproc import (
    create_func_preproc,
    create_wf_edit_func
)
from CPAC.seg_preproc.seg_preproc import create_seg_preproc

from CPAC.warp.pipeline import (
    output_to_standard,
    z_score_standardize,
    fisher_z_score_standardize,
    output_smooth,
    calc_avg
)

from CPAC.registration import (
    create_fsl_flirt_linear_reg,
    create_fsl_fnirt_nonlinear_reg,
    create_register_func_to_anat,
    create_bbregister_func_to_anat,
    create_wf_calculate_ants_warp,
    create_wf_apply_ants_warp,
    create_wf_c3d_fsl_to_itk,
    create_wf_collect_transforms
)
from CPAC.nuisance import create_nuisance_workflow, bandpass_voxels, NuisanceRegressor
from CPAC.aroma import create_aroma
from CPAC.median_angle import create_median_angle_correction
from CPAC.generate_motion_statistics import motion_power_statistics
from CPAC.scrubbing import create_scrubbing_preproc
from CPAC.timeseries import (
    create_surface_registration,
    get_roi_timeseries,
    get_voxel_timeseries,
    get_vertices_timeseries,
    get_spatial_map_timeseries
)
from CPAC.warp.pipeline import (
    ants_apply_warps_func_mni,
    ants_apply_inverse_warps_template_to_func
)

from CPAC.vmhc.vmhc import create_vmhc
from CPAC.reho.reho import create_reho
from CPAC.alff.alff import create_alff
from CPAC.sca.sca import create_sca, create_temporal_reg

from CPAC.connectome.pipeline import create_connectome

from CPAC.utils.symlinks import create_symlinks
from CPAC.utils.datasource import (
    create_func_datasource,
    create_anat_datasource,
    create_roi_mask_dataflow,
    create_spatial_map_dataflow,
    create_check_for_s3_node
)
from CPAC.utils import Configuration, Strategy, Outputs, function, find_files

from CPAC.utils.interfaces.datasink import DataSink

from CPAC.qc.pipeline import create_qc_workflow
from CPAC.qc.utils import generate_qc_pages

from CPAC.utils.utils import (
    extract_one_d,
    get_scan_params,
    get_tr,
    extract_txt,
    create_log,
    extract_output_mean,
    create_output_mean_csv,
    get_zscore,
    get_fisher_zscore,
    add_afni_prefix
)

logger = logging.getLogger('nipype.workflow')
# config.enable_debug_mode()

# TODO ASH move to somewhere else
def pick_wm(seg_prob_list):
    seg_prob_list.sort()
    return seg_prob_list[-1]


def create_log_node(workflow, logged_wf, output, index, scan_id=None):
    try:
        log_dir = workflow.config['logging']['log_directory']
        if logged_wf:
            log_wf = create_log(wf_name='log_%s' % logged_wf.name)
            log_wf.inputs.inputspec.workflow = logged_wf.name
            log_wf.inputs.inputspec.index = index
            log_wf.inputs.inputspec.log_dir = log_dir
            workflow.connect(logged_wf, output, log_wf, 'inputspec.inputs')
        else:
            log_wf = create_log(wf_name='log_done_%s' % scan_id,
                                scan_id=scan_id)
            log_wf.base_dir = log_dir
            log_wf.inputs.inputspec.workflow = 'DONE'
            log_wf.inputs.inputspec.index = index
            log_wf.inputs.inputspec.log_dir = log_dir
            log_wf.inputs.inputspec.inputs = log_dir
            return log_wf
    except Exception as e:
        print(e)


def prep_workflow(sub_dict, c, run, pipeline_timing_info=None,
                  p_name=None, plugin='MultiProc', plugin_args=None):
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

    # calculate maximum potential use of cores according to current pipeline
    # configuration
    max_core_usage = int(c.maxCoresPerParticipant) * \
        int(c.numParticipantsAtOnce)

    information = """

    C-PAC version: {cpac_version}

    Setting maximum number of cores per participant to {cores}
    Setting number of participants at once to {participants}
    Setting OMP_NUM_THREADS to {threads}
    Setting MKL_NUM_THREADS to {threads}
    Setting ANTS/ITK thread usage to {ants_threads}
    Maximum potential number of cores that might be used during this run: {max_cores}

"""

    logger.info(information.format(
        cpac_version=CPAC.__version__,
        cores=c.maxCoresPerParticipant,
        participants=c.numParticipantsAtOnce,
        threads=numThreads,
        ants_threads=c.num_ants_threads,
        max_cores=max_core_usage
    ))

    # TODO ASH temporary code, remove
    # TODO ASH maybe scheme validation/normalization
    already_skullstripped = c.already_skullstripped[0]
    if already_skullstripped == 2:
        already_skullstripped = 0
    elif already_skullstripped == 3:
        already_skullstripped = 1

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

    # TODO ASH normalize file paths with schema validator
    template_anat_keys = [
        ("anat", "template_brain_only_for_anat"),
        ("anat", "template_skull_for_anat"),
        ("anat", "ref_mask"),
        ("anat", "template_symmetric_brain_only"),
        ("anat", "template_symmetric_skull"),
        ("anat", "dilated_symmetric_brain_mask"),
        ("anat", "templateSpecificationFile"),
        ("anat", "lateral_ventricles_mask"),
        ("anat", "PRIORS_CSF"),
        ("anat", "PRIORS_GRAY"),
        ("anat", "PRIORS_WHITE"),
        ("other", "configFileTwomm"),
    ]

    for key_type, key in template_anat_keys:

        node = create_check_for_s3_node(
            key,
            getattr(c, key), key_type,
            input_creds_path, c.workingDirectory
        )

        setattr(c, key, node)

    if c.reGenerateOutputs is True:
        working_dir = os.path.join(c.workingDirectory, workflow_name)
        erasable = list(find_files(working_dir, '*sink*')) + \
            list(find_files(working_dir, '*link*')) + \
            list(find_files(working_dir, '*log*'))

        for f in erasable:
            if os.path.isfile(f):
                os.remove(f)
            else:
                shutil.rmtree(f)

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
        lesion_datasource = create_anat_datasource(
            'lesion_gather_%d' % num_strat)
        lesion_datasource.inputs.inputnode.subject = subject_id
        lesion_datasource.inputs.inputnode.anat = sub_dict['lesion_mask']
        lesion_datasource.inputs.inputnode.creds_path = input_creds_path
        lesion_datasource.inputs.inputnode.dl_dir = c.workingDirectory

        strat_initial.update_resource_pool({
            'lesion_mask': (lesion_datasource, 'outputspec.anat')
        })

    strat_list += [strat_initial]

    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        if 'anatomical_brain_mask' in strat:

            anat_preproc = create_anat_preproc(method='mask',
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

            new_strat_list += [new_strat]

            continue

        if already_skullstripped:

            anat_preproc = create_anat_preproc(method=None,
                                               already_skullstripped=True,
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
            })

            new_strat_list += [new_strat]

        else:

            if not any(o in c.skullstrip_option for o in ["AFNI", "BET"]):
                err = '\n\n[!] C-PAC says: Your skull-stripping method options ' \
                    'setting does not include either \'AFNI\' or \'BET\'.\n\n' \
                    'Options you provided:\nskullstrip_option: {0}' \
                    '\n\n'.format(str(c.skullstrip_option))
                raise Exception(err)

            if "AFNI" in c.skullstrip_option:

                anat_preproc = create_anat_preproc(method='afni',
                                                   wf_name='anat_preproc_afni_%d' % num_strat)

                anat_preproc.inputs.AFNI_options.set(
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
                })

                new_strat_list += [new_strat]

            if "BET" in c.skullstrip_option:
                anat_preproc = create_anat_preproc(method='fsl',
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

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file,
                             flirt_reg_anat_mni, 'inputspec.input_brain')

            # pass the reference files
            workflow.connect(
                c.template_brain_only_for_anat, 'local_path',
                flirt_reg_anat_mni, 'inputspec.reference_brain'
            )

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

            create_log_node(workflow, flirt_reg_anat_mni, 'outputspec.output_brain',
                            num_strat)

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
                workflow.connect(
                    c.template_brain_only_for_anat, 'local_path',
                    fnirt_reg_anat_mni, 'inputspec.reference_brain'
                )

                node, out_file = strat['anatomical_reorient']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.input_skull')

                node, out_file = strat['anatomical_to_mni_linear_xfm']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.linear_aff')

                workflow.connect(
                    c.template_skull_for_anat, 'local_path',
                    fnirt_reg_anat_mni, 'inputspec.reference_skull'
                )

                workflow.connect(
                    c.ref_mask, 'local_path',
                    fnirt_reg_anat_mni, 'inputspec.ref_mask'
                )

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

                create_log_node(workflow, fnirt_reg_anat_mni, 'outputspec.output_brain',
                                num_strat)

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
                    num_threads=num_ants_cores
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

                # get the skull-stripped anatomical from resource pool
                node, out_file = strat['anatomical_brain']

                # pass the anatomical to the workflow
                workflow.connect(node, out_file,
                                 ants_reg_anat_mni,
                                 'inputspec.anatomical_brain')

                # pass the reference file
                workflow.connect(
                    c.template_brain_only_for_anat, 'local_path',
                    ants_reg_anat_mni, 'inputspec.reference_brain'
                )

                # get the reorient skull-on anatomical from resource pool
                node, out_file = strat['anatomical_reorient']

                # pass the anatomical to the workflow
                workflow.connect(node, out_file,
                                 ants_reg_anat_mni,
                                 'inputspec.anatomical_skull')

                # pass the reference file
                workflow.connect(
                    c.template_skull_for_anat, 'local_path',
                    ants_reg_anat_mni, 'inputspec.reference_skull'
                )

            else:
                
                node, out_file = strat['anatomical_brain']

                workflow.connect(node, out_file, ants_reg_anat_mni,
                                 'inputspec.anatomical_brain')

                # pass the reference file
                workflow.connect(
                    c.template_brain_only_for_anat, 'local_path',
                    ants_reg_anat_mni, 'inputspec.reference_brain'
                )

            ants_reg_anat_mni.inputs.inputspec.set(
                dimension=3,
                use_histogram_matching=True,
                winsorize_lower_quantile=0.01,
                winsorize_upper_quantile=0.99,
                metric=['MI', 'MI', 'CC'],
                metric_weight=[1, 1, 1],
                radius_or_number_of_bins=[32, 32, 4],
                sampling_strategy=['Regular', 'Regular', None],
                sampling_percentage=[0.25, 0.25, None],
                number_of_iterations=[
                    [1000, 500, 250, 100],
                    [1000, 500, 250, 100],
                    [100, 100, 70, 20]
                ],
                convergence_threshold=[1e-8, 1e-8, 1e-9],
                convergence_window_size=[10, 10, 15],
                transforms=['Rigid', 'Affine', 'SyN'],
                transform_parameters=[[0.1], [0.1], [0.1, 3, 0]],
                shrink_factors=[
                    [8, 4, 2, 1],
                    [8, 4, 2, 1],
                    [6, 4, 2, 1]
                ],
                smoothing_sigmas=[
                    [3, 2, 1, 0],
                    [3, 2, 1, 0],
                    [3, 2, 1, 0]
                ]
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

                # Not sure to understand how log nodes work yet
                create_log_node(workflow, lesion_preproc,
                                'inputspec.lesion', num_strat)

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

            create_log_node(workflow, ants_reg_anat_mni,
                            'outputspec.normalized_output_brain', num_strat)

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

                node, out_file = strat['anatomical_brain']
                workflow.connect(node, out_file,
                                 flirt_reg_anat_symm_mni,
                                 'inputspec.input_brain')

                # pass the reference files
                workflow.connect(
                    c.template_symmetric_brain_only, 'local_path',
                    flirt_reg_anat_symm_mni, 'inputspec.reference_brain'
                )

                if 'ANTS' in c.regOption:
                    strat = strat.fork()
                    new_strat_list.append(strat)

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

                create_log_node(workflow, flirt_reg_anat_symm_mni,
                                'outputspec.output_brain',
                                num_strat)

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
                    workflow.connect(
                        c.template_brain_only_for_anat, 'local_path',
                        fnirt_reg_anat_symm_mni, 'inputspec.reference_brain'
                    )

                    node, out_file = strat['anatomical_reorient']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.input_skull')

                    node, out_file = strat['anatomical_to_mni_linear_xfm']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.linear_aff')

                    workflow.connect(
                        c.template_symmetric_skull, 'local_path',
                        fnirt_reg_anat_symm_mni, 'inputspec.reference_skull'
                    )

                    workflow.connect(
                        c.dilated_symmetric_brain_mask, 'local_path',
                        fnirt_reg_anat_symm_mni, 'inputspec.ref_mask'
                    )

                    strat.append_name(fnirt_reg_anat_symm_mni.name)
                    strat.set_leaf_properties(fnirt_reg_anat_symm_mni,
                                              'outputspec.output_brain')

                    strat.update_resource_pool({
                        'anatomical_to_symmetric_mni_nonlinear_xfm': (
                        fnirt_reg_anat_symm_mni, 'outputspec.nonlinear_xfm'),
                        'symmetric_anatomical_to_standard': (
                        fnirt_reg_anat_symm_mni, 'outputspec.output_brain')
                    }, override=True)

                    create_log_node(workflow, fnirt_reg_anat_symm_mni,
                                    'outputspec.output_brain',
                                    num_strat)

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
                        num_threads=num_ants_cores
                    )

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
                                     'inputspec.anatomical_brain')

                    # pass the reference file
                    workflow.connect(c.template_symmetric_brain_only, 'local_path',
                                    ants_reg_anat_symm_mni, 'inputspec.reference_brain')

                    # get the reorient skull-on anatomical from resource
                    # pool
                    node, out_file = strat['anatomical_reorient']

                    # pass the anatomical to the workflow
                    workflow.connect(node, out_file,
                                     ants_reg_anat_symm_mni,
                                     'inputspec.anatomical_skull')

                    # pass the reference file
                    workflow.connect(c.template_symmetric_skull, 'local_path',
                                     ants_reg_anat_symm_mni, 'inputspec.reference_skull')


                else:
                    # get the skullstripped anatomical from resource pool
                    node, out_file = strat['anatomical_brain']

                    workflow.connect(node, out_file,
                                     ants_reg_anat_symm_mni,
                                     'inputspec.anatomical_brain')

                    # pass the reference file
                    workflow.connect(c.template_symmetric_brain_only, 'local_path',
                                    ants_reg_anat_symm_mni, 'inputspec.reference_brain')

                ants_reg_anat_symm_mni.inputs.inputspec.set(
                    dimension=3,
                    use_histogram_matching=True,
                    winsorize_lower_quantile=0.01,
                    winsorize_upper_quantile=0.99,
                    metric=['MI', 'MI', 'CC'],
                    metric_weight=[1, 1, 1],
                    radius_or_number_of_bins=[32, 32, 4],
                    sampling_strategy=['Regular', 'Regular', None],
                    sampling_percentage=[0.25, 0.25, None],
                    number_of_iterations=[[1000, 500, 250, 100],
                                          [1000, 500, 250, 100],
                                          [100, 100, 70, 20]],
                    convergence_threshold=[1e-8, 1e-8, 1e-9],
                    convergence_window_size=[10, 10, 15],
                    transforms=['Rigid', 'Affine', 'SyN'],
                    transform_parameters=[[0.1], [0.1], [0.1, 3, 0]],
                    shrink_factors=[[8, 4, 2, 1],
                                    [8, 4, 2, 1],
                                    [6, 4, 2, 1]],
                    smoothing_sigmas=[[3, 2, 1, 0],
                                      [3, 2, 1, 0],
                                      [3, 2, 1, 0]]
                )

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

                    # Not sure to understand how log nodes work yet
                    create_log_node(workflow, lesion_preproc,
                                    'inputspec.lesion', num_strat)

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

                create_log_node(workflow, ants_reg_anat_symm_mni,
                                'outputspec.normalized_output_brain',
                                num_strat)

        strat_list += new_strat_list

    # Inserting Segmentation Preprocessing Workflow

    new_strat_list = []

    if 1 in c.runSegmentationPreprocessing:

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            seg_preproc = None

            # TODO ASH based on config, instead of nodes?
            if 'anat_mni_fnirt_register' in nodes or 'anat_mni_flirt_register' in nodes:
                seg_preproc = create_seg_preproc(use_ants=False,
                                                 wf_name='seg_preproc_%d' % num_strat)
            elif 'anat_mni_ants_register' in nodes:
                seg_preproc = create_seg_preproc(use_ants=True,
                                                 wf_name='seg_preproc_%d' % num_strat)

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

            create_log_node(workflow,
                            seg_preproc, 'outputspec.partial_volume_map',
                            num_strat)

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
                ("scan", func_paths_dict.keys())

            # Add in nodes to get parameters from configuration file
            # a node which checks if scan_parameters are present for each scan
            scan_params = \
                pe.Node(function.Function(input_names=['data_config_scan_params',
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
                                                        'stop_indx'],
                                        function=get_scan_params,
                                        as_module=True),
                        name='scan_params_%d' % num_strat)

            if "Selected Functional Volume" in c.func_reg_input:

                get_func_volume = pe.Node(interface=afni.Calc(),
                                        name='get_func_volume_%d' % num_strat)

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
                pipeconfig_tr=c.TR,
                pipeconfig_tpattern=c.slice_timing_pattern,
                pipeconfig_start_indx=c.startIdx,
                pipeconfig_stop_indx=c.stopIdx
            )

            # node to convert TR between seconds and milliseconds
            convert_tr = pe.Node(function.Function(input_names=['tr'],
                                                output_names=['tr'],
                                                function=get_tr,
                                                as_module=True),
                                name='convert_tr_%d' % num_strat)

            strat.update_resource_pool({
                'raw_functional': (func_wf, 'outputspec.rest'),
                'scan_id': (func_wf, 'outputspec.scan')
            })
  
            strat.set_leaf_properties(func_wf, 'outputspec.rest')

            if 1 in c.runEPI_DistCorr:
                try:
                    if (func_wf, 'outputspec.phase_diff') and (func_wf, 'outputspec.magnitude'):
                        strat.update_resource_pool({
                            "fmap_phase_diff": (func_wf, 'outputspec.phase_diff'),
                            "fmap_magnitude": (func_wf, 'outputspec.magnitude')
                        })
                except:
                    err = "\n\n[!] You have selected to run field map " \
                        "distortion correction, but at least one of your " \
                        "scans listed in your data configuration file is " \
                        "missing either a field map phase difference file " \
                        "or a field map magnitude file, or both.\n\n"
                    logger.warn(err)

            if "Selected Functional Volume" in c.func_reg_input:
                strat.update_resource_pool({
                    'selected_func_volume': (get_func_volume, 'out_file')
                })

        # Truncate scan length based on configuration information

        for num_strat, strat in enumerate(strat_list):
            trunc_wf = create_wf_edit_func(
                wf_name="edit_func_%d" % (num_strat)
            )

            # find the output data on the leaf node
            node, out_file = strat.get_leaf_properties()

            # connect the functional data from the leaf node into the wf
            workflow.connect(node, out_file,
                             trunc_wf, 'inputspec.func')

            # connect the other input parameters
            workflow.connect(scan_params, 'start_indx',
                            trunc_wf, 'inputspec.start_idx')
            workflow.connect(scan_params, 'stop_indx',
                            trunc_wf, 'inputspec.stop_idx')

            # replace the leaf node with the output from the recently added
            # workflow
            strat.set_leaf_properties(trunc_wf, 'outputspec.edited_func')

        # EPI Field-Map based Distortion Correction

        new_strat_list = []

        rp = strat.get_resource_pool()

        if 1 in c.runEPI_DistCorr and 'fmap_phase_diff' in rp.keys() and 'fmap_magnitude' in rp.keys():

            for num_strat, strat in enumerate(strat_list):

                if 'BET' in c.fmap_distcorr_skullstrip:
                    epi_distcorr = create_EPI_DistCorr(
                        use_BET=True,
                        wf_name='epi_distcorr_%d' % (num_strat)
                    )
                    epi_distcorr.inputs.bet_frac_input.bet_frac = c.fmap_distcorr_frac
                    epi_distcorr.get_node('bet_frac_input').iterables = \
                        ('bet_frac', c.fmap_distcorr_frac)
                else:
                    epi_distcorr = create_EPI_DistCorr(
                        use_BET=False,
                        wf_name='epi_distcorr_%d' % (num_strat)
                    )
                    epi_distcorr.inputs.afni_threshold_input.afni_threshold = \
                        c.fmap_distcorr_threshold

                epi_distcorr.inputs.deltaTE_input.deltaTE = c.fmap_distcorr_deltaTE
                epi_distcorr.inputs.dwellT_input.dwellT = c.fmap_distcorr_dwell_time
                epi_distcorr.inputs.dwell_asym_ratio_input.dwell_asym_ratio = c.fmap_distcorr_dwell_asym_ratio

                epi_distcorr.get_node('deltaTE_input').iterables = (
                    'deltaTE', c.fmap_distcorr_deltaTE
                )
                epi_distcorr.get_node('dwellT_input').iterables = (
                    'dwellT', c.fmap_distcorr_dwell_time
                )
                epi_distcorr.get_node('dwell_asym_ratio_input').iterables = (
                    'dwell_asym_ratio', c.fmap_distcorr_dwell_asym_ratio
                )

                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file, epi_distcorr,
                                'inputspec.func_file')

                node, out_file = strat['anatomical_reorient']
                workflow.connect(node, out_file, epi_distcorr,
                                'inputspec.anat_file')

                node, out_file = strat['fmap_phase_diff']
                workflow.connect(node, out_file, epi_distcorr,
                                'inputspec.fmap_pha')

                node, out_file = strat['fmap_magnitude']
                workflow.connect(node, out_file, epi_distcorr,
                                'inputspec.fmap_mag')

                # TODO ASH review forking
                if 0 in c.runEPI_DistCorr:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(epi_distcorr.name)

                strat.update_resource_pool({
                    'despiked_fieldmap': (epi_distcorr, 'outputspec.fmap_despiked'),
                    'fieldmap_mask': (epi_distcorr, 'outputspec.fieldmapmask'),
                    'prepared_fieldmap_map': (epi_distcorr, 'outputspec.fieldmap')
                })

        strat_list += new_strat_list


        # Slice Timing Correction Workflow

        new_strat_list = []

        if 1 in c.slice_timing_correction:

            for num_strat, strat in enumerate(strat_list):

                # create TShift AFNI node
                func_slice_timing_correction = pe.Node(
                    interface=preprocess.TShift(),
                    name='func_slice_timing_correction_%d' % (num_strat))
                func_slice_timing_correction.inputs.outputtype = 'NIFTI_GZ'

                node, out_file = strat.get_leaf_properties()

                workflow.connect(node, out_file,
                                func_slice_timing_correction, 'in_file')

                # TODO ASH normalize TR w schema validation
                # we might prefer to use the TR stored in the NIFTI header
                # if not, use the value in the scan_params node
                if c.TR:
                    if isinstance(c.TR, str):
                        if "None" in c.TR or "none" in c.TR:
                            pass
                        else:
                            workflow.connect(scan_params, 'tr',
                                            func_slice_timing_correction, 'tr')
                    else:
                        workflow.connect(scan_params, 'tr',
                                        func_slice_timing_correction, 'tr')

                if not "Use NIFTI Header" in c.slice_timing_pattern:

                    # add the @ prefix to the tpattern file going into
                    # AFNI 3dTshift - needed this so the tpattern file
                    # output from get_scan_params would be tied downstream
                    # via a connection (to avoid poofing)
                    add_prefix = pe.Node(util.Function(input_names=['tpattern'],
                                                    output_names=[
                                                        'afni_prefix'],
                                                    function=add_afni_prefix),
                                        name='func_slice_timing_correction_add_afni_prefix_%d' % num_strat)
                    workflow.connect(scan_params, 'tpattern',
                                    add_prefix, 'tpattern')
                    workflow.connect(add_prefix, 'afni_prefix',
                                    func_slice_timing_correction,
                                    'tpattern')

                # add the name of the node to the strat name
                strat.append_name(func_slice_timing_correction.name)

                # set the leaf node
                strat.set_leaf_properties(func_slice_timing_correction, 'out_file')

                # add the outputs to the resource pool
                strat.update_resource_pool({
                    'slice_time_corrected': (func_slice_timing_correction, 'out_file')
                })

        # add new strats (if forked)
        strat_list += new_strat_list


        # Functional Image Preprocessing Workflow

        new_strat_list = []

        if '3dAutoMask' in c.functionalMasking:

            for num_strat, strat in enumerate(strat_list):

                func_preproc = create_func_preproc(
                    use_bet=False,
                    wf_name='func_preproc_automask_%d' % num_strat
                )

                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file, func_preproc,
                                'inputspec.func')

                func_preproc.inputs.inputspec.twopass = \
                    getattr(c, 'functional_volreg_twopass', True)

                # TODO ASH review forking
                if 'BET' in c.functionalMasking:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(func_preproc.name)

                strat.set_leaf_properties(func_preproc, 'outputspec.preprocessed')

                # add stuff to resource pool if we need it
                strat.update_resource_pool({
                    'mean_functional': (func_preproc, 'outputspec.example_func'),
                    'functional_preprocessed_mask': (func_preproc, 'outputspec.preprocessed_mask'),
                    'movement_parameters': (func_preproc, 'outputspec.movement_parameters'),
                    'max_displacement': (func_preproc, 'outputspec.max_displacement'),
                    'functional_preprocessed': (func_preproc, 'outputspec.preprocessed'),
                    'functional_brain_mask': (func_preproc, 'outputspec.mask'),
                    'motion_correct': (func_preproc, 'outputspec.motion_correct'),
                    'coordinate_transformation': (func_preproc, 'outputspec.oned_matrix_save')
                })

                create_log_node(workflow, func_preproc,
                                'outputspec.preprocessed', num_strat)

        strat_list += new_strat_list

        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            if 'BET' in c.functionalMasking and 'func_preproc_automask' not in nodes:

                func_preproc = create_func_preproc(use_bet=True,
                                                wf_name='func_preproc_bet_%d' % num_strat)

                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file, func_preproc,
                                'inputspec.func')

                func_preproc.inputs.inputspec.twopass = \
                    getattr(c, 'functional_volreg_twopass', True)

                strat.append_name(func_preproc.name)

                strat.set_leaf_properties(func_preproc, 'outputspec.preprocessed')

                # TODO redundant with above resource pool additions?
                strat.update_resource_pool({
                    'mean_functional': (func_preproc, 'outputspec.example_func'),
                    'functional_preprocessed_mask': (func_preproc, 'outputspec.preprocessed_mask'),
                    'movement_parameters': (func_preproc, 'outputspec.movement_parameters'),
                    'max_displacement': (func_preproc, 'outputspec.max_displacement'),
                    'functional_preprocessed': (func_preproc, 'outputspec.preprocessed'),
                    'functional_brain_mask': (func_preproc, 'outputspec.mask'),
                    'motion_correct': (func_preproc, 'outputspec.motion_correct'),
                    'coordinate_transformation': (func_preproc, 'outputspec.oned_matrix_save'),
                })

                create_log_node(workflow, func_preproc, 'outputspec.preprocessed',
                                num_strat)

        strat_list += new_strat_list


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
                if 'epi_distcorr' in nodes and 1 not in c.runBBReg:
                    dist_corr = True
                    # TODO: for now, disabling dist corr when BBR is disabled
                    err = "\n\n[!] Field map distortion correction is enabled, " \
                        "but Boundary-Based Registration is off- BBR is " \
                        "required for distortion correction.\n\n"
                    raise Exception(err)

                func_to_anat = create_register_func_to_anat(dist_corr,
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

                if dist_corr:
                    # apply field map distortion correction outputs to
                    # the func->anat registration

                    func_to_anat.inputs.echospacing_input.set(
                        echospacing=c.fmap_distcorr_dwell_time[0]
                    )
                    func_to_anat.inputs.pedir_input.set(
                        pedir=c.fmap_distcorr_pedir
                    )

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
                # strat.set_leaf_properties(func_mni_warp, 'out_file')

                strat.update_resource_pool({
                    'mean_functional_in_anat': (func_to_anat, 'outputspec.anat_func_nobbreg'),
                    'functional_to_anat_linear_xfm': (func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg')
                })

                # create_log_node(workflow, func_to_anat, 'outputspec.mni_func', num_strat)

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
                if 'seg_preproc' in nodes:

                    dist_corr = False
                    if 'epi_distcorr' in nodes:
                        dist_corr = True

                    func_to_anat_bbreg = create_bbregister_func_to_anat(
                        dist_corr,
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

                    # Input segmentation probability maps for white matter
                    # segmentation
                    node, out_file = strat['seg_probability_maps']
                    workflow.connect(node, (out_file, pick_wm),
                                    func_to_anat_bbreg,
                                    'inputspec.anat_wm_segmentation')

                    if dist_corr:
                        # apply field map distortion correction outputs to
                        # the func->anat registration

                        func_to_anat_bbreg.inputs.echospacing_input.echospacing = c.fmap_distcorr_dwell_time[
                            0]
                        func_to_anat_bbreg.inputs.pedir_input.pedir = c.fmap_distcorr_pedir

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

                    # create_log_node(workflow, func_to_anat, 'outputspec.mni_func', num_strat)

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

        # Inserting Generate Motion Statistics Workflow

        for num_strat, strat in enumerate(strat_list):

            gen_motion_stats = motion_power_statistics(
                'gen_motion_stats_%d' % num_strat
            )

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
                            gen_motion_stats, 'inputspec.max_displacement')

            node, out_file = strat['functional_brain_mask']
            workflow.connect(node, out_file,
                            gen_motion_stats, 'inputspec.mask')

            node, out_file = strat['coordinate_transformation']
            workflow.connect(node, out_file,
                             gen_motion_stats, 'inputspec.transformations')

            strat.append_name(gen_motion_stats.name)

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

                    # 4D FUNCTIONAL apply warp
                    node, out_file = strat.get_leaf_properties()
                    mean_func_node, mean_func_out_file = strat["mean_functional"]
                    
                    # Insert it on the resource pool, so no need to connect externally
                    ants_apply_warps_func_mni(
                        workflow, strat, num_strat, num_ants_cores,
                        node, out_file,
                        mean_func_node, mean_func_out_file,
                        c.template_brain_only_for_func,
                        "ica_aroma_functional_to_standard",
                        "Linear", 3
                    )

                    aroma_preproc = create_aroma(tr=TR,
                                                 wf_name='create_aroma_%d'
                                                         % num_strat)

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

                    ants_apply_inverse_warps_template_to_func(
                        workflow, strat, num_strat, num_ants_cores, node,
                        out_file, mean_func_node, mean_func_out_file,
                        "ica_aroma_denoised_functional", "Linear", 3
                    )

                    node, out_file = strat["ica_aroma_denoised_functional"]
                    strat.set_leaf_properties(node, out_file)

                    if c.aroma_denoise_type == 'nonaggr':
                        create_log_node(workflow, aroma_preproc,
                                        'outputspec.nonaggr_denoised_file',
                                        num_strat)
                    elif c.aroma_denoise_type == 'aggr':
                        create_log_node(workflow, aroma_preproc,
                                        'outputspec.aggr_denoised_file',
                                        num_strat)

                    strat.append_name(aroma_preproc.name)

        strat_list += new_strat_list


        # Inserting Nuisance Workflow

        new_strat_list = []

        if 1 in c.runNuisance:

            for num_strat, strat in enumerate(strat_list):

                # for each strategy, create a new one without nuisance
                if 0 in c.runNuisance:
                    new_strat_list.append(strat.fork())

                nodes = strat.get_nodes_names()

                has_segmentation = 'seg_preproc' in nodes
                use_ants = 'anat_mni_fnirt_register' not in nodes and 'anat_mni_flirt_register' not in nodes

                for regressors_selector_i, regressors_selector in enumerate(c.Regressors):

                    new_strat = strat.fork()

                    # to guarantee immutability
                    regressors_selector = NuisanceRegressor(
                        copy.deepcopy(regressors_selector),
                        copy.deepcopy(c.Regressors)
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

                    nuisance_regression_workflow = create_nuisance_workflow(
                        regressors_selector,
                        use_ants=use_ants,
                        name='nuisance_{0}_{1}'.format(regressors_selector_i, num_strat)
                    )

                    node, out_file = new_strat['anatomical_brain']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_workflow, 'inputspec.anatomical_file_path'
                    )

                    if has_segmentation:

                        workflow.connect(
                            c.lateral_ventricles_mask, 'local_path',
                            nuisance_regression_workflow, 'inputspec.lat_ventricles_mask_file_path'
                        )

                        node, out_file = new_strat['anatomical_gm_mask']
                        workflow.connect(
                            node, out_file,
                            nuisance_regression_workflow, 'inputspec.gm_mask_file_path'
                        )

                        node, out_file = new_strat['anatomical_wm_mask']
                        workflow.connect(
                            node, out_file,
                            nuisance_regression_workflow, 'inputspec.wm_mask_file_path'
                        )

                        node, out_file = new_strat['anatomical_csf_mask']
                        workflow.connect(
                            node, out_file,
                            nuisance_regression_workflow, 'inputspec.csf_mask_file_path'
                        )

                    node, out_file = new_strat['movement_parameters']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_workflow,
                        'inputspec.motion_parameters_file_path'
                    )

                    node, out_file= new_strat['functional_to_anat_linear_xfm']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_workflow,
                        'inputspec.func_to_anat_linear_xfm_file_path'
                    )

                    node, out_file = new_strat.get_leaf_properties()
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_workflow,
                        'inputspec.functional_file_path'
                    )

                    node, out_file = new_strat['frame_wise_displacement_jenkinson']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_workflow,
                        'inputspec.fd_j_file_path'
                    )

                    node, out_file = new_strat['frame_wise_displacement_power']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_workflow,
                        'inputspec.fd_p_file_path'
                    )

                    node, out_file = new_strat['dvars']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_workflow,
                        'inputspec.dvars_file_path'
                    )

                    node, out_file = new_strat['functional_brain_mask']
                    workflow.connect(
                        node, out_file,
                        nuisance_regression_workflow,
                        'inputspec.functional_brain_mask_file_path'
                    )

                    nuisance_regression_workflow.get_node('inputspec').iterables = ([
                        ('selector', [regressors_selector]),
                    ])

                    if use_ants:

                        # pass the ants_affine_xfm to the input for the
                        # INVERSE transform, but ants_affine_xfm gets inverted
                        # within the workflow

                        node, out_file = new_strat['ants_initial_xfm']
                        workflow.connect(
                            node, out_file,
                            nuisance_regression_workflow,
                            'inputspec.anat_to_mni_initial_xfm_file_path'
                        )

                        node, out_file = new_strat['ants_rigid_xfm']
                        workflow.connect(
                            node, out_file,
                            nuisance_regression_workflow,
                            'inputspec.anat_to_mni_rigid_xfm_file_path'
                        )

                        node, out_file = new_strat['ants_affine_xfm']
                        workflow.connect(
                            node, out_file,
                            nuisance_regression_workflow,
                            'inputspec.anat_to_mni_affine_xfm_file_path'
                        )
                    else:
                        node, out_file = new_strat['mni_to_anatomical_linear_xfm']
                        workflow.connect(
                            node, out_file,
                            nuisance_regression_workflow,
                            'inputspec.mni_to_anat_linear_xfm_file_path'
                        )


                    new_strat.append_name(nuisance_regression_workflow.name)

                    new_strat.set_leaf_properties(
                        nuisance_regression_workflow,
                        'outputspec.residual_file_path'
                    )

                    new_strat.update_resource_pool({
                        'nuisance_regression_selector': regressors_selector,
                        
                        'functional_nuisance_residuals': (
                            nuisance_regression_workflow,
                            'outputspec.residual_file_path')
                        ,
                        'functional_nuisance_regressors': (
                            nuisance_regression_workflow,
                            'outputspec.regressors_file_path'
                        ),
                    })

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

                create_log_node(workflow,
                                median_angle_corr, 'outputspec.subject',
                                num_strat)

        strat_list += new_strat_list

        for num_strat, strat in enumerate(strat_list):
            # Keep non-bandpassed version of the output for ALFF
            strat.update_resource_pool({
                'functional_freq_unfiltered': strat.get_leaf_properties()
            })

        # Inserting Bandpassing Workflow
        for num_strat, strat in enumerate(strat_list):

            if 'nuisance_regression_selector' not in strat:
                continue

            if not strat['nuisance_regression_selector'].get('Bandpass'):
                continue

            bandpass_selector = strat['nuisance_regression_selector']['Bandpass']

            frequency_filter = pe.Node(
                function.Function(input_names=['realigned_file',
                                               'bandpass_freqs',
                                               'sample_period'],
                                  output_names=['bandpassed_file'],
                                  function=bandpass_voxels,
                                  as_module=True),
                name='frequency_filter_%d' % num_strat
            )

            frequency_filter.inputs.bandpass_freqs = [
                bandpass_selector.get('bottom_frequency'),
                bandpass_selector.get('top_frequency')
            ]

            node, out_file = strat.get_leaf_properties()
            workflow.connect(node, out_file,
                             frequency_filter, 'realigned_file')

            strat.append_name(frequency_filter.name)

            strat.set_leaf_properties(frequency_filter, 'bandpassed_file')
            strat.update_resource_pool({
                'functional_freq_filtered': (frequency_filter, 'bandpassed_file')
            })


        # Func -> Template, uses antsApplyTransforms (ANTS) or ApplyWarp (FSL) to
        #  apply the warp; also includes mean functional warp
        new_strat_list = []

        if 1 in c.runRegisterFuncToMNI:

            for num_strat, strat in enumerate(strat_list):

                nodes = strat.get_nodes_names()

                # Run FSL ApplyWarp
                if 'anat_mni_flirt_register' in nodes or 'anat_mni_fnirt_register' in nodes:

                    func_mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                                            name='func_mni_fsl_warp_%d' % num_strat)
                    func_mni_warp.inputs.ref_file = c.template_brain_only_for_func

                    functional_brain_mask_to_standard = pe.Node(
                        interface=fsl.ApplyWarp(),
                        name='func_mni_fsl_warp_mask_%d' % num_strat
                    )
                    functional_brain_mask_to_standard.inputs.interp = 'nn'
                    functional_brain_mask_to_standard.inputs.ref_file = c.template_skull_for_func

                    mean_functional_warp = pe.Node(
                        interface=fsl.ApplyWarp(),
                        name='mean_func_fsl_warp_%d' % num_strat
                    )
                    mean_functional_warp.inputs.ref_file = c.template_brain_only_for_func

                    motion_correct_warp = pe.Node(
                        interface=fsl.ApplyWarp(),
                        name="motion_correct_fsl_warp_%d" % num_strat
                    )
                    motion_correct_warp.inputs.ref_file = c.template_brain_only_for_func

                    if 'anat_mni_fnirt_register' in nodes:
                        node, out_file = strat['anatomical_to_mni_nonlinear_xfm']
                        workflow.connect(node, out_file,
                                         func_mni_warp, 'field_file')
                        workflow.connect(node, out_file,
                                         functional_brain_mask_to_standard, 'field_file')
                        workflow.connect(node, out_file,
                                         mean_functional_warp, 'field_file')
                        workflow.connect(node, out_file,
                                         motion_correct_warp, 'field_file')

                        node, out_file = strat['functional_to_anat_linear_xfm']
                        workflow.connect(node, out_file,
                                         func_mni_warp, 'premat')
                        workflow.connect(node, out_file,
                                         functional_brain_mask_to_standard, 'premat')
                        workflow.connect(node, out_file,
                                         mean_functional_warp, 'premat')
                        workflow.connect(node, out_file,
                                         motion_correct_warp, 'premat')

                        node, out_file = strat.get_leaf_properties()
                        workflow.connect(node, out_file,
                                         func_mni_warp, 'in_file')

                        node, out_file = strat['functional_brain_mask']
                        workflow.connect(node, out_file,
                                         functional_brain_mask_to_standard, 'in_file')

                        node, out_file = strat['mean_functional']
                        workflow.connect(node, out_file,
                                         mean_functional_warp, 'in_file')

                        node, out_file = strat['motion_correct']
                        workflow.connect(node, out_file,
                                         motion_correct_warp, 'in_file')

                    elif 'anat_mni_flirt_register' in nodes:
                        func_anat_warp = pe.Node(interface=fsl.ApplyWarp(),
	                                             name='func_anat_fsl_warp_%d' % num_strat)
                        functional_brain_mask_to_anat = pe.Node(
	                        interface=fsl.ApplyWarp(),
	                        name='func_anat_fsl_warp_mask_%d' % num_strat
	                    )
                        functional_brain_mask_to_anat.inputs.interp = 'nn'

                        mean_functional_to_anat = pe.Node(
                            interface=fsl.ApplyWarp(),
	                        name='mean_func_to_anat_fsl_warp_%d' % num_strat
	                    )

                        motion_correct_to_anat_warp = pe.Node(
	                        interface=fsl.ApplyWarp(),
	                        name="motion_correct_to_anat_fsl_warp_%d" % num_strat
	                    )

                        node, out_file = strat.get_leaf_properties()
                        workflow.connect(node, out_file,
                                         func_anat_warp, 'in_file')

                        node, out_file = strat['functional_brain_mask']
                        workflow.connect(node, out_file,
                                         functional_brain_mask_to_anat, 'in_file')

                        node, out_file = strat['mean_functional']
                        workflow.connect(node, out_file,
                                         mean_functional_to_anat, 'in_file')

                        node, out_file = strat['motion_correct']
                        workflow.connect(node, out_file,
                                         motion_correct_to_anat_warp, 'in_file')

                        node, out_file = strat['anatomical_brain']
                        workflow.connect(node, out_file,
                                         func_anat_warp, 'ref_file')
                        workflow.connect(node, out_file,
                                         functional_brain_mask_to_anat, 'ref_file')
                        workflow.connect(node, out_file,
                                         mean_functional_to_anat, 'ref_file')
                        workflow.connect(node, out_file,
                                         motion_correct_to_anat_warp, 'ref_file') 

                        node, out_file = strat['functional_to_anat_linear_xfm']
                        workflow.connect(node, out_file,
                                         func_anat_warp, 'premat')
                        workflow.connect(node, out_file,
                                         functional_brain_mask_to_anat, 'premat')
                        workflow.connect(node, out_file,
                                         mean_functional_to_anat, 'premat')
                        workflow.connect(node, out_file,
                                         motion_correct_to_anat_warp, 'premat')

                        node, out_file = strat.get_leaf_properties()
                        workflow.connect(func_anat_warp, 'out_file',
                                         func_mni_warp, 'in_file')

                        workflow.connect(functional_brain_mask_to_anat, 'out_file',
                                         functional_brain_mask_to_standard, 'in_file')

                        workflow.connect(mean_functional_to_anat, 'out_file',
                                         mean_functional_warp, 'in_file')

                        workflow.connect(motion_correct_to_anat_warp, 'out_file',
                                         motion_correct_warp, 'in_file')

                        node, out_file = strat['anatomical_to_mni_linear_xfm']
                        workflow.connect(node, out_file,
                                         func_mni_warp, 'premat')
                        workflow.connect(node, out_file,
                                         functional_brain_mask_to_standard, 'premat')
                        workflow.connect(node, out_file,
                                         mean_functional_warp, 'premat')
                        workflow.connect(node, out_file,
                                         motion_correct_warp, 'premat')

                    strat.update_resource_pool({
                        'functional_to_standard': (func_mni_warp, 'out_file'),
                        'functional_brain_mask_to_standard': (functional_brain_mask_to_standard, 'out_file'),
                        'mean_functional_to_standard': (mean_functional_warp, 'out_file'),
                        'motion_correct_to_standard': (motion_correct_warp, 'out_file')
                    })

                    strat.append_name(func_mni_warp.name)
                    create_log_node(workflow,
                                    func_mni_warp, 'out_file',
                                    num_strat)

            strat_list += new_strat_list

            for num_strat, strat in enumerate(strat_list):

                nodes = strat.get_nodes_names()

                if 'ANTS' in c.regOption and \
                    'anat_mni_flirt_register' not in nodes and \
                        'anat_mni_fnirt_register' not in nodes:

                    # ANTS warp application

                    # 4D FUNCTIONAL apply warp
                    node, out_file = strat.get_leaf_properties()
                    node2, out_file2 = \
                        strat["mean_functional"]

                    warp_func_wf = ants_apply_warps_func_mni(
                        workflow, strat, num_strat, num_ants_cores,
                        node, out_file,
                        node2, out_file2,
                        c.template_brain_only_for_func,
                        "functional_to_standard",
                        "Linear", 3
                    )

                    create_log_node(workflow, warp_func_wf,
                                    'outputspec.output_image', num_strat)

                    # 4D FUNCTIONAL MOTION-CORRECTED apply warp
                    node, out_file = \
                        strat['motion_correct']
                    node2, out_file2 = \
                        strat["mean_functional"]

                    warp_motion_wf = ants_apply_warps_func_mni(
                        workflow, strat, num_strat, num_ants_cores,
                        node, out_file,
                        node2, out_file2,
                        c.template_brain_only_for_func,
                        "motion_correct_to_standard",
                        "Linear", 3
                    )

                    create_log_node(workflow, warp_motion_wf,
                                    'outputspec.output_image', num_strat)

                    # FUNCTIONAL BRAIN MASK (binary, no timeseries) apply warp
                    node, out_file = \
                        strat["functional_brain_mask"]

                    warp_mask_wf = ants_apply_warps_func_mni(
                        workflow, strat, num_strat, num_ants_cores,
                        node, out_file,
                        node, out_file,
                        c.template_brain_only_for_func,
                        "functional_brain_mask_to_standard",
                        "NearestNeighbor", 0
                    )

                    create_log_node(workflow, warp_mask_wf,
                                    'outputspec.output_image', num_strat)

                    # FUNCTIONAL MEAN (no timeseries) apply warp
                    node, out_file = \
                        strat["mean_functional"]

                    warp_mean_wf = ants_apply_warps_func_mni(
                        workflow, strat, num_strat, num_ants_cores,
                        node, out_file,
                        node, out_file,
                        c.template_brain_only_for_func,
                        "mean_functional_to_standard",
                        "Linear", 0
                    )

                    create_log_node(workflow, warp_mean_wf,
                                    'outputspec.output_image', num_strat)

        strat_list += new_strat_list
        
        # Derivatives

        # Inserting ALFF/fALFF workflow
        #     NOTE: this is calculated using the functional time series from
        #           before frequency filtering and beyond
        new_strat_list = []

        if 1 in c.runALFF:
            for num_strat, strat in enumerate(strat_list):

                alff = create_alff('alff_falff_%d' % num_strat)

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

                create_log_node(workflow,
                                alff, 'outputspec.falff_img', num_strat)

        strat_list += new_strat_list

        # Inserting VMHC Workflow

        new_strat_list = []

        if 1 in c.runVMHC:

            for num_strat, strat in enumerate(strat_list):

                nodes = strat.get_nodes_names()

                if 'func_mni_fsl_warp' in nodes:
                    if 'anat_mni_fnirt_register' not in nodes and 'anat_mni_flirt_register' in nodes:
                        vmhc = create_vmhc(False, True, 'vmhc_%d' % num_strat)
                    elif 'anat_mni_fnirt_register' in nodes:
                        vmhc = create_vmhc(False, False, 'vmhc_%d' % num_strat)
                else:
                    vmhc = create_vmhc(True, False, 'vmhc_%d' % num_strat,
                                       int(num_ants_cores))

                vmhc.inputs.inputspec.standard_for_func = c.template_skull_for_func
                vmhc.inputs.fwhm_input.fwhm = c.fwhm
                vmhc.get_node('fwhm_input').iterables = ('fwhm', c.fwhm)

                node, out_file = strat.get_leaf_properties()
                workflow.connect(node, out_file,
                                vmhc, 'inputspec.rest_res')

                node, out_file = strat['functional_to_anat_linear_xfm']
                workflow.connect(node, out_file,
                                vmhc, 'inputspec.example_func2highres_mat')

                node, out_file = strat['functional_brain_mask']
                workflow.connect(node, out_file,
                                vmhc, 'inputspec.rest_mask')

                node, out_file = strat['mean_functional']
                workflow.connect(node, out_file,
                                vmhc, 'inputspec.mean_functional')

                node, out_file = strat['anatomical_brain']
                workflow.connect(node, out_file,
                                vmhc, 'inputspec.brain')

                # TODO ASH normalize w schema val
                if 'ANTS' in c.regOption and \
                    'anat_mni_flirt_register' not in nodes and \
                    'anat_mni_fnirt_register' not in nodes and \
                    'anat_symmetric_mni_flirt_register' not in nodes and \
                    'anat_symmetric_mni_fnirt_register' not in nodes:

                    node, out_file = strat['ants_symmetric_initial_xfm']
                    workflow.connect(node, out_file,
                                    vmhc, 'inputspec.ants_symm_initial_xfm')

                    node, out_file = strat['ants_symmetric_rigid_xfm']
                    workflow.connect(node, out_file,
                                    vmhc, 'inputspec.ants_symm_rigid_xfm')

                    node, out_file = strat['ants_symmetric_affine_xfm']
                    workflow.connect(node, out_file,
                                    vmhc, 'inputspec.ants_symm_affine_xfm')

                    node, out_file = strat['anatomical_to_symmetric_mni_nonlinear_xfm']
                    workflow.connect(node, out_file,
                                    vmhc, 'inputspec.ants_symm_warp_field')

                else:
                    if 'anat_mni_fnirt_register' in nodes:
                        node, out_file = strat['anatomical_to_symmetric_mni_nonlinear_xfm']
                        workflow.connect(node, out_file,
                                         vmhc, 'inputspec.fnirt_nonlinear_warp')
                    elif 'anat_mni_flirt_register' in nodes:
                        node, out_file = strat[
                            'anatomical_to_symmetric_mni_linear_xfm']
                        workflow.connect(node, out_file,
                                         vmhc,
                                         'inputspec.flirt_linear_aff')

                strat.update_resource_pool({
                    'vmhc_raw_score': (vmhc, 'outputspec.VMHC_FWHM_img'),
                    'vmhc_fisher_zstd': (vmhc, 'outputspec.VMHC_Z_FWHM_img'),
                    'vmhc_fisher_zstd_zstat_map': (vmhc, 'outputspec.VMHC_Z_stat_FWHM_img')
                })

                strat.append_name(vmhc.name)

                create_log_node(
                    workflow, vmhc, 'outputspec.VMHC_FWHM_img', num_strat)

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
                
                create_log_node(workflow, reho, 'outputspec.raw_reho_map',
                                num_strat)

        ts_analysis_dict = {}
        sca_analysis_dict = {}

        # TODO ASH normalize w schema val
        if c.tsa_roi_paths:

            tsa_roi_dict = c.tsa_roi_paths[0]

            # Timeseries and SCA config selections processing

            # flip the dictionary
            for roi_path in tsa_roi_dict.keys():

                ts_analysis_to_run = map(
                    lambda x: x.strip(),
                    tsa_roi_dict[roi_path].split(",")
                )

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

            # c.tsa_roi_paths and c.sca_roi_paths come in a format as such:
            # a list containing a dictionary
            # [
            #     {
            #         '/path/to/rois1.nii.gz': 'Avg, MultReg',
            #         '/path/to/rois2.nii.gz': 'Avg, MultReg',
            #         '/path/to/rois3.nii.gz': 'Avg, MultReg',
            #         '/path/to/rois4.nii.gz': 'DualReg'
            #     }
            # ]

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

                for analysis_type in sca_roi_dict[roi_path].split(","):
                    analysis_type = analysis_type.replace(" ", "")

                    if analysis_type not in sca_analysis_dict.keys():
                        sca_analysis_dict[analysis_type] = []

                    sca_analysis_dict[analysis_type].append(roi_path)

        # Section: Spatial Regression Based Time Series

        new_strat_list = []

        if "SpatialReg" in ts_analysis_dict.keys() or \
            "DualReg" in sca_analysis_dict.keys():

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

                    create_log_node(workflow, spatial_map_timeseries,
                                    'outputspec.subject_timeseries', num_strat)

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

                    create_log_node(workflow, spatial_map_timeseries_for_dr,
                                    'outputspec.subject_timeseries',
                                    num_strat)

        strat_list += new_strat_list

        if 1 in c.runROITimeseries and ("Avg" in ts_analysis_dict.keys() or \
            "Avg" in sca_analysis_dict.keys() or \
            "MultReg" in sca_analysis_dict.keys()):

            # ROI Based Time Series
            new_strat_list = []

            for num_strat, strat in enumerate(strat_list):

                if "Avg" in ts_analysis_dict.keys():
                    resample_functional_to_roi = pe.Node(interface=fsl.FLIRT(),
                                                        name='resample_functional_to_roi_%d' % num_strat)

                    resample_functional_to_roi.inputs.set(
                        interp='trilinear',
                        apply_xfm=True,
                        in_matrix_file=c.identityMatrix
                    )

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
                                     resample_functional_to_roi, 'in_file')
                    workflow.connect(roi_dataflow, 'outputspec.out_file',
                                     resample_functional_to_roi, 'reference')

                    # connect it to the roi_timeseries
                    workflow.connect(roi_dataflow, 'outputspec.out_file',
                                     roi_timeseries, 'input_roi.roi')
                    workflow.connect(resample_functional_to_roi, 'out_file',
                                     roi_timeseries, 'inputspec.rest')

                    strat.append_name(roi_timeseries.name)
                    strat.update_resource_pool({
                        'roi_timeseries': (roi_timeseries, 'outputspec.roi_outputs'),
                        'functional_to_roi': (resample_functional_to_roi, 'out_file')
                    })
                    create_log_node(workflow, roi_timeseries, 'outputspec.roi_outputs',
                                    num_strat)

                if "Avg" in sca_analysis_dict.keys():

                    # same workflow, except to run TSE and send it to the resource
                    # pool so that it will not get sent to SCA
                    resample_functional_to_roi_for_sca = pe.Node(
                        interface=fsl.FLIRT(),
                        name='resample_functional_to_roi_for_sca_%d' % num_strat
                    )

                    resample_functional_to_roi_for_sca.inputs.set(
                        interp='trilinear',
                        apply_xfm=True,
                        in_matrix_file=c.identityMatrix
                    )

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
                                     resample_functional_to_roi_for_sca,
                                     'in_file')
                    workflow.connect(roi_dataflow_for_sca,
                                     'outputspec.out_file',
                                     resample_functional_to_roi_for_sca,
                                     'reference')

                    # connect it to the roi_timeseries
                    workflow.connect(roi_dataflow_for_sca,
                                     'outputspec.out_file',
                                     roi_timeseries_for_sca, 'input_roi.roi')
                    workflow.connect(resample_functional_to_roi_for_sca,
                                     'out_file',
                                     roi_timeseries_for_sca, 'inputspec.rest')

                    strat.append_name(roi_timeseries_for_sca.name)
                    strat.update_resource_pool({
                        'roi_timeseries_for_SCA': (roi_timeseries_for_sca, 'outputspec.roi_outputs'),
                        'functional_to_roi_for_SCA': (resample_functional_to_roi, 'out_file')

                    })
                    create_log_node(workflow, roi_timeseries_for_sca,
                                    'outputspec.roi_outputs', num_strat)

                if "MultReg" in sca_analysis_dict.keys():

                    # same workflow, except to run TSE and send it to the resource
                    # pool so that it will not get sent to SCA
                    resample_functional_to_roi_for_multreg = pe.Node(
                        interface=fsl.FLIRT(),
                        name='resample_functional_to_roi_for_mult_reg_%d' % num_strat
                    )

                    resample_functional_to_roi_for_multreg.inputs.set(
                        interp='trilinear',
                        apply_xfm=True,
                        in_matrix_file=c.identityMatrix
                    )

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
                                    resample_functional_to_roi_for_multreg,
                                    'in_file')
                    workflow.connect(roi_dataflow_for_multreg,
                                    'outputspec.out_file',
                                    resample_functional_to_roi_for_multreg,
                                    'reference')

                    # connect it to the roi_timeseries
                    workflow.connect(roi_dataflow_for_multreg,
                                    'outputspec.out_file',
                                    roi_timeseries_for_multreg,
                                    'input_roi.roi')
                    workflow.connect(resample_functional_to_roi_for_multreg,
                                    'out_file',
                                    roi_timeseries_for_multreg,
                                    'inputspec.rest')

                    strat.append_name(roi_timeseries_for_multreg.name)
                    strat.update_resource_pool({
                        'roi_timeseries_for_SCA_multreg': (roi_timeseries_for_multreg, 'outputspec.roi_outputs')
                    })
                    create_log_node(workflow, roi_timeseries_for_multreg,
                                    'outputspec.roi_outputs', num_strat)

        strat_list += new_strat_list


        # Connectome
        if "PearsonCorr" in ts_analysis_dict.keys() or \
            "PartialCorr" in ts_analysis_dict.keys():

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

                resample_functional_to_mask = pe.Node(interface=fsl.FLIRT(),
                                                    name='resample_functional_to_mask_%d' % num_strat)
                resample_functional_to_mask.inputs.set(
                    interp='trilinear',
                    apply_xfm=True,
                    in_matrix_file=c.identityMatrix
                )

                mask_dataflow = create_roi_mask_dataflow(ts_analysis_dict["Voxel"],
                                                        'mask_dataflow_%d' % num_strat)

                voxel_timeseries = get_voxel_timeseries(
                    'voxel_timeseries_%d' % num_strat)
                voxel_timeseries.inputs.inputspec.output_type = c.roiTSOutputs

                node, out_file = strat['functional_to_standard']

                # resample the input functional file to mask
                workflow.connect(node, out_file,
                                resample_functional_to_mask, 'in_file')
                workflow.connect(mask_dataflow, 'outputspec.out_file',
                                resample_functional_to_mask, 'reference')

                # connect it to the voxel_timeseries
                workflow.connect(mask_dataflow, 'outputspec.out_file',
                                voxel_timeseries, 'input_mask.mask')
                workflow.connect(resample_functional_to_mask, 'out_file',
                                voxel_timeseries, 'inputspec.rest')

                strat.append_name(voxel_timeseries.name)
                strat.update_resource_pool({
                    'voxel_timeseries': (voxel_timeseries, 'outputspec.mask_outputs')
                })

                create_log_node(workflow, voxel_timeseries,
                                'outputspec.mask_outputs', num_strat)

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

                create_log_node(workflow,
                                sca_roi, 'outputspec.correlation_stack',
                                num_strat)

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

                create_log_node(workflow, dr_temp_reg,
                                'outputspec.temp_reg_map', num_strat)

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

                create_log_node(workflow, sc_temp_reg,
                                'outputspec.temp_reg_map', num_strat)

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

            if 1 in c.runRegisterFuncToMNI:
                rp = strat.get_resource_pool()
                for key in sorted(rp.keys()):
                    # connect nodes to apply warps to template
                    if key in Outputs.native_nonsmooth:
                        # smoothing happens at the end, so only the non-smooth
                        # named output labels for the native-space outputs
                        strat = output_to_standard(
                            workflow, key, strat, num_strat, c)
                    elif key in Outputs.native_nonsmooth_mult:
                        strat = output_to_standard(workflow, key, strat, num_strat, c,
                                                map_node=True)

            if "Before" in c.smoothing_order:
                # run smoothing before Z-scoring
                if 1 in c.run_smoothing:
                    rp = strat.get_resource_pool()
                    for key in sorted(rp.keys()):
                        # connect nodes for smoothing
                        if "centrality" in key:
                            # centrality needs its own mask
                            strat = output_smooth(workflow, key,
                                                c.templateSpecificationFile, c.fwhm,
                                                strat, num_strat, map_node=True)
                        elif key in Outputs.native_nonsmooth:
                            # native space
                            strat = output_smooth(workflow, key, "functional_brain_mask", c.fwhm,
                                                strat, num_strat)
                        elif key in Outputs.native_nonsmooth_mult:
                            # native space with multiple files (map nodes)
                            strat = output_smooth(workflow, key, "functional_brain_mask", c.fwhm,
                                                strat, num_strat, map_node=True)
                        elif key in Outputs.template_nonsmooth:
                            # template space
                            strat = output_smooth(workflow, key,
                                                "functional_brain_mask_to_standard", c.fwhm,
                                                strat, num_strat)
                        elif key in Outputs.template_nonsmooth_mult:
                            # template space with multiple files (map nodes)
                            strat = output_smooth(workflow, key,
                                                "functional_brain_mask_to_standard", c.fwhm,
                                                strat, num_strat, map_node=True)

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
                                                        "functional_brain_mask_to_standard",
                                                        strat, num_strat)
                        elif key in Outputs.template_raw_mult:
                            # same as above but multiple files so mapnode required
                            strat = z_score_standardize(workflow, key,
                                                        "functional_brain_mask_to_standard",
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
                                                        "functional_brain_mask_to_standard",
                                                        strat, num_strat)
                        elif key in Outputs.template_raw_mult:
                            # same as above but multiple files so mapnode required
                            strat = z_score_standardize(workflow, key,
                                                        "functional_brain_mask_to_standard",
                                                        strat, num_strat,
                                                        map_node=True)

                if 1 in c.run_smoothing:
                    rp = strat.get_resource_pool()
                    for key in sorted(rp.keys()):
                        # connect nodes for smoothing
                        if "centrality" in key:
                            # centrality needs its own mask
                            strat = output_smooth(workflow, key,
                                                c.templateSpecificationFile, c.fwhm,
                                                strat, num_strat, map_node=True)
                        elif key in Outputs.native_nonsmooth:
                            # native space
                            strat = output_smooth(workflow, key, "functional_brain_mask", c.fwhm,
                                                strat, num_strat)
                        elif key in Outputs.native_nonsmooth_mult:
                            # native space with multiple files (map nodes)
                            strat = output_smooth(workflow, key, "functional_brain_mask", c.fwhm,
                                                strat, num_strat, map_node=True)
                        elif key in Outputs.template_nonsmooth:
                            # template space
                            strat = output_smooth(workflow, key,
                                                "functional_brain_mask_to_standard", c.fwhm,
                                                strat, num_strat)
                        elif key in Outputs.template_nonsmooth_mult:
                            # template space with multiple files (map nodes)
                            strat = output_smooth(workflow, key,
                                                "functional_brain_mask_to_standard", c.fwhm,
                                                strat, num_strat, map_node=True)

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
    qc_montage_id_a = {}
    qc_montage_id_s = {}
    qc_plot_id = {}
    qc_hist_id = {}

    raising = None

    if 1 in c.generateQualityControlImages:
        qc_montage_id_a, qc_montage_id_s, qc_hist_id, qc_plot_id = \
            create_qc_workflow(workflow, c, strat_list, Outputs.qc)

    logger.info('\n\n' + 'Pipeline building completed.' + '\n\n')

    # Run the pipeline only if the user signifies.
    # otherwise, only construct the pipeline (above)
    if run == 1:

        try:
            workflow.write_graph(graph2use='hierarchical')
        except:
            pass

        # this section creates names for the different branched strategies.
        # it identifies where the pipeline has forked and then appends the
        # name of the forked nodes to the branch name in the output directory

        fork_points_labels = Strategy.get_forking_labels(strat_list)

        # DataSink
        pipeline_ids = []

        scan_ids = ['scan_anat']
        if 'func' in sub_dict:
            scan_ids += ['scan_' + str(scan_id)
                         for scan_id in sub_dict['func']]
        if 'rest' in sub_dict:
            scan_ids += ['scan_' + str(scan_id)
                         for scan_id in sub_dict['rest']]


        for num_strat, strat in enumerate(strat_list):

            if p_name is None or p_name == 'None':
                pipeline_id = c.pipelineName
            else:
                pipeline_id = p_name

            if fork_points_labels[strat]:
                pipeline_id += '_' + fork_points_labels[strat]

            pipeline_ids.append(pipeline_id)


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


            # TODO enforce value with schema validation
            try:
                encrypt_data = bool(c.s3Encryption[0])
            except:
                encrypt_data = False

            ndmg_out = False
            try:
                # let's encapsulate this inside a Try..Except block so if
                # someone doesn't have ndmg_outputs in their pipe config,
                # it will default to the regular datasink
                #     TODO: update this when we change to the optionals
                #     TODO: only pipe config
                if 1 in c.ndmg_mode:
                    ndmg_out = True
            except:
                pass

            if ndmg_out:
                # create the graphs
                from CPAC.utils.ndmg_utils import (
                    ndmg_roi_timeseries,
                    ndmg_create_graphs
                )

                atlases = []
                if 'Avg' in ts_analysis_dict.keys():
                    atlases = ts_analysis_dict['Avg']

                roi_dataflow_for_ndmg = create_roi_mask_dataflow(atlases,
                    'roi_dataflow_for_ndmg_%d' % num_strat
                )

                resample_functional_to_roi = pe.Node(interface=fsl.FLIRT(),
                                                     name='resample_functional_to_roi_ndmg_%d' % num_strat)
                resample_functional_to_roi.inputs.set(
                    interp='trilinear',
                    apply_xfm=True,
                    in_matrix_file=c.identityMatrix
                )
                workflow.connect(roi_dataflow_for_ndmg, 'outputspec.out_file',
                                    resample_functional_to_roi, 'reference')

                ndmg_ts = pe.Node(function.Function(
                    input_names=['func_file',
                                 'label_file'],
                    output_names=['roi_ts',
                                  'rois',
                                  'roits_file'],
                    function=ndmg_roi_timeseries,
                    as_module=True
                ), name='ndmg_ts_%d' % num_strat)

                node, out_file = strat['functional_to_standard']
                workflow.connect(node, out_file,
                                 resample_functional_to_roi, 'in_file')
                workflow.connect(resample_functional_to_roi, 'out_file',
                                 ndmg_ts, 'func_file')
                workflow.connect(roi_dataflow_for_ndmg, 'outputspec.out_file',
                                 ndmg_ts, 'label_file')

                ndmg_graph = pe.MapNode(function.Function(
                    input_names=['ts', 'labels'],
                    output_names=['out_file'],
                    function=ndmg_create_graphs,
                    as_module=True
                ), name='ndmg_graphs_%d' % num_strat, iterfield=['labels'])

                workflow.connect(ndmg_ts, 'roi_ts', ndmg_graph, 'ts')
                workflow.connect(roi_dataflow_for_ndmg, 'outputspec.out_file',
                                 ndmg_graph, 'labels')

                strat.update_resource_pool({
                    'ndmg_ts': (ndmg_ts, 'roits_file'),
                    'ndmg_graph': (ndmg_graph, 'out_file')
                })


            rp = strat.get_resource_pool()

            if c.write_debugging_outputs:
                workdir = os.path.join(c.workingDirectory, workflow_name)
                rp_pkl = os.path.join(workdir, 'resource_pool.pkl')
                with open(rp_pkl, 'wt') as f:
                    pickle.dump(rp, f)

            output_sink_nodes = []

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

                    if 0 not in c.runRegisterFuncToMNI:
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

                    anat_res_tag = c.resolution_for_anat.replace('mm', '')
                    func_res_tag = c.resolution_for_func_preproc.replace('mm', '')

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
                        'functional_mask_to_standard': (
                            'func',
                            'registered',
                            '{0}_bold_space-{1}_res-{2}x{2}x{2}_registered_mask'
                            .format(id_tag, func_template_tag, func_res_tag)
                        ),
                        'ndmg_ts': (
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
                    # regular datasink
                    ds = pe.Node(
                        DataSink(),
                        name='sinker_{}_{}'.format(num_strat, resource_i)
                    )
                    ds.inputs.base_directory = c.outputDirectory
                    ds.inputs.creds_path = creds_path
                    ds.inputs.encrypt_bucket_keys = encrypt_data
                    ds.inputs.container = os.path.join(
                        'pipeline_%s' % pipeline_id, subject_id
                    )
                    ds.inputs.regexp_substitutions = [
                        (r"/_sca_roi(.)*[/]", '/'),
                        (r"/_smooth_centrality_(\d)+[/]", '/'),
                        (r"/_z_score(\d)+[/]", "/"),
                        (r"/_dr_tempreg_maps_zstat_files_smooth_(\d)+[/]", "/"),
                        (r"/_sca_tempreg_maps_zstat_files_smooth_(\d)+[/]", "/"),
                        (r"/qc___", '/qc/')
                    ]

                    node, out_file = rp[resource]
                    workflow.connect(node, out_file, ds, resource)

                    output_sink_nodes += [(ds, 'out_file')]


            if 1 in c.runSymbolicLinks and not ndmg_out and \
                not c.outputDirectory.lower().startswith('s3://'):

                merge_link_node = pe.Node(
                    interface=Merge(len(output_sink_nodes)),
                    name='create_symlinks_paths_{}'.format(num_strat)
                )
                merge_link_node.inputs.ravel_inputs = True

                link_node = pe.Node(
                    interface=function.Function(
                        input_names=[
                            'output_dir',
                            'symlink_dir',
                            'pipeline_id',
                            'subject_id',
                            'paths',
                        ],
                        output_names=[],
                        function=create_symlinks,
                        as_module=True
                    ), name='create_symlinks_{}'.format(num_strat)
                )

                link_node.inputs.output_dir = c.outputDirectory
                link_node.inputs.subject_id = subject_id
                link_node.inputs.pipeline_id = 'pipeline_%s' % pipeline_id

                for i, (node, node_input) in enumerate(output_sink_nodes):
                    workflow.connect(node, node_input,
                                     merge_link_node, 'in{}'.format(i))

                workflow.connect(merge_link_node, 'out', link_node, 'paths')

            try:
                G = nx.DiGraph()
                strat_name = strat.get_name()
                G.add_edges_from([
                    (strat_name[s], strat_name[s + 1])
                    for s in range(len(strat_name) - 1)
                ])

                dotfilename = os.path.join(log_dir, 'strategy.dot')
                nx.drawing.nx_pydot.write_dot(G, dotfilename)
                format_dot(dotfilename, 'png')
            except:
                logger.warn('Cannot Create the strategy and pipeline '
                            'graph, dot or/and pygraphviz is not installed')


        forks = "\n\nStrategy forks:\n" + \
                "\n".join(["- " + pipe for pipe in sorted(set(pipeline_ids))]) + \
                "\n\n"

        logger.info(forks)


        try:

            pipeline_start_datetime = strftime("%Y-%m-%d %H:%M:%S")

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
            for node_name in workflow.list_node_names():
                node = workflow.get_node(node_name)
                cb_logger.debug(json.dumps({
                    "id": str(node),
                    "hash": node.inputs.get_hashval()[1],
                }))

            # Add status callback function that writes in callback log
            if nipype.__version__ not in ('1.1.2'):
                err_msg = "This version of Nipype may not be compatible with " \
                            "CPAC v%s, please install Nipype version 1.1.2\n" \
                            % (CPAC.__version__)
                logger.error(err_msg)
            else:
                from CPAC.utils.monitoring import log_nodes_cb
                plugin_args['status_callback'] = log_nodes_cb

                
            # Actually run the pipeline now, for the current subject
            workflow.run(plugin=plugin, plugin_args=plugin_args)


            # Dump subject info pickle file to subject log dir
            subject_info['status'] = 'Completed'

            subject_info_file = os.path.join(
                log_dir, 'subject_info_%s.pkl' % subject_id
            )
            with open(subject_info_file, 'wb') as info:
                pickle.dump(subject_info, info)

            for i, _ in enumerate(pipeline_ids):
                for scan in scan_ids:
                    create_log_node(workflow, None, None, i, scan).run()

            if 1 in c.generateQualityControlImages and not ndmg_out:
                for pip_id in pipeline_ids:
                    pipeline_base = os.path.join(c.outputDirectory,
                                                'pipeline_%s' % pip_id)

                    sub_output_dir = os.path.join(pipeline_base, subject_id)
                    qc_output_folder = os.path.join(sub_output_dir, 'qc_html')

                    generate_qc_pages(qc_output_folder,
                                    sub_output_dir,
                                    qc_montage_id_a,
                                    qc_montage_id_s,
                                    qc_plot_id,
                                    qc_hist_id)

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

#         except KeyboardInterrupt as e:

#             raising = e

#             execution_info = """

#     Stopped of subject workflow {workflow}

#     CPAC run stopped:

#         Pipeline configuration: {pipeline}
#         Subject workflow: {workflow}
#         Elapsed run time (minutes): {elapsed}
#         Timing information saved in {log_dir}/cpac_individual_timing_{pipeline}.csv
#         System time of start:      {run_start}

# """
        except:

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

            logger.info(execution_info.format(
                workflow=workflow_name,
                pipeline=c.pipelineName,
                log_dir=c.logDirectory,
                elapsed=(time.time() - pipeline_start_time) / 60,
                run_start=pipeline_start_datetime,
                run_finish=strftime("%Y-%m-%d %H:%M:%S")
            ))

            # Remove working directory when done
            if c.removeWorkingDir:
                try:
                    subject_wd = os.path.join(c.workingDirectory, workflow_name)
                    if os.path.exists(subject_wd):
                        logger.info("Removing working dir: %s" % subject_wd)
                        shutil.rmtree(subject_wd)
                except:
                    logger.warn('Could not remove subjects %s working directory',
                                workflow_name)

            # if raising:
            #     raise raising

    return workflow
