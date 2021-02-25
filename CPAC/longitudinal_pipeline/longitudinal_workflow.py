# -*- coding: utf-8 -*-
import os
import copy
import time
import shutil

from nipype import config
from nipype import logging
from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
from nipype.interfaces.utility import Merge, IdentityInterface
import nipype.interfaces.utility as util

from indi_aws import aws_utils

from CPAC.utils.utils import concat_list
from CPAC.utils.interfaces.datasink import DataSink
from CPAC.utils.interfaces.function import Function

import CPAC

from CPAC.pipeline.cpac_pipeline import initialize_nipype_wf, \
    connect_pipeline, build_anat_preproc_stack
from CPAC.pipeline.engine import initiate_rpool

from CPAC.registration import (
    create_fsl_flirt_linear_reg,
    create_fsl_fnirt_nonlinear_reg,
    create_register_func_to_anat,
    create_bbregister_func_to_anat,
    create_wf_calculate_ants_warp
)

from CPAC.registration.utils import run_ants_apply_warp

from CPAC.utils.datasource import (
    resolve_resolution,
    create_anat_datasource,
    create_check_for_s3_node
)

from CPAC.longitudinal_pipeline.longitudinal_preproc import (
    subject_specific_template
)

from CPAC.utils import Strategy, find_files, function, Outputs

from CPAC.utils.utils import (
    check_config_resources,
    check_system_deps,
    get_scan_params,
    get_tr
)

logger = logging.getLogger('nipype.workflow')


def mask_T1w_longitudinal_template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "mask_T1w_longitudinal_template",
     "config": ["longitudinal_template_generation"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["desc-brain_T1w"],
     "outputs": ["space-T1w_desc-brain_mask"]}
    '''

    brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                         name=f'longitudinal_anatomical_brain_mask_'
                              f'{pipe_num}')
    brain_mask.inputs.args = '-bin'

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, brain_mask, 'in_file')

    outputs = {
        'space-T1w_desc-brain_mask': (brain_mask, 'out_file')
    }

    return (wf, outputs)


def create_datasink(datasink_name, config, subject_id, session_id='',
                    strat_name='', map_node_iterfield=None):
    """

    Parameters
    ----------
    datasink_name
    config
    subject_id
    session_id
    strat_name
    map_node_iterfield

    Returns
    -------

    """
    try:
        encrypt_data = bool(
            config.pipeline_setup['Amazon-AWS']['s3_encryption'])
    except:
        encrypt_data = False

    # TODO Enforce value with schema validation
    # Extract credentials path for output if it exists
    try:
        # Get path to creds file
        creds_path = ''
        if config.pipeline_setup['Amazon-AWS'][
            'aws_output_bucket_credentials']:
            creds_path = str(config.pipeline_setup['Amazon-AWS'][
                                 'aws_output_bucket_credentials'])
            creds_path = os.path.abspath(creds_path)

        if config.pipeline_setup['output_directory'][
            'path'].lower().startswith('s3://'):
            # Test for s3 write access
            s3_write_access = \
                aws_utils.test_bucket_access(creds_path,
                                             config.pipeline_setup[
                                                 'output_directory']['path'])

            if not s3_write_access:
                raise Exception('Not able to write to bucket!')

    except Exception as e:
        if config.pipeline_setup['output_directory'][
            'path'].lower().startswith('s3://'):
            err_msg = 'There was an error processing credentials or ' \
                      'accessing the S3 bucket. Check and try again.\n' \
                      'Error: %s' % e
            raise Exception(err_msg)

    if map_node_iterfield is not None:
        ds = pe.MapNode(
            DataSink(infields=map_node_iterfield),
            name='sinker_{}'.format(datasink_name),
            iterfield=map_node_iterfield
        )
    else:
        ds = pe.Node(
            DataSink(),
            name='sinker_{}'.format(datasink_name)
        )

    ds.inputs.base_directory = config.pipeline_setup['output_directory'][
        'path']
    ds.inputs.creds_path = creds_path
    ds.inputs.encrypt_bucket_keys = encrypt_data
    ds.inputs.container = os.path.join(
        'pipeline_%s_%s' % (
        config.pipeline_setup['pipeline_name'], strat_name),
        subject_id, session_id
    )
    return ds


def connect_anat_preproc_inputs(strat, anat_preproc, strat_name,
                                strat_nodes_list_list, workflow):
    """
    Parameters
    ----------
    strat : Strategy
        the strategy object you want to fork
    anat_preproc : Workflow
        the anat_preproc workflow node to be connected and added to the resource pool
    strat_name : str
        name of the strategy
    strat_nodes_list_list : list
        a list of strat_nodes_list
    workflow : Workflow
        main longitudinal workflow

    Returns
    -------
    new_strat : Strategy
        the fork of strat with the resource pool updated
    strat_nodes_list_list : list
        a list of strat_nodes_list
    """

    new_strat = strat.fork()

    tmp_node, out_key = new_strat['anatomical']
    workflow.connect(tmp_node, out_key, anat_preproc, 'inputspec.anat')

    tmp_node, out_key = new_strat['template_cmass']
    workflow.connect(tmp_node, out_key, anat_preproc,
                     'inputspec.template_cmass')

    new_strat.append_name(anat_preproc.name)

    new_strat.update_resource_pool({
        'anatomical_brain': (
            anat_preproc, 'outputspec.brain'),
        'anatomical_skull_leaf': (
            anat_preproc, 'outputspec.reorient'),
        'anatomical_brain_mask': (
            anat_preproc, 'outputspec.brain_mask'),
    })

    try:
        strat_nodes_list_list[strat_name].append(new_strat)
    except KeyError:
        strat_nodes_list_list[strat_name] = [new_strat]

    return new_strat, strat_nodes_list_list


def pick_map(file_list, index, file_type):
    if isinstance(file_list, list):
        if len(file_list) == 1:
            file_list = file_list[0]
        for file_name in file_list:
            if file_name.endswith(f"{file_type}_{index}.nii.gz"):
                return file_name
    return None


def anat_longitudinal_wf(subject_id, sub_list, config):
    """
    Parameters
    ----------
    subject_id : str
        the id of the subject
    sub_list : list of dict
        this is a list of sessions for one subject and each session if the same dictionary as the one given to
        prep_workflow
    config : configuration
        a configuration object containing the information of the pipeline config. (Same as for prep_workflow)

    Returns
    -------
        None
    """

    # list of lists for every strategy
    session_id_list = []
    session_wfs = {}

    # Loop over the sessions to create the input for the longitudinal
    # algorithm
    for session in sub_list:

        unique_id = session['unique_id']
        session_id_list.append(unique_id)

        try:
            creds_path = session['creds_path']
            if creds_path and 'none' not in creds_path.lower():
                if os.path.exists(creds_path):
                    input_creds_path = os.path.abspath(creds_path)
                else:
                    err_msg = 'Credentials path: "%s" for subject "%s" ' \
                              'session "%s" was not found. Check this path ' \
                              'and try again.' % (creds_path, subject_id,
                                                  unique_id)
                    raise Exception(err_msg)
            else:
                input_creds_path = None
        except KeyError:
            input_creds_path = None

        workflow = initialize_nipype_wf(config, sub_list[0],
                                        # just grab the first one for the name
                                        name="anat_longitudinal_template")

        workflow, rpool = initiate_rpool(workflow, config, session)
        pipeline_blocks = build_anat_preproc_stack(rpool, config)
        workflow = connect_pipeline(workflow, config, rpool, pipeline_blocks)

        session_wfs[unique_id] = (workflow, rpool)

        rpool.gather_pipes(workflow, config)

        workflow.run()

    # Here we have all the anat_preproc set up for every session
    # loop over the different anat preproc strategies

    # - ingress output dir for this one subject id
    # - merge nodes to get the desc-brain and desc-preprocs into one node
    #   as a list
    # - run flirt to create the 'template' (orig space tho)
    # - register this image to template (use node blocks now?)
    #    and also run segmentation! I guess on the long-temp
    # - put the xfm's into each session's output dir, and the tissue segs
    #   too
    #
    # now we're back to one sub at a time
    # - warp the preprocced stuff to template space, but using the
    #   xfm's from the long-temp-to-standard reg!
    # - warp the tissue segs to native space
    # - put all back into output dir

    '''
    # just grab one rpool, they should all have the same pipe_idx's/strats
    strats = rpool.get_pipe_idxs('desc-brain_T1w')

    for strat in strats: #strat_name, strat_nodes_list in strat_nodes_list_list.items():

        pipe_num = rpool.get_pipe_number(strat)

        strat_nodes_list = []
        for unique_id, entry in session_wfs.items():
            session_rpool = entry[1]
            strat_nodes_list.append(session_rpool[strat])

        # Merge node to feed the anat_preproc outputs to the longitudinal
        # template generation
        brain_merge_node = pe.Node(
            interface=Merge(len(strat_nodes_list)),
            name="anat_longitudinal_brain_merge_" + pipe_num)

        skull_merge_node = pe.Node(
            interface=Merge(len(strat_nodes_list)),
            name="anat_longitudinal_skull_merge_" + pipe_num)

        # This node will generate the longitudinal template (the functions are in longitudinal_preproc)
        # Later other algorithms could be added to calculate it, like the multivariate template from ANTS
        # It would just require to change it here.
        template_node = subject_specific_template(
            workflow_name='subject_specific_anat_template_' + pipe_num
        )

        template_node.inputs.set(
            avg_method=config.longitudinal_template_generation[
                'average_method'],
            dof=config.longitudinal_template_generation['dof'],
            interp=config.longitudinal_template_generation['interp'],
            cost=config.longitudinal_template_generation['cost'],
            convergence_threshold=config.longitudinal_template_generation[
                'convergence_threshold'],
            thread_pool=config.longitudinal_template_generation[
                'thread_pool'],
            unique_id_list=session_wfs.keys()
        )

        workflow.connect(brain_merge_node, 'out', template_node,
                         'input_brain_list')
        workflow.connect(skull_merge_node, 'out', template_node,
                         'input_skull_list')

        reg_strat_list = register_anat_longitudinal_template_to_standard(
            template_node, config, workflow, strat_init, strat_name)

        # Register T1 to the standard template
        # TODO add session information in node name
        for num_reg_strat, reg_strat in enumerate(reg_strat_list):

            if reg_strat.get('registration_method') == 'FSL':

                fsl_apply_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                                            name='fsl_apply_warp_anat_longitudinal_to_standard_{0}_'.format(
                                                strat_name),
                                            iterfield=['in_file'])

                workflow.connect(template_node, "output_brain_list",
                                 fsl_apply_warp, 'in_file')

                node, out_file = reg_strat['template_brain_for_anat']
                workflow.connect(node, out_file,
                                 fsl_apply_warp, 'ref_file')

                # TODO how to include linear xfm?
                # node, out_file = reg_strat['anatomical_to_mni_linear_xfm']
                # workflow.connect(node, out_file, fsl_apply_warp, 'premat')

                node, out_file = reg_strat['anatomical_to_mni_nonlinear_xfm']
                workflow.connect(node, out_file,
                                 fsl_apply_warp, 'field_file')

                reg_strat.update_resource_pool({
                    'anatomical_to_standard': (fsl_apply_warp, 'out_file')
                })

            elif reg_strat.get('registration_method') == 'ANTS':

                ants_apply_warp = pe.MapNode(
                    util.Function(input_names=['moving_image',
                                               'reference',
                                               'initial',
                                               'rigid',
                                               'affine',
                                               'nonlinear',
                                               'interp'],
                                  output_names=['out_image'],
                                  function=run_ants_apply_warp),
                    name='ants_apply_warp_anat_longitudinal_to_standard_{0}_'.format(
                        strat_name),
                    iterfield=['moving_image'])

                workflow.connect(template_node, "output_brain_list",
                                 ants_apply_warp, 'moving_image')

                node, out_file = reg_strat['template_brain_for_anat']
                workflow.connect(node, out_file, ants_apply_warp, 'reference')

                node, out_file = reg_strat['ants_initial_xfm']
                workflow.connect(node, out_file, ants_apply_warp, 'initial')

                node, out_file = reg_strat['ants_rigid_xfm']
                workflow.connect(node, out_file, ants_apply_warp, 'rigid')

                node, out_file = reg_strat['ants_affine_xfm']
                workflow.connect(node, out_file, ants_apply_warp, 'affine')

                node, out_file = reg_strat['anatomical_to_mni_nonlinear_xfm']
                workflow.connect(node, out_file, ants_apply_warp, 'nonlinear')

                ants_apply_warp.inputs.interp = config.anatRegANTSinterpolation

                reg_strat.update_resource_pool({
                    'anatomical_to_standard': (ants_apply_warp, 'out_image')
                })

        # Register tissue segmentation from longitudinal template space to native space
        fsl_convert_xfm = pe.MapNode(interface=fsl.ConvertXFM(),
                                     name=f'fsl_xfm_longitudinal_to_native_{strat_name}',
                                     iterfield=['in_file'])

        fsl_convert_xfm.inputs.invert_xfm = True

        workflow.connect(template_node, "warp_list",
                         fsl_convert_xfm, 'in_file')

        def seg_apply_warp(strat_name, resource, type='str', file_type=None):

            if type == 'str':

                fsl_apply_xfm = pe.MapNode(interface=fsl.ApplyXFM(),
                                           name=f'fsl_apply_xfm_longitudinal_to_native_{resource}_{strat_name}',
                                           iterfield=['reference',
                                                      'in_matrix_file'])

                fsl_apply_xfm.inputs.interp = 'nearestneighbour'

                node, out_file = reg_strat[resource]
                workflow.connect(node, out_file,
                                 fsl_apply_xfm, 'in_file')

                workflow.connect(brain_merge_node, 'out',
                                 fsl_apply_xfm, 'reference')

                workflow.connect(fsl_convert_xfm, "out_file",
                                 fsl_apply_xfm, 'in_matrix_file')

                reg_strat.update_resource_pool({
                    resource: (fsl_apply_xfm, 'out_file')
                }, override=True)

            elif type == 'list':

                for index in range(3):

                    fsl_apply_xfm = pe.MapNode(interface=fsl.ApplyXFM(),
                                               name=f'fsl_apply_xfm_longitudinal_to_native_{resource}_{index}_{strat_name}',
                                               iterfield=['reference',
                                                          'in_matrix_file'])

                    fsl_apply_xfm.inputs.interp = 'nearestneighbour'

                    pick_seg_map = pe.Node(Function(
                        input_names=['file_list', 'index', 'file_type'],
                        output_names=['file_name'],
                        function=pick_map),
                                           name=f'pick_{file_type}_{index}_{strat_name}')

                    node, out_file = reg_strat[resource]

                    workflow.connect(node, out_file,
                                     pick_seg_map, 'file_list')

                    pick_seg_map.inputs.index = index
                    pick_seg_map.inputs.file_type = file_type

                    workflow.connect(pick_seg_map, 'file_name',
                                     fsl_apply_xfm, 'in_file')

                    workflow.connect(brain_merge_node, 'out',
                                     fsl_apply_xfm, 'reference')

                    workflow.connect(fsl_convert_xfm, 'out_file',
                                     fsl_apply_xfm, 'in_matrix_file')

                    concat_seg_map = pe.Node(
                        Function(input_names=['in_list1', 'in_list2'],
                                 output_names=['out_list'],
                                 function=concat_list),
                        name=f'concat_{file_type}_{index}_{strat_name}')

                    if index == 0:
                        workflow.connect(fsl_apply_xfm, 'out_file',
                                         concat_seg_map, 'in_list1')

                        reg_strat.update_resource_pool({
                            f'temporary_{resource}_list': (
                            concat_seg_map, 'out_list')
                        })

                    else:
                        workflow.connect(fsl_apply_xfm, 'out_file',
                                         concat_seg_map, 'in_list2')

                        node, out_file = reg_strat[
                            f'temporary_{resource}_list']

                        workflow.connect(node, out_file,
                                         concat_seg_map, 'in_list1')

                        reg_strat.update_resource_pool({
                            f'temporary_{resource}_list': (
                            concat_seg_map, 'out_list')
                        }, override=True)

                reg_strat.update_resource_pool({
                    resource: (concat_seg_map, 'out_list')
                }, override=True)

        for seg in ['anatomical_gm_mask', 'anatomical_csf_mask',
                    'anatomical_wm_mask',
                    'seg_mixeltype', 'seg_partial_volume_map']:
            seg_apply_warp(strat_name=strat_name, resource=seg)

        # apply warp on list
        seg_apply_warp(strat_name=strat_name, resource='seg_probability_maps',
                       type='list', file_type='prob')

        seg_apply_warp(strat_name=strat_name,
                       resource='seg_partial_volume_files', type='list',
                       file_type='pve')
    '''


# TODO check:
# 1 func alone works
# 2 anat + func works, pass anat strategy list?
def func_preproc_longitudinal_wf(subject_id, sub_list, config):
    """
    Parameters
    ----------
    subject_id : string
        the id of the subject
    sub_list : list of dict
        this is a list of sessions for one subject and each session if the same dictionary as the one given to
        prep_workflow
    config : configuration
        a configuration object containing the information of the pipeline config. (Same as for prep_workflow)

    Returns
    -------
    strat_list_ses_list : list of list
        a list of strategies; within each strategy, a list of sessions
    """

    datasink = pe.Node(nio.DataSink(), name='sinker')
    datasink.inputs.base_directory = \
    config.pipeline_setup['working_directory']['path']

    session_id_list = []
    ses_list_strat_list = {}

    workflow_name = 'func_preproc_longitudinal_' + str(subject_id)
    workflow = pe.Workflow(name=workflow_name)
    workflow.base_dir = config.pipeline_setup['working_directory']['path']
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(
            config.pipeline_setup['crash_directory']['path'])
    }

    for sub_dict in sub_list:
        if 'func' in sub_dict or 'rest' in sub_dict:
            if 'func' in sub_dict:
                func_paths_dict = sub_dict['func']
            else:
                func_paths_dict = sub_dict['rest']

            unique_id = sub_dict['unique_id']
            session_id_list.append(unique_id)

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

            strat = Strategy()
            strat_list = [strat]
            node_suffix = '_'.join([subject_id, unique_id])

            # Functional Ingress Workflow
            # add optional flag
            workflow, diff, blip, fmap_rp_list = connect_func_ingress(
                workflow,
                strat_list,
                config,
                sub_dict,
                subject_id,
                input_creds_path,
                node_suffix)

            # Functional Initial Prep Workflow
            workflow, strat_list = connect_func_init(workflow, strat_list,
                                                     config, node_suffix)

            # Functional Image Preprocessing Workflow
            workflow, strat_list = connect_func_preproc(workflow, strat_list,
                                                        config, node_suffix)

            # Distortion Correction
            workflow, strat_list = connect_distortion_correction(workflow,
                                                                 strat_list,
                                                                 config,
                                                                 diff,
                                                                 blip,
                                                                 fmap_rp_list,
                                                                 node_suffix)

            ses_list_strat_list[node_suffix] = strat_list

    # Here we have all the func_preproc set up for every session of the subject

    # TODO create a list of list ses_list_strat_list
    # a list of skullstripping strategies,
    # a list of sessions within each strategy list
    # TODO rename and reorganize dict
    # TODO update strat name
    strat_list_ses_list = {}
    strat_list_ses_list['func_default'] = []

    for sub_ses_id, strat_nodes_list in ses_list_strat_list.items():
        strat_list_ses_list['func_default'].append(strat_nodes_list[0])

    workflow.run()

    return strat_list_ses_list


def merge_func_preproc(working_directory):
    """
    Parameters
    ----------
    working_directory : string
        a path to the working directory

    Returns
    -------
    brain_list : list
        a list of func preprocessed brain
    skull_list : list
        a list of func preprocessed skull
    """

    brain_list = []
    skull_list = []

    for dirpath, dirnames, filenames in os.walk(working_directory):
        for f in filenames:
            if 'func_get_preprocessed_median' in dirpath and '.nii.gz' in f:
                filepath = os.path.join(dirpath, f)
                brain_list.append(filepath)
            if 'func_get_motion_correct_median' in dirpath and '.nii.gz' in f:
                filepath = os.path.join(dirpath, f)
                skull_list.append(filepath)

    brain_list.sort()
    skull_list.sort()

    return brain_list, skull_list


def register_func_longitudinal_template_to_standard(
        longitudinal_template_node, c, workflow, strat_init, strat_name):
    sub_mem_gb, num_cores_per_sub, num_ants_cores, num_omp_cores = \
        check_config_resources(c)

    strat_init_new = strat_init.fork()

    strat_init_new.update_resource_pool({
        'functional_preprocessed_median': (
        longitudinal_template_node, 'brain_template'),
        'motion_correct_median': (
        longitudinal_template_node, 'skull_template')
    })

    strat_list = [strat_init_new]

    new_strat_list = []

    regOption = c.anatomical_preproc[
        'registration_workflow'
    ]['registration']['using']

    if 'FSL' in regOption:

        for num_strat, strat in enumerate(strat_list):

            flirt_reg_func_mni = create_fsl_flirt_linear_reg(
                'func_mni_flirt_register_%s_%d' % (strat_name, num_strat)
            )

            if c.functional_registration['2-func_registration_to_template'][
                'FNIRT_pipelines']['interpolation'] not in ["trilinear",
                                                            "sinc", "spline"]:
                err_msg = 'The selected FSL interpolation method may be in the list of values: "trilinear", "sinc", "spline"'
                raise Exception(err_msg)

            # Input registration parameters
            flirt_reg_func_mni.inputs.inputspec.interp = \
            c.functional_registration['2-func_registration_to_template'][
                'FNIRT_pipelines']['interpolation']

            node, out_file = strat['functional_preprocessed_median']
            workflow.connect(node, out_file,
                             flirt_reg_func_mni, 'inputspec.input_brain')

            # pass the reference files
            node, out_file = strat['template_brain_for_func_preproc']
            workflow.connect(node, out_file, flirt_reg_func_mni,
                             'inputspec.reference_brain')

            if 'ANTS' in regOption:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(flirt_reg_func_mni.name)

            strat.update_resource_pool({
                'registration_method': 'FSL',
                'func_longitudinal_to_mni_linear_xfm': (
                flirt_reg_func_mni, 'outputspec.linear_xfm'),
                'mni_to_func_longitudinal_linear_xfm': (
                flirt_reg_func_mni, 'outputspec.invlinear_xfm'),
                'func_longitudinal_template_to_standard': (
                flirt_reg_func_mni, 'outputspec.output_brain')
            })

    strat_list += new_strat_list

    new_strat_list = []

    try:
        fsl_linear_reg_only = c.fsl_linear_reg_only
    except AttributeError:
        fsl_linear_reg_only = [0]

    if 'FSL' in regOption and 0 in fsl_linear_reg_only:

        for num_strat, strat in enumerate(strat_list):

            if strat.get('registration_method') == 'FSL':

                fnirt_reg_func_mni = create_fsl_fnirt_nonlinear_reg(
                    'func_mni_fnirt_register_%s_%d' % (strat_name, num_strat)
                )

                # brain input
                node, out_file = strat['functional_preprocessed_median']
                workflow.connect(node, out_file,
                                 fnirt_reg_func_mni, 'inputspec.input_brain')

                # brain reference
                node, out_file = strat['template_brain_for_func_preproc']
                workflow.connect(node, out_file,
                                 fnirt_reg_func_mni,
                                 'inputspec.reference_brain')

                # skull input
                node, out_file = strat['motion_correct_median']
                workflow.connect(node, out_file,
                                 fnirt_reg_func_mni, 'inputspec.input_skull')

                # skull reference
                node, out_file = strat['template_skull_for_func_preproc']
                workflow.connect(node, out_file,
                                 fnirt_reg_func_mni,
                                 'inputspec.reference_skull')

                node, out_file = strat['func_longitudinal_to_mni_linear_xfm']
                workflow.connect(node, out_file,
                                 fnirt_reg_func_mni, 'inputspec.linear_aff')

                node, out_file = strat['template_ref_mask']
                workflow.connect(node, out_file,
                                 fnirt_reg_func_mni, 'inputspec.ref_mask')

                # assign the FSL FNIRT config file specified in pipeline
                # config.yml
                fnirt_reg_func_mni.inputs.inputspec.fnirt_config = \
                c.anatomical_preproc['registration_workflow']['registration'][
                    'FSL-FNIRT']['fnirt_config']

                if 1 in fsl_linear_reg_only:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(fnirt_reg_func_mni.name)

                strat.update_resource_pool({
                    'func_longitudinal_to_mni_nonlinear_xfm': (
                    fnirt_reg_func_mni, 'outputspec.nonlinear_xfm'),
                    'func_longitudinal_template_to_standard': (
                    fnirt_reg_func_mni, 'outputspec.output_brain')
                }, override=True)

    strat_list += new_strat_list

    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        # or run ANTS anatomical-to-MNI registration instead
        if 'ANTS' in regOption and \
                        strat.get('registration_method') != 'FSL':

            ants_reg_func_mni = \
                create_wf_calculate_ants_warp(
                    'func_mni_ants_register_%s_%d' % (strat_name, num_strat),
                    num_threads=num_ants_cores,
                    reg_ants_skull=
                    c.anatomical_preproc['registration_workflow'][
                        'reg_with_skull']
                )

            if c.functional_registration['2-func_registration_to_template'][
                'ANTs_pipelines']['interpolation'] not in ['Linear',
                                                           'BSpline',
                                                           'LanczosWindowedSinc']:
                err_msg = 'The selected ANTS interpolation method may be in the list of values: "Linear", "BSpline", "LanczosWindowedSinc"'
                raise Exception(err_msg)

            # Input registration parameters
            ants_reg_func_mni.inputs.inputspec.interp = \
            c.functional_registration['2-func_registration_to_template'][
                'ANTs_pipelines']['interpolation']

            # calculating the transform with the skullstripped is
            # reported to be better, but it requires very high
            # quality skullstripping. If skullstripping is imprecise
            # registration with skull is preferred
            if c.anatomical_preproc['registration_workflow'][
                'reg_with_skull']:

                # get the skull-stripped anatomical from resource pool
                node, out_file = strat['functional_preprocessed_median']

                # pass the anatomical to the workflow
                workflow.connect(node, out_file,
                                 ants_reg_func_mni, 'inputspec.moving_brain')

                # get the reorient skull-on anatomical from resource pool
                node, out_file = strat['motion_correct_median']

                # pass the anatomical to the workflow
                workflow.connect(node, out_file,
                                 ants_reg_func_mni, 'inputspec.moving_skull')

                # pass the reference file
                node, out_file = strat['template_brain_for_func_preproc']
                workflow.connect(node, out_file,
                                 ants_reg_func_mni,
                                 'inputspec.reference_brain')

                # pass the reference file
                node, out_file = strat['template_skull_for_func_preproc']
                workflow.connect(node, out_file,
                                 ants_reg_func_mni,
                                 'inputspec.reference_skull')

            else:

                node, out_file = strat['functional_preprocessed_median']

                workflow.connect(node, out_file,
                                 ants_reg_func_mni, 'inputspec.moving_brain')

                # pass the reference file
                node, out_file = strat['template_brain_for_func_preproc']
                workflow.connect(node, out_file,
                                 ants_reg_func_mni,
                                 'inputspec.reference_brain')

            # pass the reference mask file
            node, out_file = strat['template_brain_mask_for_func_preproc']
            workflow.connect(
                node, out_file,
                ants_reg_func_mni, 'inputspec.reference_mask'
            )

            # pass the reference mask file
            node, out_file = strat['functional_brain_mask']
            workflow.connect(
                node, out_file,
                ants_reg_func_mni, 'inputspec.moving_mask'
            )

            ants_reg_func_mni.inputs.inputspec.ants_para = \
            c.anatomical_preproc['registration_workflow']['registration'][
                'ANTs']['T1_registration']
            ants_reg_func_mni.inputs.inputspec.fixed_image_mask = None

            strat.append_name(ants_reg_func_mni.name)

            strat.update_resource_pool({
                'registration_method': 'ANTS',
                'ants_initial_xfm': (
                ants_reg_func_mni, 'outputspec.ants_initial_xfm'),
                'ants_rigid_xfm': (
                ants_reg_func_mni, 'outputspec.ants_rigid_xfm'),
                'ants_affine_xfm': (
                ants_reg_func_mni, 'outputspec.ants_affine_xfm'),
                'func_longitudinal_to_mni_nonlinear_xfm': (
                ants_reg_func_mni, 'outputspec.warp_field'),
                'mni_to_func_longitudinal_nonlinear_xfm': (
                ants_reg_func_mni, 'outputspec.inverse_warp_field'),
                'func_longitudinal_to_mni_ants_composite_xfm': (
                ants_reg_func_mni, 'outputspec.composite_transform'),
                'func_longitudinal_template_to_standard': (
                ants_reg_func_mni, 'outputspec.normalized_output_brain')
            })

    strat_list += new_strat_list

    '''
    # Func -> T1 Registration (Initial Linear Reg)
    workflow, strat_list, diff_complete = connect_func_to_anat_init_reg(workflow, strat_list, c)

    # Func -> T1 Registration (BBREG)
    workflow, strat_list = connect_func_to_anat_bbreg(workflow, strat_list, c, diff_complete)

    # Func -> T1/EPI Template
    workflow, strat_list = connect_func_to_template_reg(workflow, strat_list, c)
    '''

    return workflow, strat_list


def func_longitudinal_template_wf(subject_id, strat_list, config):
    '''
    Parameters
    ----------
    subject_id : string
        the id of the subject
    strat_list : list of list
        first level strategy, second level session
    config : configuration
        a configuration object containing the information of the pipeline config.

    Returns
    -------
        None
    '''

    workflow_name = 'func_longitudinal_template_' + str(subject_id)
    workflow = pe.Workflow(name=workflow_name)
    workflow.base_dir = config.pipeline_setup['working_directory']['path']
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(
            config.pipeline_setup['crash_directory']['path'])
    }

    # strat_nodes_list = strat_list['func_default']
    strat_init = Strategy()

    templates_for_resampling = [
        (config.resolution_for_func_preproc,
         config.template_brain_only_for_func,
         'template_brain_for_func_preproc', 'resolution_for_func_preproc'),
        (config.resolution_for_func_preproc, config.template_skull_for_func,
         'template_skull_for_func_preproc', 'resolution_for_func_preproc'),
        (config.resolution_for_func_preproc, config.ref_mask_for_func,
         'template_ref_mask', 'resolution_for_func_preproc'),
        # TODO check float resolution
        (config.resolution_for_func_preproc,
         config.functional_registration['2-func_registration_to_template'][
             'target_template']['EPI_template']['template_epi'],
         'template_epi', 'resolution_for_func_preproc'),
        (config.resolution_for_func_derivative,
         config.functional_registration['2-func_registration_to_template'][
             'target_template']['EPI_template']['template_epi'],
         'template_epi_derivative', 'resolution_for_func_derivative'),
        (config.resolution_for_func_derivative,
         config.template_brain_only_for_func,
         'template_brain_for_func_derivative', 'resolution_for_func_preproc'),
        (
        config.resolution_for_func_derivative, config.template_skull_for_func,
        'template_skull_for_func_derivative', 'resolution_for_func_preproc'),
    ]

    for resolution, template, template_name, tag in templates_for_resampling:
        resampled_template = pe.Node(Function(
            input_names=['resolution', 'template', 'template_name', 'tag'],
            output_names=['resampled_template'],
            function=resolve_resolution,
            as_module=True),
                                     name='resampled_' + template_name)

        resampled_template.inputs.resolution = resolution
        resampled_template.inputs.template = template
        resampled_template.inputs.template_name = template_name
        resampled_template.inputs.tag = tag

        strat_init.update_resource_pool({
            template_name: (resampled_template, 'resampled_template')
        })

    merge_func_preproc_node = pe.Node(
        Function(input_names=['working_directory'],
                 output_names=['brain_list', 'skull_list'],
                 function=merge_func_preproc,
                 as_module=True),
        name='merge_func_preproc')

    merge_func_preproc_node.inputs.working_directory = \
    config.pipeline_setup['working_directory']['path']

    template_node = subject_specific_template(
        workflow_name='subject_specific_func_template_' + subject_id
    )

    template_node.inputs.set(
        avg_method=config.longitudinal_template_average_method,
        dof=config.longitudinal_template_dof,
        interp=config.longitudinal_template_interp,
        cost=config.longitudinal_template_cost,
        convergence_threshold=config.longitudinal_template_convergence_threshold,
        thread_pool=config.longitudinal_template_thread_pool,
    )

    workflow.connect(merge_func_preproc_node, 'brain_list',
                     template_node, 'input_brain_list')

    workflow.connect(merge_func_preproc_node, 'skull_list',
                     template_node, 'input_skull_list')

    workflow, strat_list = register_func_longitudinal_template_to_standard(
        template_node,
        config,
        workflow,
        strat_init,
        'default'
    )

    workflow.run()

    return
