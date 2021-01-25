# -*- coding: utf-8 -*-
import os
import copy
import time
import shutil

from nipype import config
from nipype import logging
import nipype.pipeline.engine as pe
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


def register_anat_longitudinal_template_to_standard(
        longitudinal_template_node, c, workflow, strat_init, strat_name):
    brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                         name=f'longitudinal_anatomical_brain_mask_{strat_name}')

    brain_mask.inputs.args = '-bin'

    workflow.connect(longitudinal_template_node, 'brain_template',
                     brain_mask, 'in_file')

    strat_init_new = strat_init.fork()

    strat_init_new.update_resource_pool({
        'anatomical_brain': (longitudinal_template_node, 'brain_template'),
        'anatomical_skull_leaf': (
        longitudinal_template_node, 'skull_template'),
        'anatomical_brain_mask': (brain_mask, 'out_file')
    })

    strat_list = [strat_init_new]

    # only need to run once for each subject
    already_skullstripped = c.anatomical_preproc['brain_extraction'][
        'already_skullstripped']
    if already_skullstripped == 2:
        already_skullstripped = 0
    elif already_skullstripped == 3:
        already_skullstripped = 1

    sub_mem_gb, num_cores_per_sub, num_ants_cores, num_omp_cores = \
        check_config_resources(c)

    new_strat_list = []

    regOption = c.anatomical_preproc[
        'registration_workflow'
    ]['registration']['using']

    # either run FSL anatomical-to-MNI registration, or...
    if 'FSL' in regOption:
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
                'anat_mni_flirt_register_%s_%d' % (strat_name, num_strat)
            )

            # Input registration parameters
            flirt_reg_anat_mni.inputs.inputspec.interp = \
            c.anatomical_preproc['registration_workflow']['registration'][
                'FSL-FNIRT']['interpolation']

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file,
                             flirt_reg_anat_mni, 'inputspec.input_brain')

            # pass the reference files
            node, out_file = strat['template_brain_for_anat']
            workflow.connect(node, out_file, flirt_reg_anat_mni,
                             'inputspec.reference_brain')

            if 'ANTS' in regOption:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(flirt_reg_anat_mni.name)

            strat.update_resource_pool({
                'registration_method': 'FSL',
                'anatomical_to_mni_linear_xfm': (
                flirt_reg_anat_mni, 'outputspec.linear_xfm'),
                'mni_to_anatomical_linear_xfm': (
                flirt_reg_anat_mni, 'outputspec.invlinear_xfm'),
                'anat_longitudinal_template_to_standard': (
                flirt_reg_anat_mni, 'outputspec.output_brain')
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

                fnirt_reg_anat_mni = create_fsl_fnirt_nonlinear_reg(
                    'anat_mni_fnirt_register_%s_%d' % (strat_name, num_strat)
                )

                # brain input
                node, out_file = strat['anatomical_brain']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.input_brain')

                # brain reference
                node, out_file = strat['template_brain_for_anat']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni,
                                 'inputspec.reference_brain')

                # skull input
                node, out_file = strat['anatomical_skull_leaf']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.input_skull')

                # skull reference
                node, out_file = strat['template_skull_for_anat']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni,
                                 'inputspec.reference_skull')

                node, out_file = strat['anatomical_to_mni_linear_xfm']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.linear_aff')

                node, out_file = strat['template_ref_mask']
                workflow.connect(node, out_file,
                                 fnirt_reg_anat_mni, 'inputspec.ref_mask')

                # assign the FSL FNIRT config file specified in pipeline
                # config.yml
                fnirt_reg_anat_mni.inputs.inputspec.fnirt_config = \
                c.anatomical_preproc['registration_workflow']['registration'][
                    'FSL-FNIRT']['fnirt_config']

                if 1 in fsl_linear_reg_only:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(fnirt_reg_anat_mni.name)

                strat.update_resource_pool({
                    'anatomical_to_mni_nonlinear_xfm': (
                    fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm'),
                    'anat_longitudinal_template_to_standard': (
                    fnirt_reg_anat_mni, 'outputspec.output_brain')
                }, override=True)

    strat_list += new_strat_list

    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        # or run ANTS anatomical-to-MNI registration instead
        if 'ANTS' in regOption and \
                        strat.get('registration_method') != 'FSL':

            ants_reg_anat_mni = \
                create_wf_calculate_ants_warp(
                    'anat_mni_ants_register_%s_%d' % (strat_name, num_strat),
                    num_threads=num_ants_cores,
                    reg_ants_skull=
                    c.anatomical_preproc['registration_workflow'][
                        'reg_with_skull']
                )

            # if someone doesn't have anatRegANTSinterpolation in their pipe config,
            # it will default to LanczosWindowedSinc
            if not hasattr(c, 'anatRegANTSinterpolation'):
                setattr(c, 'anatRegANTSinterpolation', 'LanczosWindowedSinc')

            if c.anatRegANTSinterpolation not in ['Linear', 'BSpline',
                                                  'LanczosWindowedSinc']:
                err_msg = 'The selected ANTS interpolation method may be in the list of values: "Linear", "BSpline", "LanczosWindowedSinc"'
                raise Exception(err_msg)

            # Input registration parameters
            ants_reg_anat_mni.inputs.inputspec.interp = c.anatRegANTSinterpolation

            # calculating the transform with the skullstripped is
            # reported to be better, but it requires very high
            # quality skullstripping. If skullstripping is imprecise
            # registration with skull is preferred
            if c.anatomical_preproc['registration_workflow'][
                'reg_with_skull']:

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
                                 ants_reg_anat_mni, 'inputspec.moving_brain')

                # get the reorient skull-on anatomical from resource pool
                node, out_file = strat['anatomical_skull_leaf']

                # pass the anatomical to the workflow
                workflow.connect(node, out_file,
                                 ants_reg_anat_mni, 'inputspec.moving_skull')

                # pass the reference file
                node, out_file = strat['template_brain_for_anat']
                workflow.connect(node, out_file,
                                 ants_reg_anat_mni,
                                 'inputspec.reference_brain')

                # pass the reference file
                node, out_file = strat['template_skull_for_anat']
                workflow.connect(node, out_file,
                                 ants_reg_anat_mni,
                                 'inputspec.reference_skull')

            else:

                node, out_file = strat['anatomical_brain']

                workflow.connect(node, out_file,
                                 ants_reg_anat_mni, 'inputspec.moving_brain')

                # pass the reference file
                node, out_file = strat['template_brain_for_anat']
                workflow.connect(node, out_file,
                                 ants_reg_anat_mni,
                                 'inputspec.reference_brain')

            # pass the reference mask file
            node, out_file = strat['template_brain_mask_for_anat']
            workflow.connect(
                node, out_file,
                ants_reg_anat_mni, 'inputspec.reference_mask'
            )

            # pass the reference mask file
            node, out_file = strat['anatomical_brain_mask']
            workflow.connect(
                node, out_file,
                ants_reg_anat_mni, 'inputspec.moving_mask'
            )

            ants_reg_anat_mni.inputs.inputspec.ants_para = \
            c.anatomical_preproc['registration_workflow']['registration'][
                'ANTs']['T1_registration']
            ants_reg_anat_mni.inputs.inputspec.fixed_image_mask = None

            strat.append_name(ants_reg_anat_mni.name)

            strat.update_resource_pool({
                'registration_method': 'ANTS',
                'ants_initial_xfm': (
                ants_reg_anat_mni, 'outputspec.ants_initial_xfm'),
                'ants_rigid_xfm': (
                ants_reg_anat_mni, 'outputspec.ants_rigid_xfm'),
                'ants_affine_xfm': (
                ants_reg_anat_mni, 'outputspec.ants_affine_xfm'),
                'anatomical_to_mni_nonlinear_xfm': (
                ants_reg_anat_mni, 'outputspec.warp_field'),
                'mni_to_anatomical_nonlinear_xfm': (
                ants_reg_anat_mni, 'outputspec.inverse_warp_field'),
                'anat_to_mni_ants_composite_xfm': (
                ants_reg_anat_mni, 'outputspec.composite_transform'),
                'anat_longitudinal_template_to_standard': (
                ants_reg_anat_mni, 'outputspec.normalized_output_brain')
            })

    strat_list += new_strat_list

    # [SYMMETRIC] T1 -> Symmetric Template, Non-linear registration (FNIRT/ANTS)

    new_strat_list = []

    if c.voxel_mirrored_homotopic_connectivity['run'] and \
            c.functional_preproc['run']:

        for num_strat, strat in enumerate(strat_list):

            if 'FSL' in regOption and \
                            strat.get('registration_method') != 'ANTS':

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
                    'anat_symmetric_mni_flirt_register_%s_%d' % (
                    strat_name, num_strat)
                )

                flirt_reg_anat_symm_mni.inputs.inputspec.interp = \
                c.anatomical_preproc['registration_workflow']['registration'][
                    'FSL-FNIRT']['interpolation']

                node, out_file = strat['anatomical_brain']
                workflow.connect(node, out_file,
                                 flirt_reg_anat_symm_mni,
                                 'inputspec.input_brain')

                node, out_file = strat['template_symmetric_brain']
                workflow.connect(node, out_file,
                                 flirt_reg_anat_symm_mni,
                                 'inputspec.reference_brain')

                # if 'ANTS' in regOption:
                #    strat = strat.fork()
                #    new_strat_list.append(strat)

                strat.append_name(flirt_reg_anat_symm_mni.name)

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

        if 'FSL' in regOption and 0 in fsl_linear_reg_only:

            for num_strat, strat in enumerate(strat_list):

                if strat.get('registration_method') == 'FSL':
                    fnirt_reg_anat_symm_mni = create_fsl_fnirt_nonlinear_reg(
                        'anat_symmetric_mni_fnirt_register_%s_%d' % (
                        strat_name, num_strat)
                    )

                    node, out_file = strat['anatomical_brain']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.input_brain')

                    node, out_file = strat['anatomical_skull_leaf']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.input_skull')

                    node, out_file = strat['template_brain_for_anat']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.reference_brain')

                    node, out_file = strat['template_symmetric_skull']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.reference_skull')

                    node, out_file = strat['anatomical_to_mni_linear_xfm']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.linear_aff')

                    node, out_file = strat[
                        'template_dilated_symmetric_brain_mask']
                    workflow.connect(node, out_file,
                                     fnirt_reg_anat_symm_mni,
                                     'inputspec.ref_mask')

                    strat.append_name(fnirt_reg_anat_symm_mni.name)

                    strat.update_resource_pool({
                        'anatomical_to_symmetric_mni_nonlinear_xfm': (
                            fnirt_reg_anat_symm_mni,
                            'outputspec.nonlinear_xfm'),
                        'symmetric_anatomical_to_standard': (
                            fnirt_reg_anat_symm_mni,
                            'outputspec.output_brain')
                    }, override=True)

        strat_list += new_strat_list

        new_strat_list = []

        for num_strat, strat in enumerate(strat_list):

            if 'ANTS' in regOption and \
                            strat.get('registration_method') != 'FSL':

                ants_reg_anat_symm_mni = \
                    create_wf_calculate_ants_warp(
                        'anat_symmetric_mni_ants_register_%s_%d' % (
                        strat_name, num_strat),
                        num_threads=num_ants_cores,
                        reg_ants_skull=
                        c.anatomical_preproc['registration_workflow'][
                            'reg_with_skull']
                    )

                # Input registration parameters
                ants_reg_anat_symm_mni.inputs.inputspec.interp = c.anatRegANTSinterpolation

                # calculating the transform with the skullstripped is
                # reported to be better, but it requires very high
                # quality skullstripping. If skullstripping is imprecise
                # registration with skull is preferred
                if 1 in c.anatomical_preproc['registration_workflow'][
                    'reg_with_skull']:

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
                                     ants_reg_anat_symm_mni,
                                     'inputspec.reference_brain')

                    # get the reorient skull-on anatomical from resource pool
                    node, out_file = strat['anatomical_skull_leaf']

                    # pass the anatomical to the workflow
                    workflow.connect(node, out_file,
                                     ants_reg_anat_symm_mni,
                                     'inputspec.moving_skull')

                    # pass the reference file
                    node, out_file = strat['template_symmetric_skull']
                    workflow.connect(node, out_file,
                                     ants_reg_anat_symm_mni,
                                     'inputspec.reference_skull')

                else:
                    # get the skullstripped anatomical from resource pool
                    node, out_file = strat['anatomical_brain']

                    workflow.connect(node, out_file,
                                     ants_reg_anat_symm_mni,
                                     'inputspec.moving_brain')

                    # pass the reference file
                    node, out_file = strat['template_symmetric_brain']
                    workflow.connect(node, out_file,
                                     ants_reg_anat_symm_mni,
                                     'inputspec.reference_brain')

                # pass the reference mask file
                node, out_file = strat['template_brain_mask_for_anat']
                workflow.connect(
                    node, out_file,
                    ants_reg_anat_symm_mni, 'inputspec.reference_mask'
                )

                # pass the reference mask file
                node, out_file = strat['anatomical_brain_mask']
                workflow.connect(
                    node, out_file,
                    ants_reg_anat_symm_mni, 'inputspec.moving_mask'
                )

                ants_reg_anat_symm_mni.inputs.inputspec.ants_para = \
                c.anatomical_preproc['registration_workflow']['registration'][
                    'ANTs']['T1_registration']

                ants_reg_anat_symm_mni.inputs.inputspec.fixed_image_mask = None

                strat.append_name(ants_reg_anat_symm_mni.name)

                strat.update_resource_pool({
                    'ants_symmetric_initial_xfm': (
                    ants_reg_anat_symm_mni, 'outputspec.ants_initial_xfm'),
                    'ants_symmetric_rigid_xfm': (
                    ants_reg_anat_symm_mni, 'outputspec.ants_rigid_xfm'),
                    'ants_symmetric_affine_xfm': (
                    ants_reg_anat_symm_mni, 'outputspec.ants_affine_xfm'),
                    'anatomical_to_symmetric_mni_nonlinear_xfm': (
                    ants_reg_anat_symm_mni, 'outputspec.warp_field'),
                    'symmetric_mni_to_anatomical_nonlinear_xfm': (
                        ants_reg_anat_symm_mni,
                        'outputspec.inverse_warp_field'),
                    'anat_to_symmetric_mni_ants_composite_xfm': (
                        ants_reg_anat_symm_mni,
                        'outputspec.composite_transform'),
                    'symmetric_anatomical_to_standard': (
                    ants_reg_anat_symm_mni,
                    'outputspec.normalized_output_brain')
                })

        strat_list += new_strat_list

    # Inserting Segmentation Preprocessing Workflow
    workflow, strat_list = connect_anat_segmentation(workflow, strat_list, c,
                                                     strat_name)

    return strat_list


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
    workflow: Workflow
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

    workflow = pe.Workflow(
        name="anat_longitudinal_template_" + str(subject_id))
    workflow.base_dir = config.pipeline_setup['working_directory']['path']
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(
            config.pipeline_setup['log_directory']['path'])
    }

    resampled_template = pe.Node(Function(
        input_names=['resolution', 'template', 'template_name', 'tag'],
        output_names=['resampled_template'],
        function=resolve_resolution,
        as_module=True),
                                 name='template_skull_for_anat')
    resampled_template.inputs.resolution = \
    config.anatomical_preproc['registration_workflow']['resolution_for_anat']
    resampled_template.inputs.template = \
    config.anatomical_preproc['registration_workflow'][
        'template_skull_for_anat']
    resampled_template.inputs.template_name = 'template_skull_for_anat'
    resampled_template.inputs.tag = 'resolution_for_anat'

    # Node to calculate the center of mass of the standard template to align the images with it.
    template_center_of_mass = pe.Node(
        interface=afni.CenterMass(),
        name='template_skull_for_anat_center_of_mass'
    )

    template_center_of_mass.inputs.cm_file = "template_center_of_mass.txt"

    workflow.connect(resampled_template, 'resampled_template',
                     template_center_of_mass, 'in_file')

    # list of lists for every strategy
    strat_nodes_list_list = {}

    # list of the data config dictionaries to be updated during the preprocessing
    # creds_list = []

    session_id_list = []

    # Loop over the sessions to create the input for the longitudinal algorithm
    for session in sub_list:

        unique_id = session['unique_id']
        session_id_list.append(unique_id)

        try:
            creds_path = session['creds_path']
            if creds_path and 'none' not in creds_path.lower():
                if os.path.exists(creds_path):
                    input_creds_path = os.path.abspath(creds_path)
                else:
                    err_msg = 'Credentials path: "%s" for subject "%s" session "%s" ' \
                              'was not found. Check this path and try ' \
                              'again.' % (creds_path, subject_id, unique_id)
                    raise Exception(err_msg)
            else:
                input_creds_path = None
        except KeyError:
            input_creds_path = None

        template_keys = [
            ("anat",
             ["anatomical_preproc", "segmentation_workflow", "2-use_priors",
              "CSF_path"]),
            ("anat",
             ["anatomical_preproc", "segmentation_workflow", "2-use_priors",
              "GM_path"]),
            ("anat",
             ["anatomical_preproc", "segmentation_workflow", "2-use_priors",
              "WM_path"]),
            ("anat",
             ["anatomical_preproc", "segmentation_workflow", "1-segmentation",
              "Template_Based", "CSF"]),
            ("anat",
             ["anatomical_preproc", "segmentation_workflow", "1-segmentation",
              "Template_Based", "GRAY"]),
            ("anat",
             ["anatomical_preproc", "segmentation_workflow", "1-segmentation",
              "Template_Based", "WHITE"]),
            ("other", ["voxel_mirrored_homotopic_connectivity",
                       "symmetric_registration", "FNIRT_pipelines",
                       "config_file"]),
        ]

        for key_type, key in template_keys:

            attr = config.get_nested(config, key)

            if isinstance(attr, str) or attr == None:
                node = create_check_for_s3_node(
                    key[-1],
                    attr, key_type,
                    input_creds_path,
                    config.pipeline_setup['working_directory']['path'],
                    map_node=False
                )

                config.set_nested(config, key, node)

        strat = Strategy()
        strat_list = []
        node_suffix = '_'.join([subject_id, unique_id])

        anat_rsc = create_anat_datasource('anat_gather_%s' % node_suffix)

        anat_rsc.inputs.inputnode.set(
            subject=subject_id,
            anat=session['anat'],
            creds_path=input_creds_path,
            dl_dir=config.pipeline_setup['working_directory']['path'],
            img_type='anat'
        )

        strat.update_resource_pool({
            'anatomical': (anat_rsc, 'outputspec.anat')
        })

        strat.update_resource_pool({
            'template_cmass': (template_center_of_mass, 'cm')
        })

        # Here we have the same strategies for the skull stripping as in prep_workflow
        if 'brain_mask' in session.keys() and session['brain_mask'] and \
                        session['brain_mask'].lower() != 'none':

            brain_rsc = create_anat_datasource(
                'brain_gather_%s' % unique_id)
            brain_rsc.inputs.inputnode.set(
                subject=subject_id,
                anat=session['brain_mask'],
                creds_path=input_creds_path,
                dl_dir=config.pipeline_setup['working_directory']['path'],
                img_type='anat'
            )

            skullstrip_method = 'mask'
            preproc_wf_name = 'anat_preproc_mask_%s' % node_suffix

            strat.append_name(brain_rsc.name)
            strat.update_resource_pool({
                'anatomical_brain_mask': (brain_rsc, 'outputspec.anat')
            })

            anat_preproc = create_anat_preproc(
                method=skullstrip_method,
                config=config,
                wf_name=preproc_wf_name)

            workflow.connect(brain_rsc, 'outputspec.brain_mask',
                             anat_preproc, 'inputspec.brain_mask')

            new_strat, strat_nodes_list_list = connect_anat_preproc_inputs(
                strat, anat_preproc, skullstrip_method + "_skullstrip",
                strat_nodes_list_list, workflow)
            strat_list.append(new_strat)

        elif already_skullstripped:
            skullstrip_method = None
            preproc_wf_name = 'anat_preproc_already_%s' % node_suffix
            anat_preproc = create_anat_preproc(
                method=skullstrip_method,
                already_skullstripped=True,
                config=config,
                wf_name=preproc_wf_name
            )
            new_strat, strat_nodes_list_list = connect_anat_preproc_inputs(
                strat, anat_preproc, 'already_skullstripped',
                strat_nodes_list_list, workflow)
            strat_list.append(new_strat)

        else:
            # TODO add other SS methods
            if "3dSkullStrip" in \
                    config.anatomical_preproc['brain_extraction'][
                        'extraction']['using']:
                skullstrip_method = 'afni'
                preproc_wf_name = 'anat_preproc_afni_%s' % node_suffix

                anat_preproc = create_anat_preproc(
                    method=skullstrip_method,
                    config=config,
                    wf_name=preproc_wf_name)
                '''
                anat_preproc.inputs.AFNI_options.set(
                    mask_vol=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['mask_vol'],
                    shrink_factor=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['shrink_factor'],
                    var_shrink_fac=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['var_shrink_fac'],
                    shrink_fac_bot_lim=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['shrink_factor_bot_lim'],
                    avoid_vent=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['avoid_vent'],
                    niter=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['n_iterations'],
                    pushout=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['pushout'],
                    touchup=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['touchup'],
                    fill_hole=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['fill_hole'],
                    NN_smooth=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['NN_smooth'],
                    smooth_final=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['smooth_final'],
                    avoid_eyes=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['avoid_eyes'],
                    use_edge=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['use_edge'],
                    exp_frac=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['exp_frac'],
                    push_to_edge=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['push_to_edge'],
                    use_skull=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['use_skull'],
                    perc_int=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['perc_int'],
                    max_inter_iter=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['max_inter_iter'],
                    fac=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['fac'],
                    blur_fwhm=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['blur_fwhm'],
                    monkey=config.anatomical_preproc['brain_extraction']['extraction']['AFNI-3dSkullStrip']['monkey'],
                )
                '''
                new_strat, strat_nodes_list_list = connect_anat_preproc_inputs(
                    strat, anat_preproc, skullstrip_method + "_skullstrip",
                    strat_nodes_list_list, workflow)
                strat_list.append(new_strat)

            if "BET" in config.anatomical_preproc['brain_extraction'][
                'extraction']['using']:
                skullstrip_method = 'fsl'
                preproc_wf_name = 'anat_preproc_fsl_%s' % node_suffix

                anat_preproc = create_anat_preproc(
                    method=skullstrip_method,
                    config=config,
                    wf_name=preproc_wf_name)
                '''
                anat_preproc.inputs.BET_options.set(
                    frac=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['frac'],
                    mask_boolean=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['mask_boolean'],
                    mesh_boolean=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['mesh_boolean'],
                    outline=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['outline'],
                    padding=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['padding'],
                    radius=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['radius'],
                    reduce_bias=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['reduce_bias'],
                    remove_eyes=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['remove_eyes'],
                    robust=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['robust'],
                    skull=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['skull'],
                    surfaces=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['surfaces'],
                    threshold=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['threshold'],
                    vertical_gradient=config.anatomical_preproc['brain_extraction']['extraction']['FSL-BET']['vertical_gradient'],
                )
                '''
                new_strat, strat_nodes_list_list = connect_anat_preproc_inputs(
                    strat, anat_preproc, skullstrip_method + "_skullstrip",
                    strat_nodes_list_list, workflow)
                strat_list.append(new_strat)

            if not any(o in config.anatomical_preproc['brain_extraction'][
                'extraction']['using'] for o in
                       ["3dSkullStrip", "BET"]):
                err = '\n\n[!] C-PAC says: Your skull-stripping ' \
                      'method options setting does not include either' \
                      ' \'3dSkullStrip\' or \'BET\'.\n\n Options you ' \
                      'provided:\nbrain_extraction: \n' \
                      'extraction: \n' \
                      'using:[{0}]\n\n'.format(
                    str(config.anatomical_preproc['brain_extraction'][
                            'extraction']['using']))
                raise Exception(err)

    # Here we have all the anat_preproc set up for every session of the subject

    strat_init = Strategy()

    templates_for_resampling = [
        (config.anatomical_preproc['registration_workflow'][
             'resolution_for_anat'],
         config.anatomical_preproc['registration_workflow'][
             'template_brain_only_for_anat'], 'template_brain_for_anat',
         'resolution_for_anat'),
        (config.anatomical_preproc['registration_workflow'][
             'resolution_for_anat'],
         config.anatomical_preproc['registration_workflow'][
             'template_skull_for_anat'], 'template_skull_for_anat',
         'resolution_for_anat'),
        (config.anatomical_preproc['registration_workflow'][
             'resolution_for_anat'],
         config.voxel_mirrored_homotopic_connectivity[
             'symmetric_registration']['template_symmetric_brain_only'],
         'template_symmetric_brain', 'resolution_for_anat'),
        (config.anatomical_preproc['registration_workflow'][
             'resolution_for_anat'],
         config.voxel_mirrored_homotopic_connectivity[
             'symmetric_registration']['template_symmetric_skull'],
         'template_symmetric_skull', 'resolution_for_anat'),
        (config.anatomical_preproc['registration_workflow'][
             'resolution_for_anat'],
         config.voxel_mirrored_homotopic_connectivity[
             'symmetric_registration']['dilated_symmetric_brain_mask'],
         'template_dilated_symmetric_brain_mask', 'resolution_for_anat'),
        (config.anatomical_preproc['registration_workflow'][
             'resolution_for_anat'],
         config.anatomical_preproc['registration_workflow']['registration'][
             'FSL-FNIRT']['ref_mask'], 'template_ref_mask',
         'resolution_for_anat'),
        (config.functional_registration['2-func_registration_to_template'][
             'output_resolution']['func_preproc_outputs'],
         config.functional_registration['2-func_registration_to_template'][
             'target_template']['T1_template']['template_brain'],
         'template_brain_for_func_preproc', 'resolution_for_func_preproc'),
        (config.functional_registration['2-func_registration_to_template'][
             'output_resolution']['func_preproc_outputs'],
         config.functional_registration['2-func_registration_to_template'][
             'target_template']['T1_template']['template_skull'],
         'template_skull_for_func_preproc', 'resolution_for_func_preproc'),
        (config.functional_registration['2-func_registration_to_template'][
             'output_resolution']['func_derivative_outputs'],
         config.functional_registration['2-func_registration_to_template'][
             'target_template']['T1_template']['template_brain'],
         'template_brain_for_func_derivative', 'resolution_for_func_preproc'),
        (config.functional_registration['2-func_registration_to_template'][
             'output_resolution']['func_derivative_outputs'],
         config.functional_registration['2-func_registration_to_template'][
             'target_template']['T1_template']['template_skull'],
         'template_skull_for_func_derivative', 'resolution_for_func_preproc'),
    ]

    # update resampled template to resource pool
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

    # loop over the different skull stripping strategies
    for strat_name, strat_nodes_list in strat_nodes_list_list.items():

        node_suffix = '_'.join([strat_name, subject_id])

        # Merge node to feed the anat_preproc outputs to the longitudinal template generation
        brain_merge_node = pe.Node(
            interface=Merge(len(strat_nodes_list)),
            name="anat_longitudinal_brain_merge_" + node_suffix)

        skull_merge_node = pe.Node(
            interface=Merge(len(strat_nodes_list)),
            name="anat_longitudinal_skull_merge_" + node_suffix)

        # This node will generate the longitudinal template (the functions are in longitudinal_preproc)
        # Later other algorithms could be added to calculate it, like the multivariate template from ANTS
        # It would just require to change it here.
        template_node = subject_specific_template(
            workflow_name='subject_specific_anat_template_' + node_suffix
        )

        unique_id_list = [i.get_name()[0].split('_')[-1] for i in
                          strat_nodes_list]

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
            unique_id_list=unique_id_list
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

        # Update resource pool
        # longitudinal template
        rsc_key = 'anatomical_longitudinal_template_'
        ds_template = create_datasink(rsc_key + node_suffix, config,
                                      subject_id,
                                      strat_name='longitudinal_' + strat_name)
        workflow.connect(template_node, 'brain_template',
                         ds_template, rsc_key)

        # T1 to longitudinal template warp
        rsc_key = 'anatomical_to_longitudinal_template_warp_'
        ds_warp_list = create_datasink(rsc_key + node_suffix, config,
                                       subject_id,
                                       strat_name='longitudinal_' + strat_name,
                                       map_node_iterfield=[
                                           'anatomical_to_longitudinal_template_warp'])
        workflow.connect(template_node, "warp_list",
                         ds_warp_list,
                         'anatomical_to_longitudinal_template_warp')

        # T1 in longitudinal template space
        rsc_key = 'anatomical_to_longitudinal_template_'
        t1_list = create_datasink(rsc_key + node_suffix, config, subject_id,
                                  strat_name='longitudinal_' + strat_name,
                                  map_node_iterfield=[
                                      'anatomical_to_longitudinal_template'])
        workflow.connect(template_node, "output_brain_list",
                         t1_list, 'anatomical_to_longitudinal_template')

        # longitudinal to standard registration items
        for num_strat, strat in enumerate(reg_strat_list):
            for rsc_key in strat.resource_pool.keys():
                rsc_nodes_suffix = '_'.join(
                    ['_longitudinal_to_standard', strat_name, str(num_strat)])
                if rsc_key in Outputs.any:
                    node, rsc_name = strat[rsc_key]
                    ds = create_datasink(rsc_key + rsc_nodes_suffix, config,
                                         subject_id,
                                         strat_name='longitudinal_' + strat_name)
                    workflow.connect(node, rsc_name, ds, rsc_key)

        # individual minimal preprocessing items
        for i in range(len(strat_nodes_list)):
            rsc_nodes_suffix = "_%s_%d" % (node_suffix, i)
            for rsc_key in strat_nodes_list[i].resource_pool.keys():
                if rsc_key in Outputs.any:
                    node, rsc_name = strat_nodes_list[i][rsc_key]
                    ds = create_datasink(rsc_key + rsc_nodes_suffix, config,
                                         subject_id,
                                         session_id_list[i],
                                         'longitudinal_' + strat_name)
                    workflow.connect(node, rsc_name, ds, rsc_key)

            rsc_key = 'anatomical_brain'
            anat_preproc_node, rsc_name = strat_nodes_list[i][rsc_key]
            workflow.connect(anat_preproc_node,
                             rsc_name, brain_merge_node,
                             'in{}'.format(
                                 i + 1))  # the in{}.format take i+1 because the Merge nodes inputs starts at 1

            rsc_key = 'anatomical_skull_leaf'
            anat_preproc_node, rsc_name = strat_nodes_list[i][rsc_key]
            workflow.connect(anat_preproc_node,
                             rsc_name, skull_merge_node,
                             'in{}'.format(i + 1))

    workflow.run()

    return reg_strat_list  # strat_nodes_list_list # for func wf?


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
