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
# pylint: disable=too-many-lines,ungrouped-imports,wrong-import-order
from typing import Optional
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nodeblock import nodeblock
from nipype.interfaces import afni, ants, c3, fsl, utility as util
from nipype.interfaces.afni import utils as afni_utils

from CPAC.anat_preproc.lesion_preproc import create_lesion_preproc
from CPAC.func_preproc.utils import chunk_ts, split_ts_chunks
from CPAC.registration.utils import seperate_warps_list, \
                                    check_transforms, \
                                    generate_inverse_transform_flags, \
                                    single_ants_xfm_to_list, \
                                    interpolation_string, \
                                    change_itk_transform_type, \
                                    hardcoded_reg, \
                                    one_d_to_mat, \
                                    run_c3d, \
                                    run_c4d
from CPAC.utils.interfaces.fsl import Merge as fslMerge
from CPAC.utils.typing import LIST_OR_STR, TUPLE
from CPAC.utils.utils import check_prov_for_motion_tool, check_prov_for_regtool


def apply_transform(wf_name, reg_tool, time_series=False, multi_input=False,
                    num_cpus=1, num_ants_cores=1):

    if not reg_tool:
        raise Exception("\n[!] Developer info: the 'reg_tool' parameter sent "
                        f"to the 'apply_transform' node for '{wf_name}' is "
                        f"empty.\n")

    wf = pe.Workflow(name=wf_name)

    inputNode = pe.Node(
        util.IdentityInterface(fields=['input_image',
                                       'reference',
                                       'transform',
                                       'interpolation']),
        name='inputspec')

    outputNode = pe.Node(
        util.IdentityInterface(fields=['output_image']),
        name='outputspec')

    if int(num_cpus) > 1 and time_series:
        # parallelize time series warp application
        # we need the node to be a MapNode to feed in the list of functional
        # time series chunks
        multi_input = True

    if reg_tool == 'ants':

        if multi_input:
            apply_warp = pe.MapNode(interface=ants.ApplyTransforms(),
                                    name=f'apply_warp_{wf_name}',
                                    iterfield=['input_image'],
                                    mem_gb=0.7,
                                    mem_x=(1708448960473801 /
                                           151115727451828646838272,
                                           'input_image'))
        else:
            apply_warp = pe.Node(interface=ants.ApplyTransforms(),
                                 name=f'apply_warp_{wf_name}',
                                 mem_gb=0.7,
                                 mem_x=(1708448960473801 /
                                        151115727451828646838272,
                                        'input_image'))

        apply_warp.inputs.dimension = 3
        apply_warp.interface.num_threads = int(num_ants_cores)

        if time_series:
            apply_warp.inputs.input_image_type = 3

        wf.connect(inputNode, 'reference', apply_warp, 'reference_image')

        interp_string = pe.Node(util.Function(input_names=['interpolation',
                                                           'reg_tool'],
                                              output_names=['interpolation'],
                                              function=interpolation_string),
                                name=f'interp_string',
                                mem_gb=2.5)
        interp_string.inputs.reg_tool = reg_tool

        wf.connect(inputNode, 'interpolation', interp_string, 'interpolation')
        wf.connect(interp_string, 'interpolation',
                   apply_warp, 'interpolation')

        ants_xfm_list = \
            pe.Node(util.Function(input_names=['transform'],
                                  output_names=['transform_list'],
                                  function=single_ants_xfm_to_list),
                    name=f'single_ants_xfm_to_list',
                    mem_gb=2.5)

        wf.connect(inputNode, 'transform', ants_xfm_list, 'transform')
        wf.connect(ants_xfm_list, 'transform_list', apply_warp, 'transforms')

        # parallelize the apply warp, if multiple CPUs, and it's a time
        # series!
        if int(num_cpus) > 1 and time_series:

            chunk_imports = ['import nibabel as nb']
            chunk = pe.Node(util.Function(input_names=['func_file',
                                                       'n_chunks',
                                                       'chunk_size'],
                                     output_names=['TR_ranges'],
                                     function=chunk_ts,
                                     imports=chunk_imports),
                            name=f'chunk_{wf_name}',
                            mem_gb=2.5)

            #chunk.inputs.n_chunks = int(num_cpus)

            # 10-TR sized chunks
            chunk.inputs.chunk_size = 10

            wf.connect(inputNode, 'input_image', chunk, 'func_file')

            split_imports = ['import os', 'import subprocess']
            split = pe.Node(util.Function(input_names=['func_file',
                                                       'tr_ranges'],
                                     output_names=['split_funcs'],
                                     function=split_ts_chunks,
                                     imports=split_imports),
                            name=f'split_{wf_name}',
                            mem_gb=2.5)

            wf.connect(inputNode, 'input_image', split, 'func_file')
            wf.connect(chunk, 'TR_ranges', split, 'tr_ranges')

            wf.connect(split, 'split_funcs', apply_warp, 'input_image')

            func_concat = pe.Node(interface=afni_utils.TCat(),
                                  name=f'func_concat_{wf_name}',
                                  mem_gb=2.5)
            func_concat.inputs.outputtype = 'NIFTI_GZ'

            wf.connect(apply_warp, 'output_image', func_concat, 'in_files')

            wf.connect(func_concat, 'out_file', outputNode, 'output_image')

        else:
            wf.connect(inputNode, 'input_image', apply_warp, 'input_image')
            wf.connect(apply_warp, 'output_image', outputNode, 'output_image')

    elif reg_tool == 'fsl':

        if multi_input:
            apply_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                                    name=f'fsl_apply_warp',
                                    iterfield=['in_file'],
                                    mem_gb=2.5)
        else:
            apply_warp = pe.Node(interface=fsl.ApplyWarp(),
                                 name='fsl_apply_warp',
                                 mem_gb=2.5)

        interp_string = pe.Node(util.Function(input_names=['interpolation',
                                                           'reg_tool'],
                                              output_names=['interpolation'],
                                              function=interpolation_string),
                                name=f'interp_string',
                                mem_gb=2.5)
        interp_string.inputs.reg_tool = reg_tool

        wf.connect(inputNode, 'interpolation', interp_string, 'interpolation')
        wf.connect(interp_string, 'interpolation', apply_warp, 'interp')

        # mni to t1
        wf.connect(inputNode, 'reference', apply_warp, 'ref_file')

        # NOTE: C-PAC now converts all FSL xfm's to .nii, so even if the
        #       inputNode 'transform' is a linear xfm, it's a .nii and must
        #       go in as a warpfield file
        wf.connect(inputNode, 'transform', apply_warp, 'field_file')

        # parallelize the apply warp, if multiple CPUs, and it's a time
        # series!
        if int(num_cpus) > 1 and time_series:

            chunk_imports = ['import nibabel as nb']
            chunk = pe.Node(util.Function(input_names=['func_file',
                                                       'n_chunks',
                                                       'chunk_size'],
                                     output_names=['TR_ranges'],
                                     function=chunk_ts,
                                     imports=chunk_imports),
                            name=f'chunk_{wf_name}',
                            mem_gb=2.5)

            #chunk.inputs.n_chunks = int(num_cpus)

            # 10-TR sized chunks
            chunk.inputs.chunk_size = 10

            wf.connect(inputNode, 'input_image', chunk, 'func_file')

            split_imports = ['import os', 'import subprocess']
            split = pe.Node(util.Function(input_names=['func_file',
                                                       'tr_ranges'],
                                     output_names=['split_funcs'],
                                     function=split_ts_chunks,
                                     imports=split_imports),
                            name=f'split_{wf_name}',
                            mem_gb=2.5)

            wf.connect(inputNode, 'input_image', split, 'func_file')
            wf.connect(chunk, 'TR_ranges', split, 'tr_ranges')

            wf.connect(split, 'split_funcs', apply_warp, 'in_file')

            func_concat = pe.Node(interface=afni_utils.TCat(),
                                  name=f'func_concat{wf_name}')
            func_concat.inputs.outputtype = 'NIFTI_GZ'

            wf.connect(apply_warp, 'out_file', func_concat, 'in_files')

            wf.connect(func_concat, 'out_file', outputNode, 'output_image')

        else:
            wf.connect(inputNode, 'input_image', apply_warp, 'in_file')
            wf.connect(apply_warp, 'out_file', outputNode, 'output_image')

    return wf


def transform_derivative(wf_name, label, reg_tool, num_cpus, num_ants_cores,
                         ants_interp=None, fsl_interp=None, opt=None):
    '''Transform output derivatives to template space.

    This function is designed for use with the NodeBlock connection engine.
    '''

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['in_file',
                                                       'reference',
                                                       'transform']),
                        name='inputspec')

    multi_input = False
    if 'statmap' in label:
        multi_input = True

    stack = False
    if 'correlations' in label:
        stack = True

    apply_xfm = apply_transform(f'warp_{label}_to_template', reg_tool,
                                time_series=stack,
                                multi_input=multi_input,
                                num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)

    if reg_tool == 'ants':
        apply_xfm.inputs.inputspec.interpolation = ants_interp
    elif reg_tool == 'fsl':
        apply_xfm.inputs.inputspec.interpolation = fsl_interp

    wf.connect(inputnode, 'in_file', apply_xfm, 'inputspec.input_image')
    wf.connect(inputnode, 'reference', apply_xfm, 'inputspec.reference')
    wf.connect(inputnode, 'transform', apply_xfm, 'inputspec.transform')

    outputnode = pe.Node(util.IdentityInterface(fields=['out_file']),
                         name='outputspec')

    wf.connect(apply_xfm, 'outputspec.output_image', outputnode, 'out_file')

    return wf


def convert_pedir(pedir, convert='xyz_to_int'):
    '''FSL Flirt requires pedir input encoded as an int'''
    if convert == 'xyz_to_int':
        conv_dct = {'x': 1, 'y': 2, 'z': 3, 'x-': -1, 'y-': -2, 'z-': -3,
                    'i': 1, 'j': 2, 'k': 3, 'i-': -1, 'j-': -2, 'k-': -3,
                    '-x': -1, '-i': -1, '-y': -2,
                    '-j': -2, '-z': -3, '-k': -3}
    elif convert == 'ijk_to_xyz':
        conv_dct = {'i': 'x', 'j': 'y', 'k': 'z',
                    'i-': 'x-', 'j-': 'y-', 'k-': 'z-'}

    if isinstance(pedir, bytes):
        pedir = pedir.decode()
    if not isinstance(pedir, str):
        raise Exception("\n\nPhase-encoding direction must be a "
                        "string value.\n\nValue: {0}"
                        "\n\n".format(pedir))
    if pedir not in conv_dct.keys():
        raise Exception("\n\nInvalid phase-encoding direction "
                        "entered: {0}\n\n".format(pedir))
    pedir = conv_dct[pedir]
    return pedir


def create_fsl_flirt_linear_reg(name='fsl_flirt_linear_reg'):

    linear_register = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['input_brain',
                                                       'reference_brain',
                                                       'interp',
                                                       'ref_mask']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['output_brain',
                                                        'linear_xfm',
                                                        'invlinear_xfm']),
                         name='outputspec')

    linear_reg = pe.Node(interface=fsl.FLIRT(), name='linear_reg_0')
    linear_reg.inputs.cost = 'corratio'

    inv_flirt_xfm = pe.Node(interface=fsl.utils.ConvertXFM(),
                            name='inv_linear_reg0_xfm')
    inv_flirt_xfm.inputs.invert_xfm = True

    linear_register.connect(inputspec, 'input_brain',
                            linear_reg, 'in_file')

    linear_register.connect(inputspec, 'reference_brain',
                            linear_reg, 'reference')

    linear_register.connect(inputspec, 'interp',
                            linear_reg, 'interp')

    linear_register.connect(linear_reg, 'out_file',
                            outputspec, 'output_brain')

    linear_register.connect(linear_reg, 'out_matrix_file',
                            inv_flirt_xfm, 'in_file')

    linear_register.connect(inv_flirt_xfm, 'out_file',
                            outputspec, 'invlinear_xfm')

    linear_register.connect(linear_reg, 'out_matrix_file',
                            outputspec, 'linear_xfm')

    return linear_register


def create_fsl_fnirt_nonlinear_reg(name='fsl_fnirt_nonlinear_reg'):
    """
    Performs non-linear registration of an input file to a reference file
    using FSL FNIRT.

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    nonlinear_register : nipype.pipeline.engine.Workflow

    Notes
    -----

    Workflow Inputs::

        inputspec.input_skull : string (nifti file)
            File of input brain with skull
        inputspec.reference_skull : string (nifti file)
            Target brain with skull to normalize to
        inputspec.fnirt_config : string (fsl fnirt config file)
            Configuration file containing parameters that can be
            specified in fnirt
    Workflow Outputs::

        outputspec.output_brain : string (nifti file)
            Normalizion of input brain file
        outputspec.nonlinear_xfm : string
            Nonlinear field coefficients file of nonlinear transformation

    Registration Procedure:

    1. Perform a nonlinear registration on an input file to the
       reference file utilizing affine transformation from the previous
       step as a starting point.
    2. Invert the affine transformation to provide the user a
       transformation (affine only) from the space of the reference
       file to the input file.

    Workflow Graph:

    .. image:: ../images/nonlinear_register.dot.png
        :width: 500

    Detailed Workflow Graph:
    .. image:: ../images/nonlinear_register_detailed.dot.png
        :width: 500
    """

    nonlinear_register = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['input_brain',
                                                       'input_skull',
                                                       'reference_brain',
                                                       'reference_skull',
                                                       'interp',
                                                       'ref_mask',
                                                       'linear_aff',
                                                       'fnirt_config']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['output_brain',
                                                        'nonlinear_xfm']),
                         name='outputspec')

    nonlinear_reg = pe.Node(interface=fsl.FNIRT(),
                            name='nonlinear_reg_1')

    nonlinear_reg.inputs.fieldcoeff_file = True
    nonlinear_reg.inputs.jacobian_file = True

    brain_warp = pe.Node(interface=fsl.ApplyWarp(),
                         name='brain_warp')

    nonlinear_register.connect(inputspec, 'input_skull',
                               nonlinear_reg, 'in_file')

    nonlinear_register.connect(inputspec, 'reference_skull',
                               nonlinear_reg, 'ref_file')

    nonlinear_register.connect(inputspec, 'interp',
                               brain_warp, 'interp')

    nonlinear_register.connect(inputspec, 'ref_mask',
                               nonlinear_reg, 'refmask_file')

    # FNIRT parameters are specified by FSL config file
    # ${FSLDIR}/etc/flirtsch/TI_2_MNI152_2mm.cnf (or user-specified)
    nonlinear_register.connect(inputspec, 'fnirt_config',
                               nonlinear_reg, 'config_file')

    nonlinear_register.connect(inputspec, 'linear_aff',
                               nonlinear_reg, 'affine_file')

    nonlinear_register.connect(nonlinear_reg, 'fieldcoeff_file',
                               outputspec, 'nonlinear_xfm')

    nonlinear_register.connect(inputspec, 'input_brain',
                               brain_warp, 'in_file')

    nonlinear_register.connect(nonlinear_reg, 'fieldcoeff_file',
                               brain_warp, 'field_file')

    nonlinear_register.connect(inputspec, 'reference_brain',
                               brain_warp, 'ref_file')

    nonlinear_register.connect(brain_warp, 'out_file',
                               outputspec, 'output_brain')

    return nonlinear_register


def create_fsl_fnirt_nonlinear_reg_nhp(name='fsl_fnirt_nonlinear_reg_nhp'):
    """
    Performs non-linear registration of an input file to a reference file
    using FSL FNIRT.

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    nonlinear_register : nipype.pipeline.engine.Workflow

    Notes
    -----

    Workflow Inputs::

        inputspec.input_skull : string (nifti file)
            File of input brain with skull
        inputspec.reference_skull : string (nifti file)
            Target brain with skull to normalize to
        inputspec.fnirt_config : string (fsl fnirt config file)
            Configuration file containing parameters that can be
            specified in fnirt
    Workflow Outputs::

        outputspec.output_brain : string (nifti file)
            Normalizion of input brain file
        outputspec.nonlinear_xfm : string
            Nonlinear field coefficients file of nonlinear transformation
        outputspec.nonlinear_warp : string
            Nonlinear output file with warp field

    Registration Procedure:

    1. Perform a nonlinear registration on an input file to the
       reference file utilizing affine transformation from the previous
       step as a starting point.
    2. Invert the affine transformation to provide the user a
       transformation (affine only) from the space of the reference
       file to the input file.

    Workflow Graph:

    .. image:: ../images/nonlinear_register.dot.png
        :width: 500

    Detailed Workflow Graph:
    .. image:: ../images/nonlinear_register_detailed.dot.png
        :width: 500
    """

    nonlinear_register = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['input_brain',
                                                       'input_skull',
                                                       'reference_brain',
                                                       'reference_skull',
                                                       'interp',
                                                       'ref_mask',
                                                       'linear_aff',
                                                       'fnirt_config']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['output_brain',
                                                        'output_head',
                                                        'output_mask',
                                                        'output_biasfield',
                                                        'nonlinear_xfm',
                                                        'nonlinear_warp']),
                         name='outputspec')

    nonlinear_reg = pe.Node(interface=fsl.FNIRT(),
                            name='nonlinear_reg_1')

    nonlinear_reg.inputs.fieldcoeff_file = True
    nonlinear_reg.inputs.jacobian_file = True
    nonlinear_reg.inputs.field_file = True

    nonlinear_register.connect(inputspec, 'input_skull',
                               nonlinear_reg, 'in_file')

    nonlinear_register.connect(inputspec, 'reference_skull',
                               nonlinear_reg, 'ref_file')

    nonlinear_register.connect(inputspec, 'ref_mask',
                               nonlinear_reg, 'refmask_file')

    nonlinear_register.connect(inputspec, 'fnirt_config',
                               nonlinear_reg, 'config_file')

    nonlinear_register.connect(inputspec, 'linear_aff',
                               nonlinear_reg, 'affine_file')

    brain_warp = pe.Node(interface=fsl.ApplyWarp(),
                         name='brain_warp')
    brain_warp.inputs.interp = 'nn'
    brain_warp.inputs.relwarp = True

    nonlinear_register.connect(inputspec, 'input_brain',
                               brain_warp, 'in_file')

    nonlinear_register.connect(nonlinear_reg, 'field_file',
                               brain_warp, 'field_file')

    nonlinear_register.connect(inputspec, 'reference_skull',
                               brain_warp, 'ref_file')

    head_warp = pe.Node(interface=fsl.ApplyWarp(),
                         name='head_warp')
    head_warp.inputs.interp = 'spline'
    head_warp.inputs.relwarp = True

    nonlinear_register.connect(inputspec, 'input_brain',
                               head_warp, 'in_file')

    nonlinear_register.connect(nonlinear_reg, 'field_file',
                               head_warp, 'field_file')

    nonlinear_register.connect(inputspec, 'reference_skull',
                               head_warp, 'ref_file')

    mask_warp = pe.Node(interface=fsl.ApplyWarp(),
                         name='mask_warp')
    mask_warp.inputs.interp = 'nn'
    mask_warp.inputs.relwarp = True

    nonlinear_register.connect(inputspec, 'input_brain',
                               mask_warp, 'in_file')

    nonlinear_register.connect(nonlinear_reg, 'field_file',
                               mask_warp, 'field_file')

    nonlinear_register.connect(inputspec, 'reference_skull',
                               mask_warp, 'ref_file')

    biasfield_warp = pe.Node(interface=fsl.ApplyWarp(),
                         name='biasfield_warp')
    biasfield_warp.inputs.interp = 'spline'
    biasfield_warp.inputs.relwarp = True

    nonlinear_register.connect(inputspec, 'input_brain',
                               biasfield_warp, 'in_file')

    nonlinear_register.connect(nonlinear_reg, 'field_file',
                               biasfield_warp, 'field_file')

    nonlinear_register.connect(inputspec, 'reference_skull',
                               biasfield_warp, 'ref_file')

    nonlinear_register.connect(nonlinear_reg, 'fieldcoeff_file',
                               outputspec, 'nonlinear_xfm')

    nonlinear_register.connect(nonlinear_reg, 'field_file',
                               outputspec, 'nonlinear_warp')

    nonlinear_register.connect(brain_warp, 'out_file',
                               outputspec, 'output_brain')

    nonlinear_register.connect(head_warp, 'out_file',
                               outputspec, 'output_head')

    nonlinear_register.connect(mask_warp, 'out_file',
                               outputspec, 'output_mask')

    nonlinear_register.connect(biasfield_warp, 'out_file',
                               outputspec, 'output_biasfield')

    return nonlinear_register


def create_register_func_to_anat(config, phase_diff_distcor=False,
                                 name='register_func_to_anat'):

    """
    Registers a functional scan in native space to anatomical space using a
    linear transform and does not include bbregister.

    Parameters
    ----------
    config : configuration, mandatory
        Pipeline configuration.
    fieldmap_distortion : bool, optional
        If field map-based distortion correction is being run, FLIRT should
        take in the appropriate field map-related inputs.
    name : string, optional
        Name of the workflow.

    Returns
    -------
    create_register_func_to_anat : nipype.pipeline.engine.Workflow

    Notes
    -----

    Workflow Inputs::

        inputspec.func : string (nifti file)
            Input functional scan to be registered to anatomical space
        inputspec.anat : string (nifti file)
            Corresponding anatomical scan of subject
        inputspec.interp : string
            Type of interpolation to use
            ('trilinear' or 'nearestneighbour' or 'sinc')

    Workflow Outputs::

        outputspec.func_to_anat_linear_xfm_nobbreg : string (mat file)
            Affine transformation from functional to anatomical native space
        outputspec.anat_func_nobbreg : string (nifti file)
            Functional scan registered to anatomical space
    """
    register_func_to_anat = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['func',
                                                       'anat',
                                                       'dof',
                                                       'interp',
                                                       'fieldmap',
                                                       'fieldmapmask']),
                        name='inputspec')

    inputNode_echospacing = pe.Node(
        util.IdentityInterface(fields=['echospacing']),
        name='echospacing_input')

    inputNode_pedir = pe.Node(util.IdentityInterface(fields=['pedir']),
                              name='pedir_input')

    outputspec = pe.Node(util.IdentityInterface(
        fields=['func_to_anat_linear_xfm_nobbreg', 'anat_func_nobbreg']),
        name='outputspec')

    linear_reg = pe.Node(interface=fsl.FLIRT(),
                         name='linear_func_to_anat')

    linear_reg.inputs.interp = config.registration_workflows['functional_registration']['coregistration']['interpolation']
    linear_reg.inputs.cost = config.registration_workflows['functional_registration']['coregistration']['cost']
    linear_reg.inputs.dof = config.registration_workflows['functional_registration']['coregistration']['dof']
    if config.registration_workflows['functional_registration']['coregistration']['arguments'] is not None:
        linear_reg.inputs.args = config.registration_workflows['functional_registration']['coregistration']['arguments']

    if phase_diff_distcor:
        conv_pedir = \
            pe.Node(interface=util.Function(input_names=['pedir',
                                                         'convert'],
                                            output_names=['pedir'],
                                            function=convert_pedir),
                    name='coreg_convert_pedir')
        conv_pedir.inputs.convert = 'xyz_to_int'

        register_func_to_anat.connect(inputNode_pedir, 'pedir',
                                      conv_pedir, 'pedir')
        register_func_to_anat.connect(conv_pedir, 'pedir', 
                                      linear_reg, 'pedir')
        register_func_to_anat.connect(inputspec, 'fieldmap',
                                      linear_reg, 'fieldmap')
        register_func_to_anat.connect(inputspec, 'fieldmapmask',
                                      linear_reg, 'fieldmapmask')
        register_func_to_anat.connect(inputNode_echospacing, 'echospacing',
                                      linear_reg, 'echospacing')

    register_func_to_anat.connect(inputspec, 'func', linear_reg, 'in_file')

    register_func_to_anat.connect(inputspec, 'anat', linear_reg, 'reference')

    register_func_to_anat.connect(inputspec, 'dof', linear_reg, 'dof')

    register_func_to_anat.connect(inputspec, 'interp', linear_reg, 'interp')

    register_func_to_anat.connect(linear_reg, 'out_matrix_file',
                                  outputspec,
                                  'func_to_anat_linear_xfm_nobbreg')

    register_func_to_anat.connect(linear_reg, 'out_file',
                                  outputspec, 'anat_func_nobbreg')

    return register_func_to_anat


def create_register_func_to_anat_use_T2(config, name='register_func_to_anat_use_T2'):
    # for monkey data
    # ref: https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh#L287-L295
    # https://github.com/HechengJin0/dcan-macaque-pipeline/blob/master/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh#L524-L535

    """
    Registers a functional scan in native space to anatomical space using a
    linear transform and does not include bbregister, use T1 and T2 image.

    Parameters
    ----------
    config : configuration, mandatory
        Pipeline configuration.
    name : string, optional
        Name of the workflow.

    Returns
    -------
    create_register_func_to_anat_use_T2 : nipype.pipeline.engine.Workflow

    Notes
    -----

    Workflow Inputs::

        inputspec.func : string (nifti file)
            Input functional scan to be registered to anatomical space
        inputspec.anat : string (nifti file)
            Corresponding anatomical scan of subject

    Workflow Outputs::

        outputspec.func_to_anat_linear_xfm_nobbreg : string (mat file)
            Affine transformation from functional to anatomical native space
        outputspec.anat_func_nobbreg : string (nifti file)
            Functional scan registered to anatomical space
    """


    register_func_to_anat_use_T2 = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['func',
                                                       'T1_brain',
                                                       'T2_head',
                                                       'T2_brain']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['func_to_anat_linear_xfm_nobbreg',
                                                        'func_to_anat_linear_warp_nobbreg',
                                                        'anat_func_nobbreg']),
                                                        name='outputspec')

    # ${FSLDIR}/bin/flirt -interp spline -dof 6 -in ${fMRIFolder}/${ScoutName}_gdc -ref ${T1wFolder}/${T2wRestoreImage} -omat "$fMRIFolder"/Scout2T2w.mat -out ${fMRIFolder}/Scout2T2w.nii.gz -searchrx -30 30 -searchry -30 30 -searchrz -30 30 -cost mutualinfo
    linear_reg_func_to_t2 = pe.Node(interface=fsl.FLIRT(),
                         name='linear_reg_func_to_t2')
    linear_reg_func_to_t2.inputs.interp = 'spline'
    linear_reg_func_to_t2.inputs.cost = 'mutualinfo'
    linear_reg_func_to_t2.inputs.dof = 6
    linear_reg_func_to_t2.inputs.searchr_x = [30, 30]
    linear_reg_func_to_t2.inputs.searchr_y = [30, 30]
    linear_reg_func_to_t2.inputs.searchr_z = [30, 30]

    register_func_to_anat_use_T2.connect(inputspec, 'func', linear_reg_func_to_t2, 'in_file')

    register_func_to_anat_use_T2.connect(inputspec, 'T2_head', linear_reg_func_to_t2, 'reference')

    # ${FSLDIR}/bin/convert_xfm -omat "$fMRIFolder"/T2w2Scout.mat -inverse "$fMRIFolder"/Scout2T2w.mat
    invt = pe.Node(interface=fsl.ConvertXFM(), name='convert_xfm')
    invt.inputs.invert_xfm = True

    register_func_to_anat_use_T2.connect(linear_reg_func_to_t2, 'out_matrix_file', invt, 'in_file')

    # ${FSLDIR}/bin/applywarp --interp=nn -i ${T1wFolder}/${T2wRestoreImageBrain} -r ${fMRIFolder}/${ScoutName}_gdc --premat="$fMRIFolder"/T2w2Scout.mat -o ${fMRIFolder}/Scout_brain_mask.nii.gz
    anat_to_func = pe.Node(interface=fsl.ApplyWarp(),
                           name='anat_to_func')
    anat_to_func.inputs.interp = 'nn'

    register_func_to_anat_use_T2.connect(inputspec, 'T2_brain', anat_to_func, 'in_file')
    register_func_to_anat_use_T2.connect(inputspec, 'func', anat_to_func, 'ref_file')
    register_func_to_anat_use_T2.connect(invt, 'out_file', anat_to_func, 'premat')

    # ${FSLDIR}/bin/fslmaths ${fMRIFolder}/Scout_brain_mask.nii.gz -bin ${fMRIFolder}/Scout_brain_mask.nii.gz
    func_brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                                  name=f'func_brain_mask')
    func_brain_mask.inputs.args = '-bin'

    register_func_to_anat_use_T2.connect(anat_to_func, 'out_file', func_brain_mask, 'in_file')

    # ${FSLDIR}/bin/fslmaths ${fMRIFolder}/${ScoutName}_gdc -mas ${fMRIFolder}/Scout_brain_mask.nii.gz ${fMRIFolder}/Scout_brain_dc.nii.gz
    func_brain = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='func_brain')
    func_brain.inputs.op_string = "-mas %s "

    register_func_to_anat_use_T2.connect(inputspec, 'func', func_brain, 'in_file')
    register_func_to_anat_use_T2.connect(func_brain_mask, 'out_file', func_brain, 'operand_files')

    # ## re-registering the maked brain to the T1 brain:
    # ${FSLDIR}/bin/flirt -interp spline -dof 6 -in ${fMRIFolder}/Scout_brain_dc.nii.gz -ref ${T1wFolder}/${T1wRestoreImageBrain} -omat "$fMRIFolder"/${ScoutName}_gdc2T1w_init.mat -out ${fMRIFolder}/${ScoutName}_gdc2T1w_brain_init -searchrx -30 30 -searchry -30 30 -searchrz -30 30 -cost mutualinfo
    linear_reg_func_to_t1 = pe.Node(interface=fsl.FLIRT(),
                         name='linear_reg_func_to_t1')
    linear_reg_func_to_t1.inputs.interp = 'spline'
    linear_reg_func_to_t1.inputs.cost = 'mutualinfo'
    linear_reg_func_to_t1.inputs.dof = 6
    linear_reg_func_to_t1.inputs.searchr_x = [30, 30]
    linear_reg_func_to_t1.inputs.searchr_y = [30, 30]
    linear_reg_func_to_t1.inputs.searchr_z = [30, 30]

    register_func_to_anat_use_T2.connect(func_brain, 'out_file', linear_reg_func_to_t1, 'in_file')

    register_func_to_anat_use_T2.connect(inputspec, 'T1_brain', linear_reg_func_to_t1, 'reference')

    # #taking out warpfield as it is not being made without a fieldmap.
    # ${FSLDIR}/bin/convertwarp --relout --rel -r ${T1wFolder}/${T2wRestoreImage} --postmat=${fMRIFolder}/${ScoutName}_gdc2T1w_init.mat -o ${fMRIFolder}/${ScoutName}_gdc2T1w_init_warp
    convert_warp = pe.Node(interface=fsl.ConvertWarp(), name='convert_warp')

    convert_warp.inputs.out_relwarp = True
    convert_warp.inputs.relwarp = True

    register_func_to_anat_use_T2.connect(linear_reg_func_to_t1, 'out_matrix_file', convert_warp, 'postmat')

    register_func_to_anat_use_T2.connect(inputspec, 'T2_head', convert_warp, 'reference')


    register_func_to_anat_use_T2.connect(linear_reg_func_to_t1, 'out_matrix_file',
                                         outputspec,
                                         'func_to_anat_linear_xfm_nobbreg')

    register_func_to_anat_use_T2.connect(convert_warp, 'out_file',
                                         outputspec,
                                         'func_to_anat_linear_warp_nobbreg')

    register_func_to_anat_use_T2.connect(linear_reg_func_to_t1, 'out_file',
                                         outputspec, 'anat_func_nobbreg')

    return register_func_to_anat_use_T2


def create_bbregister_func_to_anat(phase_diff_distcor=False,
                                   name='bbregister_func_to_anat'):

    """
    Registers a functional scan in native space to structural.  This is
    meant to be used after create_nonlinear_register() has been run and
    relies on some of its outputs.

    Parameters
    ----------
    fieldmap_distortion : bool, optional
        If field map-based distortion correction is being run, FLIRT should
        take in the appropriate field map-related inputs.
    name : string, optional
        Name of the workflow.

    Returns
    -------
    register_func_to_anat : nipype.pipeline.engine.Workflow

    Notes
    -----

    Workflow Inputs::

        inputspec.func : string (nifti file)
            Input functional scan to be registered to MNI space
        inputspec.anat : string (nifti file)
            Corresponding full-head or brain scan of subject
        inputspec.linear_reg_matrix : string (mat file)
            Affine matrix from linear functional to anatomical registration
        inputspec.anat_wm_segmentation : string (nifti file)
            White matter segmentation probability mask in anatomical space
        inputspec.bbr_schedule : string (.sch file)
            Boundary based registration schedule file for flirt command

    Workflow Outputs::

        outputspec.func_to_anat_linear_xfm : string (mat file)
            Affine transformation from functional to anatomical native space
        outputspec.anat_func : string (nifti file)
            Functional data in anatomical space
    """

    register_bbregister_func_to_anat = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['func',
                                                       'anat',
                                                       'linear_reg_matrix',
                                                       'anat_wm_segmentation',
                                                       'bbr_schedule',
                                                       'bbr_wm_mask_args',
                                                       'fieldmap',
                                                       'fieldmapmask']),
                        name='inputspec')

    inputNode_echospacing = pe.Node(
        util.IdentityInterface(fields=['echospacing']),
        name='echospacing_input')

    inputNode_pedir = pe.Node(util.IdentityInterface(fields=['pedir']),
                              name='pedir_input')

    outputspec = pe.Node(util.IdentityInterface(
        fields=['func_to_anat_linear_xfm', 'anat_func']), name='outputspec')

    wm_bb_mask = pe.Node(interface=fsl.ImageMaths(),
                         name='wm_bb_mask')

    register_bbregister_func_to_anat.connect(
        inputspec, 'bbr_wm_mask_args',
        wm_bb_mask, 'op_string')

    register_bbregister_func_to_anat.connect(inputspec,
                                             'anat_wm_segmentation',
                                             wm_bb_mask, 'in_file')

    def bbreg_args(bbreg_target):
        return '-cost bbr -wmseg ' + bbreg_target

    bbreg_func_to_anat = pe.Node(interface=fsl.FLIRT(),
                                 name='bbreg_func_to_anat')
    bbreg_func_to_anat.inputs.dof = 6

    register_bbregister_func_to_anat.connect(
        inputspec, 'bbr_schedule',
        bbreg_func_to_anat, 'schedule')

    register_bbregister_func_to_anat.connect(
        wm_bb_mask, ('out_file', bbreg_args),
        bbreg_func_to_anat, 'args')

    register_bbregister_func_to_anat.connect(
        inputspec, 'func',
        bbreg_func_to_anat, 'in_file')

    register_bbregister_func_to_anat.connect(
        inputspec, 'anat',
        bbreg_func_to_anat, 'reference')

    register_bbregister_func_to_anat.connect(
        inputspec, 'linear_reg_matrix',
        bbreg_func_to_anat, 'in_matrix_file')

    if phase_diff_distcor:
        conv_pedir = \
            pe.Node(interface=util.Function(input_names=['pedir',
                                                         'convert'],
                                            output_names=['pedir'],
                                            function=convert_pedir),
                    name='bbreg_convert_pedir')
        conv_pedir.inputs.convert = 'xyz_to_int'

        register_bbregister_func_to_anat.connect(inputNode_pedir, 'pedir',
                                                 conv_pedir, 'pedir')
        register_bbregister_func_to_anat.connect(conv_pedir, 'pedir', 
                                                 bbreg_func_to_anat, 'pedir')
        register_bbregister_func_to_anat.connect(
            inputspec, 'fieldmap',
            bbreg_func_to_anat, 'fieldmap')
        register_bbregister_func_to_anat.connect(
            inputspec, 'fieldmapmask',
            bbreg_func_to_anat, 'fieldmapmask')
        register_bbregister_func_to_anat.connect(
            inputNode_echospacing, 'echospacing',
            bbreg_func_to_anat, 'echospacing')

    register_bbregister_func_to_anat.connect(
        bbreg_func_to_anat, 'out_matrix_file',
        outputspec, 'func_to_anat_linear_xfm')

    register_bbregister_func_to_anat.connect(
        bbreg_func_to_anat, 'out_file',
        outputspec, 'anat_func')

    return register_bbregister_func_to_anat


def create_wf_calculate_ants_warp(
    name='create_wf_calculate_ants_warp', num_threads=1, reg_ants_skull=1
):
    '''
    Calculates the nonlinear ANTS registration transform. This workflow
    employs the antsRegistration tool:

    http://stnava.github.io/ANTs/


    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    calc_ants_warp_wf : nipype.pipeline.engine.Workflow

    Notes
    -----

    Some of the inputs listed below are lists or lists of lists. This is
    because antsRegistration can perform multiple stages of calculations
    depending on how the user configures their registration.

    For example, if one wants to employ a different metric (with different
    parameters) at each stage, the lists would be configured like this:

    warp_wf.inputs.inputspec.transforms = ['Rigid','Affine','SyN']
    warp_wf.inputs.inputspec.transform_parameters = [[0.1],[0.1],[0.1,3,0]]

    ..where each element in the first list is a metric to be used at each
    stage, 'Rigid' being for stage 1, 'Affine' for stage 2, etc. The lists
    within the list for transform_parameters would then correspond to each
    stage's metric, with [0.1] applying to 'Rigid' and 'Affine' (stages 1 and
    2), and [0.1,3,0] applying to 'SyN' of stage 3.

    In some cases, when a parameter is not needed for a stage, 'None' must be
    entered in its place if there are other parameters for other stages.


    Workflow Inputs::

        inputspec.moving_brain : string (nifti file)
            File of brain to be normalized (registered)
        inputspec.reference_brain : string (nifti file)
            Target brain file to normalize to
        inputspec.dimension : integer
            Dimension of the image (default: 3)
        inputspec.use_histogram_matching : boolean
            Histogram match the images before registration
        inputspec.winsorize_lower_quantile : float
            Winsorize data based on quantiles (lower range)
        inputspec.winsorize_higher_quantile : float
            Winsorize data based on quantiles (higher range)
        inputspec.metric : list of strings
            Image metric(s) to be used at each stage
        inputspec.metric_weight : list of floats
            Modulate the per-stage weighting of the corresponding metric
        inputspec.radius_or_number_of_bins : list of integers
            Number of bins in each stage for the MI and Mattes metric, the
            radius for other metrics
        inputspec.sampling_strategy : list of strings
            Sampling strategy (or strategies) to use for the metrics
            {None, Regular, or Random}
        inputspec.sampling_percentage : list of floats
            Defines the sampling strategy
            {float value, or None}
        inputspec.number_of_iterations : list of lists of integers
            Determines the convergence
        inputspec.convergence_threshold : list of floats
            Threshold compared to the slope of the line fitted in convergence
        inputspec.convergence_window_size : list of integers
            Window size of convergence calculations
        inputspec.transforms : list of strings
            Selection of transform options. See antsRegistration documentation
            for a full list of options and their descriptions
        inputspec.transform_parameters : list of lists of floats
            Fine-tuning for the different transform options
        inputspec.shrink_factors : list of lists of integers
            Specify the shrink factor for the virtual domain (typically the
            fixed image) at each level
        inputspec.smoothing_sigmas : list of lists of floats
            Specify the sigma of gaussian smoothing at each level
        inputspec.fixed_image_mask: (an existing file name)
            Mask used to limit metric sampling region of the fixed imagein all
            stages
        inputspec.interp : string
            Type of interpolation to use
            ('Linear' or 'BSpline' or 'LanczosWindowedSinc')

    Workflow Outputs::

        outputspec.warp_field : string (nifti file)
            Output warp field of registration
        outputspec.inverse_warp_field : string (nifti file)
            Inverse of the warp field of the registration
        outputspec.ants_affine_xfm : string (.mat file)
            The affine matrix of the registration
        outputspec.ants_inverse_affine_xfm : string (.mat file)
            The affine matrix of the reverse registration
        outputspec.composite_transform : string (nifti file)
            The combined transform including the warp field and rigid & affine
            linear warps
        outputspec.normalized_output_brain : string (nifti file)
            Template-registered version of input brain

    Registration Procedure:

    1. Calculates a nonlinear anatomical-to-template registration.

    .. exec::
        from CPAC.registration import create_wf_calculate_ants_warp
        wf = create_wf_calculate_ants_warp()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/calculate_ants_warp.dot'
        )

    Workflow Graph:
    .. image::
        :width: 500

    Detailed Workflow Graph:

    .. image::
        :width: 500
    '''

    calc_ants_warp_wf = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(
        fields=['moving_brain',
                'reference_brain',
                'moving_skull',
                'reference_skull',
                'reference_mask',
                'moving_mask',
                'fixed_image_mask',
                'ants_para',
                'interp']),
                name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(
        fields=['ants_initial_xfm',
                'ants_rigid_xfm',
                'ants_affine_xfm',
                'warp_field',
                'inverse_warp_field',
                'composite_transform',
                'wait',
                'normalized_output_brain']), name='outputspec')

    # use ANTS to warp the masked anatomical image to a template image
    '''
    calculate_ants_warp = pe.Node(interface=ants.Registration(),
            name='calculate_ants_warp')

    calculate_ants_warp.inputs.output_warped_image = True
    calculate_ants_warp.inputs.initial_moving_transform_com = 0
    '''
    reg_imports = ['import os', 'import subprocess']
    calculate_ants_warp = \
        pe.Node(interface=util.Function(input_names=['moving_brain',
                                                     'reference_brain',
                                                     'moving_skull',
                                                     'reference_skull',
                                                     'ants_para',
                                                     'moving_mask',
                                                     'reference_mask',
                                                     'fixed_image_mask',
                                                     'interp',
                                                     'reg_with_skull'],
                                        output_names=['warp_list',
                                                      'warped_image'],
                                        function=hardcoded_reg,
                                        imports=reg_imports),
                name='calc_ants_warp',
                mem_gb=2.8,
                mem_x=(2e-7, 'moving_brain', 'xyz'),
                throttle=True)

    calculate_ants_warp.interface.num_threads = num_threads

    select_forward_initial = pe.Node(util.Function(
        input_names=['warp_list', 'selection'],
        output_names=['selected_warp'],
        function=seperate_warps_list), name='select_forward_initial')

    select_forward_initial.inputs.selection = "Initial"

    select_forward_rigid = pe.Node(util.Function(
        input_names=['warp_list', 'selection'],
        output_names=['selected_warp'],
        function=seperate_warps_list), name='select_forward_rigid')

    select_forward_rigid.inputs.selection = "Rigid"

    select_forward_affine = pe.Node(util.Function(
        input_names=['warp_list', 'selection'],
        output_names=['selected_warp'],
        function=seperate_warps_list), name='select_forward_affine')

    select_forward_affine.inputs.selection = "Affine"

    select_forward_warp = pe.Node(util.Function(
        input_names=['warp_list', 'selection'],
        output_names=['selected_warp'],
        function=seperate_warps_list), name='select_forward_warp')

    select_forward_warp.inputs.selection = "Warp"

    select_inverse_warp = pe.Node(util.Function(
        input_names=['warp_list', 'selection'],
        output_names=['selected_warp'],
        function=seperate_warps_list), name='select_inverse_warp')

    select_inverse_warp.inputs.selection = "Inverse"

    calc_ants_warp_wf.connect(
        inputspec, 'moving_brain',
        calculate_ants_warp, 'moving_brain')

    calc_ants_warp_wf.connect(
        inputspec, 'reference_brain',
        calculate_ants_warp, 'reference_brain')

    if reg_ants_skull == 1:

        calculate_ants_warp.inputs.reg_with_skull = 1

        calc_ants_warp_wf.connect(
            inputspec, 'moving_skull',
            calculate_ants_warp, 'moving_skull')

        calc_ants_warp_wf.connect(
            inputspec, 'reference_skull',
            calculate_ants_warp, 'reference_skull')

    else:
        calc_ants_warp_wf.connect(
            inputspec, 'moving_brain',
            calculate_ants_warp, 'moving_skull')

        calc_ants_warp_wf.connect(
            inputspec, 'reference_brain',
            calculate_ants_warp, 'reference_skull')

    calc_ants_warp_wf.connect(
        inputspec, 'fixed_image_mask',
        calculate_ants_warp, 'fixed_image_mask')

    calc_ants_warp_wf.connect(inputspec, 'reference_mask',
            calculate_ants_warp, 'reference_mask')

    calc_ants_warp_wf.connect(inputspec, 'moving_mask',
            calculate_ants_warp, 'moving_mask')

    calc_ants_warp_wf.connect(inputspec, 'ants_para',
            calculate_ants_warp, 'ants_para')

    calc_ants_warp_wf.connect(
        inputspec, 'interp',
        calculate_ants_warp, 'interp')

    # inter-workflow connections

    calc_ants_warp_wf.connect(
        calculate_ants_warp, 'warp_list',
        select_forward_initial, 'warp_list')

    calc_ants_warp_wf.connect(
        calculate_ants_warp, 'warp_list',
        select_forward_rigid, 'warp_list')

    calc_ants_warp_wf.connect(
        calculate_ants_warp, 'warp_list',
        select_forward_affine, 'warp_list')

    calc_ants_warp_wf.connect(
        calculate_ants_warp, 'warp_list',
        select_forward_warp, 'warp_list')

    calc_ants_warp_wf.connect(
        calculate_ants_warp, 'warp_list',
        select_inverse_warp, 'warp_list')

    # connections to outputspec

    calc_ants_warp_wf.connect(
        select_forward_initial, 'selected_warp',
        outputspec, 'ants_initial_xfm')

    calc_ants_warp_wf.connect(
        select_forward_rigid, 'selected_warp',
        outputspec, 'ants_rigid_xfm')

    calc_ants_warp_wf.connect(
        select_forward_affine, 'selected_warp',
        outputspec, 'ants_affine_xfm')

    calc_ants_warp_wf.connect(
        select_forward_warp, 'selected_warp',
        outputspec, 'warp_field')

    calc_ants_warp_wf.connect(
        select_inverse_warp, 'selected_warp',
        outputspec, 'inverse_warp_field')

    calc_ants_warp_wf.connect(
        calculate_ants_warp, 'warped_image',
        outputspec, 'normalized_output_brain')

    return calc_ants_warp_wf


def FSL_registration_connector(wf_name, cfg, orig="T1w", opt=None,
                               symmetric=False, template="T1w"):

    wf = pe.Workflow(name=wf_name)

    inputNode = pe.Node(
        util.IdentityInterface(fields=['input_brain',
                                       'reference_brain',
                                       'input_head',
                                       'reference_head',
                                       'input_mask',
                                       'reference_mask',
                                       'transform',
                                       'interpolation',
                                       'fnirt_config']),
        name='inputspec')

    sym = ''
    symm = ''
    if symmetric:
        sym = 'sym'
        symm = '_symmetric'

    tmpl = ''
    if template == 'EPI':
        tmpl = 'EPI'

    if opt == 'FSL' or opt == 'FSL-linear':

        flirt_reg_anat_mni = create_fsl_flirt_linear_reg(
            f'anat_mni_flirt_register{symm}'
        )

        # Input registration parameters
        wf.connect(inputNode, 'interpolation',
                   flirt_reg_anat_mni, 'inputspec.interp')

        wf.connect(inputNode, 'input_brain',
                   flirt_reg_anat_mni, 'inputspec.input_brain')

        wf.connect(inputNode, 'reference_brain', flirt_reg_anat_mni,
                   'inputspec.reference_brain')

        write_lin_composite_xfm = pe.Node(interface=fsl.ConvertWarp(),
                                          name=f'fsl_lin-warp_to_nii{symm}')

        wf.connect(inputNode, 'reference_brain',
                   write_lin_composite_xfm, 'reference')

        wf.connect(flirt_reg_anat_mni, 'outputspec.linear_xfm',
                   write_lin_composite_xfm, 'premat')

        write_invlin_composite_xfm = pe.Node(interface=fsl.ConvertWarp(),
                                             name=f'fsl_invlin-warp_to_'
                                                  f'nii{symm}')

        wf.connect(inputNode, 'reference_brain',
                   write_invlin_composite_xfm, 'reference')

        wf.connect(flirt_reg_anat_mni, 'outputspec.invlinear_xfm',
                   write_invlin_composite_xfm, 'premat')

        outputs = {
            f'space-{sym}template_desc-preproc_{orig}': (
                flirt_reg_anat_mni, 'outputspec.output_brain'),
            f'from-{orig}_to-{sym}{tmpl}template_mode-image_desc-linear_xfm': (
                write_lin_composite_xfm, 'out_file'),
            f'from-{sym}{tmpl}template_to-{orig}_mode-image_desc-linear_xfm': (
                write_invlin_composite_xfm, 'out_file'),
            f'from-{orig}_to-{sym}{tmpl}template_mode-image_xfm': (
                write_lin_composite_xfm, 'out_file')
        }


    if opt == 'FSL':
        if cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['ref_resolution'] ==  \
            cfg.registration_workflows['anatomical_registration']['resolution_for_anat']:
            fnirt_reg_anat_mni = create_fsl_fnirt_nonlinear_reg(
                f'anat_mni_fnirt_register{symm}'
            )
        else:
            fnirt_reg_anat_mni = create_fsl_fnirt_nonlinear_reg_nhp(
                f'anat_mni_fnirt_register{symm}'
            )

        wf.connect(inputNode, 'input_brain',
                   fnirt_reg_anat_mni, 'inputspec.input_brain')

        wf.connect(inputNode, 'reference_brain',
                   fnirt_reg_anat_mni, 'inputspec.reference_brain')

        wf.connect(inputNode, 'input_head',
                   fnirt_reg_anat_mni, 'inputspec.input_skull')

        # NOTE: crossover from above opt block
        wf.connect(flirt_reg_anat_mni, 'outputspec.linear_xfm',
                   fnirt_reg_anat_mni, 'inputspec.linear_aff')

        wf.connect(inputNode, 'reference_head',
                   fnirt_reg_anat_mni, 'inputspec.reference_skull')

        wf.connect(inputNode, 'reference_mask',
                   fnirt_reg_anat_mni, 'inputspec.ref_mask')

        # assign the FSL FNIRT config file specified in pipeline config.yml
        wf.connect(inputNode, 'fnirt_config',
                   fnirt_reg_anat_mni, 'inputspec.fnirt_config')

        if cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['ref_resolution'] ==  \
            cfg.registration_workflows['anatomical_registration']['resolution_for_anat']:
            # NOTE: this is an UPDATE because of the opt block above
            added_outputs = {
                f'space-{sym}template_desc-preproc_{orig}': (
                    fnirt_reg_anat_mni, 'outputspec.output_brain'),
                f'from-{orig}_to-{sym}{tmpl}template_mode-image_xfm': (
                    fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm')
            }
            outputs.update(added_outputs)
        else:
            # NOTE: this is an UPDATE because of the opt block above
            added_outputs = {
                f'space-{sym}template_desc-preproc_{orig}': (
                    fnirt_reg_anat_mni, 'outputspec.output_brain'),
                f'space-{sym}template_desc-head_{orig}': (
                    fnirt_reg_anat_mni, 'outputspec.output_head'),
                f'space-{sym}template_desc-{orig}_mask': (
                    fnirt_reg_anat_mni, 'outputspec.output_mask'),
                f'space-{sym}template_desc-T1wT2w_biasfield': (
                    fnirt_reg_anat_mni, 'outputspec.output_biasfield'),
                f'from-{orig}_to-{sym}{tmpl}template_mode-image_xfm': (
                    fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm'),
                f'from-{orig}_to-{sym}{tmpl}template_mode-image_warp': (
                    fnirt_reg_anat_mni, 'outputspec.nonlinear_warp')
            }
            outputs.update(added_outputs)

    return (wf, outputs)


def ANTs_registration_connector(wf_name, cfg, params, orig="T1w",
                                symmetric=False, template="T1w"):

    wf = pe.Workflow(name=wf_name)

    inputNode = pe.Node(
        util.IdentityInterface(fields=['input_brain',
                                       'reference_brain',
                                       'input_head',
                                       'reference_head',
                                       'input_mask',
                                       'reference_mask',
                                       'transform',
                                       'interpolation']),
        name='inputspec')

    sym = ''
    symm = ''
    if symmetric:
        sym = 'sym'
        symm = '_symmetric'

    tmpl = ''
    if template == 'EPI':
        tmpl = 'EPI'

    if params is None:
        err_msg = '\n\n[!] C-PAC says: \nYou have selected ANTs as your ' \
                  'anatomical registration method.\n' \
                  'However, no ANTs parameters were specified.\n' \
                  'Please specify ANTs parameters properly and try again.'
        raise Exception(err_msg)

    ants_reg_anat_mni = \
        create_wf_calculate_ants_warp(
            f'anat_mni_ants_register{symm}',
            num_threads=cfg.pipeline_setup['system_config'][
                'num_ants_threads'],
            reg_ants_skull=cfg['registration_workflows'][
                'anatomical_registration']['reg_with_skull']
        )
    ants_reg_anat_mni.inputs.inputspec.ants_para = params

    wf.connect(inputNode, 'interpolation',
               ants_reg_anat_mni, 'inputspec.interp')

    # calculating the transform with the skullstripped is
    # reported to be better, but it requires very high
    # quality skullstripping. If skullstripping is imprecise
    # registration with skull is preferred

    wf.connect(inputNode, 'input_brain',
               ants_reg_anat_mni, 'inputspec.moving_brain')

    wf.connect(inputNode, 'reference_brain',
               ants_reg_anat_mni, 'inputspec.reference_brain')

    wf.connect(inputNode, 'input_head',
               ants_reg_anat_mni, 'inputspec.moving_skull')

    wf.connect(inputNode, 'reference_head',
               ants_reg_anat_mni, 'inputspec.reference_skull')

    wf.connect(inputNode, 'input_mask',
               ants_reg_anat_mni, 'inputspec.moving_mask')

    wf.connect(inputNode, 'reference_mask',
               ants_reg_anat_mni, 'inputspec.reference_mask')

    ants_reg_anat_mni.inputs.inputspec.fixed_image_mask = None

    if orig == 'T1w':
        if cfg.registration_workflows['anatomical_registration'][
            'registration']['ANTs']['use_lesion_mask']:
            # Create lesion preproc node to apply afni Refit and Resample
            lesion_preproc = create_lesion_preproc(
                wf_name=f'lesion_preproc{symm}'
            )
            wf.connect(inputNode, 'lesion_mask',
                       lesion_preproc, 'inputspec.lesion')
            wf.connect(lesion_preproc, 'outputspec.reorient',
                       ants_reg_anat_mni, 'inputspec.fixed_image_mask')

    # combine the linear xfm's into one - makes it easier downstream
    write_composite_linear_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_linear{symm}_xfm',
        mem_gb=1.155,
        mem_x=(1708448960473801 / 1208925819614629174706176, 'input_image'))
    write_composite_linear_xfm.inputs.print_out_composite_warp_file = True
    write_composite_linear_xfm.inputs.output_image = \
        f"from-{orig}_to-{sym}{tmpl}template_mode-image_desc-linear_xfm.nii.gz"

    wf.connect(inputNode, 'input_brain',
               write_composite_linear_xfm, 'input_image')

    wf.connect(inputNode, 'reference_brain',
               write_composite_linear_xfm, 'reference_image')

    wf.connect(inputNode, 'interpolation',
               write_composite_linear_xfm, 'interpolation')

    write_composite_linear_xfm.inputs.input_image_type = 0
    write_composite_linear_xfm.inputs.dimension = 3

    collect_transforms = pe.Node(util.Merge(3),
                                 name=f'collect_transforms{symm}',
                                 mem_gb=0.8,
                                 mem_x=(263474863123069 /
                                        37778931862957161709568,
                                        'in1'))

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
               collect_transforms, 'in3')

    # check transform list to exclude Nonetype (missing) init/rig/affine
    check_transform = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['checked_transform_list',
                                    'list_length'],
                      function=check_transforms),
        name=f'check_transforms',
        mem_gb=6)

    wf.connect(collect_transforms, 'out', check_transform, 'transform_list')

    wf.connect(check_transform, 'checked_transform_list',
               write_composite_linear_xfm, 'transforms')

    # combine the linear xfm's into one - makes it easier downstream
    write_composite_invlinear_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_invlinear{symm}_xfm',
        mem_gb=1.05,
        mem_x=(1367826948979337 / 151115727451828646838272, 'input_image'))
    write_composite_invlinear_xfm.inputs.print_out_composite_warp_file = True
    write_composite_invlinear_xfm.inputs.output_image = \
        f"from-{sym}{tmpl}template_to-{orig}_mode-image_desc-linear_xfm.nii.gz"

    wf.connect(inputNode, 'reference_brain',
               write_composite_invlinear_xfm, 'input_image')

    wf.connect(inputNode, 'input_brain',
               write_composite_invlinear_xfm, 'reference_image')

    wf.connect(inputNode, 'interpolation',
               write_composite_invlinear_xfm, 'interpolation')

    write_composite_invlinear_xfm.inputs.input_image_type = 0
    write_composite_invlinear_xfm.inputs.dimension = 3

    collect_inv_transforms = pe.Node(util.Merge(3),
                                     name='collect_inv_transforms'
                                          f'{symm}')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
               collect_inv_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_inv_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_inv_transforms, 'in3')

    # check transform list to exclude Nonetype (missing) init/rig/affine
    check_invlinear_transform = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['checked_transform_list',
                                    'list_length'],
                      function=check_transforms),
        name=f'check_inv_transforms')

    wf.connect(collect_inv_transforms, 'out',
               check_invlinear_transform, 'transform_list')

    wf.connect(check_invlinear_transform, 'checked_transform_list',
               write_composite_invlinear_xfm, 'transforms')

    # generate inverse transform flags, which depends on the
    # number of transforms
    inverse_transform_flags = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['inverse_transform_flags'],
                      function=generate_inverse_transform_flags),
        name=f'inverse_transform_flags')

    wf.connect(check_invlinear_transform, 'checked_transform_list',
               inverse_transform_flags, 'transform_list')

    wf.connect(inverse_transform_flags, 'inverse_transform_flags',
               write_composite_invlinear_xfm, 'invert_transform_flags')

    # combine ALL xfm's into one - makes it easier downstream
    write_composite_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_{symm}xfm',
        mem_gb=1.5)
    write_composite_xfm.inputs.print_out_composite_warp_file = True
    write_composite_xfm.inputs.output_image = \
        f"from-{orig}_to-{sym}{tmpl}template_mode-image_xfm.nii.gz"

    wf.connect(inputNode, 'input_brain', write_composite_xfm, 'input_image')

    wf.connect(inputNode, 'reference_brain',
               write_composite_xfm, 'reference_image')

    wf.connect(inputNode, 'interpolation',
               write_composite_xfm, 'interpolation')

    write_composite_xfm.inputs.input_image_type = 0
    write_composite_xfm.inputs.dimension = 3

    collect_all_transforms = pe.Node(util.Merge(4),
                                 name=f'collect_all_transforms'
                                      f'{symm}')

    wf.connect(ants_reg_anat_mni, 'outputspec.warp_field',
               collect_all_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_all_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_all_transforms, 'in3')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
               collect_all_transforms, 'in4')

    # check transform list to exclude Nonetype (missing) init/rig/affine
    check_all_transform = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['checked_transform_list',
                                    'list_length'],
                      function=check_transforms),
        name=f'check_all_transforms')

    wf.connect(collect_all_transforms, 'out',
               check_all_transform, 'transform_list')

    wf.connect(check_all_transform, 'checked_transform_list',
               write_composite_xfm, 'transforms')

    # combine ALL xfm's into one - makes it easier downstream
    write_composite_inv_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_inv_{symm}xfm',
        mem_gb=0.3,
        mem_x=(6278549929741219 / 604462909807314587353088, 'input_image'))
    write_composite_inv_xfm.inputs.print_out_composite_warp_file = True
    write_composite_inv_xfm.inputs.output_image = \
        f"from-{sym}{tmpl}template_to-{orig}_mode-image_xfm.nii.gz"

    wf.connect(inputNode, 'reference_brain',
               write_composite_inv_xfm, 'input_image')

    wf.connect(inputNode, 'input_brain',
               write_composite_inv_xfm, 'reference_image')

    wf.connect(inputNode, 'interpolation',
               write_composite_inv_xfm, 'interpolation')

    write_composite_inv_xfm.inputs.input_image_type = 0
    write_composite_inv_xfm.inputs.dimension = 3

    collect_all_inv_transforms = pe.Node(util.Merge(4),
                                         name=f'collect_all_inv_transforms'
                                         f'{symm}')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
               collect_all_inv_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_all_inv_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_all_inv_transforms, 'in3')

    wf.connect(ants_reg_anat_mni, 'outputspec.inverse_warp_field',
               collect_all_inv_transforms, 'in4')

    # check transform list to exclude Nonetype (missing) init/rig/affine
    check_all_inv_transform = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['checked_transform_list',
                                    'list_length'],
                      function=check_transforms),
        name=f'check_all_inv_transforms')

    wf.connect(collect_all_inv_transforms, 'out',
               check_all_inv_transform, 'transform_list')

    wf.connect(check_all_inv_transform, 'checked_transform_list',
               write_composite_inv_xfm, 'transforms')

    # generate inverse transform flags, which depends on the
    # number of transforms
    inverse_all_transform_flags = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['inverse_transform_flags'],
                      function=generate_inverse_transform_flags),
        name=f'inverse_all_transform_flags')

    wf.connect(check_all_inv_transform, 'checked_transform_list',
               inverse_all_transform_flags, 'transform_list')

    wf.connect(inverse_all_transform_flags, 'inverse_transform_flags',
               write_composite_inv_xfm, 'invert_transform_flags')

    outputs = {
        f'space-{sym}template_desc-preproc_{orig}': (
            ants_reg_anat_mni, 'outputspec.normalized_output_brain'),
        f'from-{orig}_to-{sym}{tmpl}template_mode-image_xfm': (
            write_composite_xfm, 'output_image'),
        f'from-{sym}{tmpl}template_to-{orig}_mode-image_xfm': (
            write_composite_inv_xfm, 'output_image'),
        f'from-{orig}_to-{sym}{tmpl}template_mode-image_desc-linear_xfm': (
            write_composite_linear_xfm, 'output_image'),
        f'from-{sym}{tmpl}template_to-{orig}_mode-image_desc-linear_xfm': (
            write_composite_invlinear_xfm, 'output_image'),
        f'from-{orig}_to-{sym}{tmpl}template_mode-image_desc-nonlinear_xfm': (
            ants_reg_anat_mni, 'outputspec.warp_field'),
        f'from-{sym}{tmpl}template_to-{orig}_mode-image_desc-nonlinear_xfm': (
            ants_reg_anat_mni, 'outputspec.inverse_warp_field')
    }

    return (wf, outputs)


def bold_to_T1template_xfm_connector(wf_name, cfg, reg_tool, symmetric=False,
                                     blip=False):

    wf = pe.Workflow(name=wf_name)

    inputNode = pe.Node(
        util.IdentityInterface(fields=['input_brain',
                                       'mean_bold',
                                       'coreg_xfm',
                                       'T1w-brain-template_funcreg',
                                       'T1w_to_template_xfm',
                                       'template_to_T1w_xfm',
                                       'blip_warp']),
        name='inputspec')

    sym = ''
    if symmetric:
        sym = 'sym'

    if reg_tool == 'ants':
        fsl_reg_2_itk = pe.Node(c3.C3dAffineTool(), name='fsl_reg_2_itk')
        fsl_reg_2_itk.inputs.itk_transform = True
        fsl_reg_2_itk.inputs.fsl2ras = True

        # convert the .mat from linear Func->Anat to
        # ANTS format
        wf.connect(inputNode, 'coreg_xfm', fsl_reg_2_itk, 'transform_file')

        wf.connect(inputNode, 'input_brain', fsl_reg_2_itk, 'reference_file')

        wf.connect(inputNode, 'mean_bold', fsl_reg_2_itk, 'source_file')

        itk_imports = ['import os']
        change_transform = pe.Node(util.Function(
            input_names=['input_affine_file'],
            output_names=['updated_affine_file'],
            function=change_itk_transform_type,
            imports=itk_imports),
            name='change_transform_type')

        wf.connect(fsl_reg_2_itk, 'itk_transform',
                         change_transform, 'input_affine_file')

        # combine ALL xfm's into one - makes it easier downstream
        write_composite_xfm = pe.Node(
            interface=ants.ApplyTransforms(),
            name=f'write_composite_xfm',
            mem_gb=1.5)
        write_composite_xfm.inputs.print_out_composite_warp_file = True
        write_composite_xfm.inputs.output_image = \
            f"from-bold_to-{sym}template_mode-image_xfm.nii.gz"

        wf.connect(inputNode, 'mean_bold',
                   write_composite_xfm, 'input_image')

        wf.connect(inputNode, 'T1w-brain-template_funcreg',
                   write_composite_xfm, 'reference_image')

        write_composite_xfm.inputs.input_image_type = 0
        write_composite_xfm.inputs.dimension = 3
        write_composite_xfm.inputs.interpolation = \
            cfg.registration_workflows['anatomical_registration'][
                'registration']['ANTs']['interpolation']

        if not blip:
            collect_all_transforms = pe.Node(util.Merge(2),
                                             name='collect_all_transforms')
        else:
            collect_all_transforms = pe.Node(util.Merge(3),
                                             name='collect_all_transforms')

            wf.connect(inputNode, 'blip_warp',
                       collect_all_transforms, 'in3')

        wf.connect(inputNode, 'T1w_to_template_xfm',
                  collect_all_transforms, 'in1')

        wf.connect(change_transform, 'updated_affine_file',
                   collect_all_transforms, 'in2')

        wf.connect(collect_all_transforms, 'out',
                   write_composite_xfm, 'transforms')

        write_composite_inv_xfm = pe.Node(
            interface=ants.ApplyTransforms(),
            name=f'write_composite_inv_xfm',
            mem_gb=1.5)
        write_composite_inv_xfm.inputs.print_out_composite_warp_file = True
        write_composite_inv_xfm.inputs.invert_transform_flags = [True, False]
        write_composite_inv_xfm.inputs.output_image = \
            f"from-{sym}template_to-bold_mode-image_xfm.nii.gz"

        wf.connect(inputNode, 'T1w-brain-template_funcreg',
                   write_composite_inv_xfm, 'input_image')

        wf.connect(inputNode, 'mean_bold',
                   write_composite_inv_xfm, 'reference_image')

        write_composite_inv_xfm.inputs.input_image_type = 0
        write_composite_inv_xfm.inputs.dimension = 3
        write_composite_inv_xfm.inputs.interpolation = \
            cfg.registration_workflows['anatomical_registration'][
                'registration']['ANTs']['interpolation']

        collect_inv_transforms = pe.Node(util.Merge(2),
                                         name='collect_inv_transforms')

        wf.connect(change_transform, 'updated_affine_file',
                   collect_inv_transforms, 'in1')

        wf.connect(inputNode, 'template_to_T1w_xfm',
                   collect_inv_transforms, 'in2')

        wf.connect(collect_inv_transforms, 'out',
                   write_composite_inv_xfm, 'transforms')

        outputs = {
            f'from-bold_to-{sym}template_mode-image_xfm':
                (write_composite_xfm, 'output_image'),
            f'from-{sym}template_to-bold_mode-image_xfm':
                (write_composite_inv_xfm, 'output_image')
        }

    elif reg_tool == 'fsl':

        write_composite_xfm = pe.Node(interface=fsl.ConvertWarp(),
                                      name='combine_fsl_warps')

        wf.connect(inputNode, 'T1w-brain-template_funcreg',
                   write_composite_xfm, 'reference')

        if blip:
            wf.connect(inputNode, 'coreg_xfm', 
                       write_composite_xfm, 'postmat')
            wf.connect(inputNode, 'blip_warp', 
                       write_composite_xfm, 'warp1')
            wf.connect(inputNode, 'T1w_to_template_xfm', 
                       write_composite_xfm, 'warp2')
        else:
            wf.connect(inputNode, 'coreg_xfm', 
                       write_composite_xfm, 'premat')
            wf.connect(inputNode, 'T1w_to_template_xfm',
                       write_composite_xfm, 'warp1')

        outputs = {
            f'from-bold_to-{sym}template_mode-image_xfm':
                (write_composite_xfm, 'out_file'),
        }

    return (wf, outputs)


@nodeblock(
    name="register_FSL_anat_to_template",
    config=["registration_workflows", "anatomical_registration"],
    switch=["run"],
    option_key=["registration", "using"],
    option_val=["FSL", "FSL-linear"],
    inputs=[
        (
            ["desc-preproc_T1w", "space-longitudinal_desc-reorient_T1w"],
            ["desc-brain_T1w", "space-longitudinal_desc-brain_T1w"],
        ),
        "T1w-template",
        "T1w-brain-template",
        "FNIRT-T1w-template",
        "FNIRT-T1w-brain-template",
        "template-ref-mask",
    ],
    outputs={
        "space-template_desc-preproc_T1w": {"Template": "T1w-brain-template"},
        "space-template_desc-head_T1w": {"Template": "T1w-template"},
        "space-template_desc-T1w_mask": {"Template": "T1w-template"},
        "space-template_desc-T1wT2w_biasfield": {"Template": "T1w-template"},
        "from-T1w_to-template_mode-image_desc-linear_xfm": {
            "Template": "T1w-template"},
        "from-template_to-T1w_mode-image_desc-linear_xfm": {
            "Template": "T1w-template"},
        "from-T1w_to-template_mode-image_xfm": {"Template": "T1w-template"},
        "from-T1w_to-template_mode-image_warp": {"Template": "T1w-template"},
        "from-longitudinal_to-template_mode-image_desc-linear_xfm": {
            "Template": "T1w-template"
        },
        "from-template_to-longitudinal_mode-image_desc-linear_xfm": {
            "Template": "T1w-template"
        },
        "from-longitudinal_to-template_mode-image_xfm": {
            "Template": "T1w-template"},
    },
)
def register_FSL_anat_to_template(wf, cfg, strat_pool, pipe_num, opt=None):

    fsl, outputs = FSL_registration_connector(f'register_{opt}_anat_to_'
                                              f'template_{pipe_num}', cfg,
                                              orig='T1w', opt=opt)

    fsl.inputs.inputspec.interpolation = cfg.registration_workflows[
        'anatomical_registration']['registration']['FSL-FNIRT'][
        'interpolation']

    fsl.inputs.inputspec.fnirt_config = cfg.registration_workflows[
        'anatomical_registration']['registration']['FSL-FNIRT'][
        'fnirt_config']

    connect, brain = \
        strat_pool.get_data(['desc-brain_T1w',
                             'space-longitudinal_desc-brain_T1w'],
                            report_fetched=True)
    node, out = connect
    wf.connect(node, out, fsl, 'inputspec.input_brain')

    if cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['ref_resolution'] ==  \
        cfg.registration_workflows['anatomical_registration']['resolution_for_anat']:

        node, out = strat_pool.get_data('T1w-brain-template')
        wf.connect(node, out, fsl, 'inputspec.reference_brain')

        node, out = strat_pool.get_data('T1w-template')
        wf.connect(node, out, fsl, 'inputspec.reference_head')
    else:
        node, out = strat_pool.get_data('FNIRT-T1w-brain-template')
        wf.connect(node, out, fsl, 'inputspec.reference_brain')

        node, out = strat_pool.get_data('FNIRT-T1w-template')
        wf.connect(node, out, fsl, 'inputspec.reference_head')

    node, out = strat_pool.get_data(["desc-preproc_T1w",
                                     "space-longitudinal_desc-reorient_T1w"])
    wf.connect(node, out, fsl, 'inputspec.input_head')

    node, out = strat_pool.get_data('template-ref-mask')
    wf.connect(node, out, fsl, 'inputspec.reference_mask')

    if 'space-longitudinal' in brain:
        for key in outputs.keys():
            if 'from-T1w' in key:
                new_key = key.replace('from-T1w', 'from-longitudinal')
                outputs[new_key] = outputs[key]
                del outputs[key]
            if 'to-T1w' in key:
                new_key = key.replace('to-T1w', 'to-longitudinal')
                outputs[new_key] = outputs[key]
                del outputs[key]

    return (wf, outputs)


@nodeblock(
    name="register_symmetric_FSL_anat_to_template",
    config=["registration_workflows", "anatomical_registration"],
    switch=["run"],
    option_key=["registration", "using"],
    option_val=["FSL", "FSL-linear"],
    inputs=[
        (
            ["desc-preproc_T1w", "space-longitudinal_desc-reorient_T1w"],
            ["desc-brain_T1w", "space-longitudinal_desc-brain_T1w"],
        ),
        "T1w-template-symmetric",
        "T1w-brain-template-symmetric",
        "dilated-symmetric-brain-mask",
    ],
    outputs={
        "space-symtemplate_desc-preproc_T1w": {
            "Template": "T1w-brain-template-symmetric"
        },
        "from-T1w_to-symtemplate_mode-image_desc-linear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-symtemplate_to-T1w_mode-image_desc-linear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-T1w_to-symtemplate_mode-image_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-longitudinal_to-symtemplate_mode-image_desc-linear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-symtemplate_to-longitudinal_mode-image_desc-linear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-longitudinal_to-symtemplate_mode-image_xfm": {
            "Template": "T1w-template-symmetric"
        },
    },
)
def register_symmetric_FSL_anat_to_template(wf, cfg, strat_pool, pipe_num,
                                            opt=None):

    fsl, outputs = FSL_registration_connector(f'register_{opt}_anat_to_'
                                              f'template_symmetric_'
                                              f'{pipe_num}', cfg, orig='T1w',
                                              opt=opt, symmetric=True)

    fsl.inputs.inputspec.interpolation = cfg.registration_workflows[
        'anatomical_registration']['registration']['FSL-FNIRT'][
        'interpolation']

    fsl.inputs.inputspec.fnirt_config = cfg.registration_workflows[
        'anatomical_registration']['registration']['FSL-FNIRT'][
        'fnirt_config']

    connect, brain = \
        strat_pool.get_data(['desc-brain_T1w',
                             'space-longitudinal_desc-brain_T1w'],
                            report_fetched=True)
    node, out = connect
    wf.connect(node, out, fsl, 'inputspec.input_brain')

    node, out = strat_pool.get_data('T1w-brain-template-symmetric')
    wf.connect(node, out, fsl, 'inputspec.reference_brain')

    node, out = strat_pool.get_data(["desc-preproc_T1w",
                                     "space-longitudinal_desc-reorient_T1w"])
    wf.connect(node, out, fsl, 'inputspec.input_head')

    node, out = strat_pool.get_data('T1w-template-symmetric')
    wf.connect(node, out, fsl, 'inputspec.reference_head')

    node, out = strat_pool.get_data('dilated-symmetric-brain-mask')
    wf.connect(node, out, fsl, 'inputspec.reference_mask')

    if 'space-longitudinal' in brain:
        for key in outputs.keys():
            if 'from-T1w' in key:
                new_key = key.replace('from-T1w', 'from-longitudinal')
                outputs[new_key] = outputs[key]
                del outputs[key]
            if 'to-T1w' in key:
                new_key = key.replace('to-T1w', 'to-longitudinal')
                outputs[new_key] = outputs[key]
                del outputs[key]

    return (wf, outputs)


@nodeblock(
    name="register_FSL_EPI_to_template",
    config=["registration_workflows", "functional_registration",
            "EPI_registration"],
    switch=["run"],
    option_key="using",
    option_val=["FSL", "FSL-linear"],
    inputs=[
        ("sbref", "space-bold_desc-brain_mask"),
        "EPI-template",
        "EPI-template-mask",
    ],
    outputs={
        "space-template_desc-preproc_bold": {"Template": "EPI-template"},
        "from-bold_to-EPItemplate_mode-image_desc-linear_xfm": {
            "Template": "EPI-template"
        },
        "from-EPItemplate_to-bold_mode-image_desc-linear_xfm": {
            "Template": "EPI-template"
        },
        "from-bold_to-EPItemplate_mode-image_xfm": {
            "Template": "EPI-template"},
    },
)
def register_FSL_EPI_to_template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Directly register the mean functional to an EPI template. No T1w
    involved.
    '''

    fsl, outputs = FSL_registration_connector(f'register_{opt}_EPI_to_'
                                              f'template_{pipe_num}', cfg,
                                              orig='bold', opt=opt,
                                              template='EPI')

    fsl.inputs.inputspec.interpolation = cfg['registration_workflows'][
        'functional_registration']['EPI_registration']['FSL-FNIRT'][
        'interpolation']

    fsl.inputs.inputspec.fnirt_config = cfg['registration_workflows'][
        'functional_registration']['EPI_registration']['FSL-FNIRT'][
        'fnirt_config']

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, fsl, 'inputspec.input_brain')

    node, out = strat_pool.get_data('EPI-template')
    wf.connect(node, out, fsl, 'inputspec.reference_brain')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, fsl, 'inputspec.input_head')

    node, out = strat_pool.get_data('EPI-template')
    wf.connect(node, out, fsl, 'inputspec.reference_head')

    node, out = strat_pool.get_data('EPI-template-mask')
    wf.connect(node, out, fsl, 'inputspec.reference_mask')

    return (wf, outputs)


@nodeblock(
    name="register_ANTs_anat_to_template",
    config=["registration_workflows", "anatomical_registration"],
    switch=["run"],
    option_key=["registration", "using"],
    option_val="ANTS",
    inputs=[
        (
            ["desc-preproc_T1w", "space-longitudinal_desc-brain_T1w"],
            [
                "space-T1w_desc-brain_mask",
                "space-longitudinal_desc-brain_mask",
                "space-T1w_desc-acpcbrain_mask",
            ],
            [
                "desc-restore_T1w",
                "desc-head_T1w",
                "desc-preproc_T1w",
                "space-longitudinal_desc-reorient_T1w",
            ],
            "space-template_desc-head_T1w",
            "space-template_desc-preproc_T1w",
        ),
        "T1w-template",
        "T1w-brain-template",
        "T1w-brain-template-mask",
        "label-lesion_mask",
    ],
    outputs={
        "space-template_desc-preproc_T1w": {
            "Description": "The preprocessed T1w brain transformed to "
                           "template space.",
            "Template": "T1w-template",
        },
        "from-T1w_to-template_mode-image_desc-linear_xfm": {
            "Description": "Linear (affine) transform from T1w native space "
                           "to T1w-template space.",
            "Template": "T1w-template",
        },
        "from-template_to-T1w_mode-image_desc-linear_xfm": {
            "Description": "Linear (affine) transform from T1w-template space "
                           "to T1w native space.",
            "Template": "T1w-template",
        },
        "from-T1w_to-template_mode-image_desc-nonlinear_xfm": {
            "Description": "Nonlinear (warp field) transform from T1w native "
                           "space to T1w-template space.",
            "Template": "T1w-template",
        },
        "from-template_to-T1w_mode-image_desc-nonlinear_xfm": {
            "Description": "Nonlinear (warp field) transform from "
                           "T1w-template space to T1w native space.",
            "Template": "T1w-template",
        },
        "from-T1w_to-template_mode-image_xfm": {
            "Description": "Composite (affine + warp field) transform from "
                           "T1w native space to T1w-template space.",
            "Template": "T1w-template",
        },
        "from-template_to-T1w_mode-image_xfm": {
            "Description": "Composite (affine + warp field) transform from "
                           "T1w-template space to T1w native space.",
            "Template": "T1w-template",
        },
        "from-longitudinal_to-template_mode-image_desc-linear_xfm": {
            "Description": "Linear (affine) transform from "
                           "longitudinal-template space to T1w-template "
                           "space.",
            "Template": "T1w-template",
        },
        "from-template_to-longitudinal_mode-image_desc-linear_xfm": {
            "Description": "Linear (affine) transform from T1w-template "
                           "space to longitudinal-template space.",
            "Template": "T1w-template",
        },
        "from-longitudinal_to-template_mode-image_desc-nonlinear_xfm": {
            "Description": "Nonlinear (warp field) transform from "
                           "longitudinal-template space to T1w-template "
                           "space.",
            "Template": "T1w-template",
        },
        "from-template_to-longitudinal_mode-image_desc-nonlinear_xfm": {
            "Description": "Nonlinear (warp field) transform from "
                           "T1w-template space to longitudinal-template "
                           "space.",
            "Template": "T1w-template",
        },
        "from-longitudinal_to-template_mode-image_xfm": {
            "Description": "Composite (affine + warp field) transform from "
                           "longitudinal-template space to T1w-template "
                           "space.",
            "Template": "T1w-template",
        },
        "from-template_to-longitudinal_mode-image_xfm": {
            "Description": "Composite (affine + warp field) transform from "
                           "T1w-template space to longitudinal-template "
                           "space.",
            "Template": "T1w-template",
        },
    },
)
def register_ANTs_anat_to_template(wf, cfg, strat_pool, pipe_num, opt=None):

    params = cfg.registration_workflows['anatomical_registration'][
        'registration']['ANTs']['T1_registration']

    ants_rc, outputs = ANTs_registration_connector('ANTS_T1_to_template_'
                                                   f'{pipe_num}', cfg,
                                                   params, orig='T1w')

    ants_rc.inputs.inputspec.interpolation = cfg.registration_workflows[
        'anatomical_registration']['registration']['ANTs']['interpolation']

    connect, brain = \
        strat_pool.get_data(['desc-preproc_T1w',
                             'space-longitudinal_desc-brain_T1w'],
                            report_fetched=True)
    node, out = connect
    wf.connect(node, out, ants_rc, 'inputspec.input_brain')

    t1w_brain_template = strat_pool.node_data('T1w-brain-template')
    wf.connect(t1w_brain_template.node, t1w_brain_template.out,
               ants_rc, 'inputspec.reference_brain')

    # TODO check the order of T1w
    node, out = strat_pool.get_data(["desc-restore_T1w", "desc-head_T1w",
                                     "desc-preproc_T1w",
                                     "space-longitudinal_desc-reorient_T1w"])
    wf.connect(node, out, ants_rc, 'inputspec.input_head')


    t1w_template = strat_pool.node_data('T1w-template')
    wf.connect(t1w_template.node, t1w_template.out,
               ants_rc, 'inputspec.reference_head')

    brain_mask = strat_pool.node_data(["space-T1w_desc-brain_mask",
                                       "space-longitudinal_desc-brain_mask",
                                       "space-T1w_desc-acpcbrain_mask"])
    wf.connect(brain_mask.node, brain_mask.out,
               ants_rc, 'inputspec.input_mask')

    if strat_pool.check_rpool('T1w-brain-template-mask'):
        node, out = strat_pool.get_data('T1w-brain-template-mask')
        wf.connect(node, out, ants_rc, 'inputspec.reference_mask')

    if strat_pool.check_rpool('label-lesion_mask'):
        node, out = strat_pool.get_data('label-lesion_mask')
        wf.connect(node, out, ants_rc, 'inputspec.lesion_mask')

    if 'space-longitudinal' in brain:
        for key in outputs:
            for direction in ['from', 'to']:
                if f'{direction}-T1w' in key:
                    new_key = key.replace(f'{direction}-T1w',
                                          f'{direction}-longitudinal')
                    outputs[new_key] = outputs[key]
                    del outputs[key]

    return (wf, outputs)


@nodeblock(
    name="register_symmetric_ANTs_anat_to_template",
    config=["registration_workflows", "anatomical_registration"],
    switch=["run"],
    option_key=["registration", "using"],
    option_val="ANTS",
    inputs=[
        (
            ["desc-preproc_T1w", "space-longitudinal_desc-brain_T1w"],
            ["space-T1w_desc-brain_mask",
             "space-longitudinal_desc-brain_mask"],
            [
                "desc-head_T1w",
                "desc-preproc_T1w",
                "space-longitudinal_desc-reorient_T1w",
            ],
        ),
        "T1w-template-symmetric",
        "T1w-brain-template-symmetric",
        "dilated-symmetric-brain-mask",
        "label-lesion_mask",
    ],
    outputs={
        "space-symtemplate_desc-preproc_T1w": {
            "Template": "T1w-brain-template-symmetric"
        },
        "from-T1w_to-symtemplate_mode-image_desc-linear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-symtemplate_to-T1w_mode-image_desc-linear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-T1w_to-symtemplate_mode-image_desc-nonlinear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-symtemplate_to-T1w_mode-image_desc-nonlinear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-T1w_to-symtemplate_mode-image_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-symtemplate_to-T1w_mode-image_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-longitudinal_to-symtemplate_mode-image_desc-linear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-symtemplate_to-longitudinal_mode-image_desc-linear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-longitudinal_to-symtemplate_mode-image_desc-nonlinear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-symtemplate_to-longitudinal_mode-image_desc-nonlinear_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-longitudinal_to-symtemplate_mode-image_xfm": {
            "Template": "T1w-template-symmetric"
        },
        "from-symtemplate_to-longitudinal_mode-image_xfm": {
            "Template": "T1w-template-symmetric"
        },
    },
)
def register_symmetric_ANTs_anat_to_template(wf, cfg, strat_pool, pipe_num,
                                             opt=None):

    params = cfg.registration_workflows['anatomical_registration'][
        'registration']['ANTs']['T1_registration']

    ants, outputs = ANTs_registration_connector('ANTS_T1_to_template_'
                                                f'symmetric_{pipe_num}', cfg,
                                                params, orig='T1w',
                                                symmetric=True)

    ants.inputs.inputspec.interpolation = cfg.registration_workflows[
        'anatomical_registration']['registration']['ANTs']['interpolation']

    connect, brain = \
        strat_pool.get_data(['desc-preproc_T1w',
                             'space-longitudinal_desc-brain_T1w'],
                            report_fetched=True)
    node, out = connect
    wf.connect(node, out, ants, 'inputspec.input_brain')

    node, out = strat_pool.get_data('T1w-brain-template-symmetric')
    wf.connect(node, out, ants, 'inputspec.reference_brain')

    node, out = strat_pool.get_data(["desc-head_T1w", "desc-preproc_T1w",
                                     "space-longitudinal_desc-reorient_T1w"])
    wf.connect(node, out, ants, 'inputspec.input_head')

    node, out = strat_pool.get_data('T1w-template-symmetric')
    wf.connect(node, out, ants, 'inputspec.reference_head')

    node, out = strat_pool.get_data(["space-T1w_desc-brain_mask",
                                     "space-longitudinal_desc-brain_mask"])
    wf.connect(node, out, ants, 'inputspec.input_mask')

    node, out = strat_pool.get_data('dilated-symmetric-brain-mask')
    wf.connect(node, out, ants, 'inputspec.reference_mask')

    if strat_pool.check_rpool('label-lesion_mask'):
        node, out = strat_pool.get_data('label-lesion_mask')
        wf.connect(node, out, ants, 'inputspec.lesion_mask')

    if 'space-longitudinal' in brain:
        for key in outputs.keys():
            if 'from-T1w' in key:
                new_key = key.replace('from-T1w', 'from-longitudinal')
                outputs[new_key] = outputs[key]
                del outputs[key]
            if 'to-T1w' in key:
                new_key = key.replace('to-T1w', 'to-longitudinal')
                outputs[new_key] = outputs[key]
                del outputs[key]

    return (wf, outputs)


@nodeblock(
    name="register_ANTs_EPI_to_template",
    config=["registration_workflows", "functional_registration", "EPI_registration"],
    switch=["run"],
    option_key="using",
    option_val="ANTS",
    inputs=[
        ("sbref", "space-bold_desc-brain_mask"),
        "EPI-template",
        "EPI-template-mask",
    ],
    outputs={
        "space-template_desc-preproc_bold": {"Template": "EPI-template"},
        "from-bold_to-EPItemplate_mode-image_desc-linear_xfm": {
            "Template": "EPI-template"
        },
        "from-EPItemplate_to-bold_mode-image_desc-linear_xfm": {
            "Template": "EPI-template"
        },
        "from-bold_to-EPItemplate_mode-image_desc-nonlinear_xfm": {
            "Template": "EPI-template"
        },
        "from-EPItemplate_to-bold_mode-image_desc-nonlinear_xfm": {
            "Template": "EPI-template"
        },
        "from-bold_to-EPItemplate_mode-image_xfm": {
            "Template": "EPI-template"},
        "from-EPItemplate_to-bold_mode-image_xfm": {
            "Template": "EPI-template"},
    },
)
def register_ANTs_EPI_to_template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Directly register the mean functional to an EPI template. No T1w
    involved.
    '''
    params = cfg.registration_workflows['functional_registration'][
        'EPI_registration']['ANTs']['parameters']

    ants, outputs = ANTs_registration_connector('ANTS_bold_to_EPI-template'
                                                f'_{pipe_num}', cfg, params,
                                                orig='bold', template='EPI')

    ants.inputs.inputspec.interpolation = cfg.registration_workflows[
        'functional_registration']['EPI_registration']['ANTs'][
        'interpolation']

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, ants, 'inputspec.input_brain')

    node, out = strat_pool.get_data('EPI-template')
    wf.connect(node, out, ants, 'inputspec.reference_brain')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, ants, 'inputspec.input_head')

    node, out = strat_pool.get_data('EPI-template')
    wf.connect(node, out, ants, 'inputspec.reference_head')

    node, out = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out, ants, 'inputspec.input_mask')

    if strat_pool.check_rpool('EPI-template-mask'):
        node, out = strat_pool.get_data('EPI-template-mask')
        wf.connect(node, out, ants, 'inputspec.reference_mask')

    return (wf, outputs)


@nodeblock(
    name="overwrite_transform_anat_to_template",
    switch=[
        ["registration_workflows", "anatomical_registration", "run"],
        [
            "registration_workflows",
            "anatomical_registration",
            "overwrite_transform",
            "run",
        ],
    ],
    option_key=[
        "registration_workflows",
        "anatomical_registration",
        "overwrite_transform",
        "using",
    ],
    option_val="FSL",
    inputs=[
        (
            "desc-restore-brain_T1w",
            ["desc-preproc_T1w", "space-longitudinal_desc-brain_T1w"],
            ["desc-restore_T1w", "desc-preproc_T1w", "desc-reorient_T1w",
             "T1w"],
            ["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
            "space-T1w_desc-brain_mask",
            "T1w-template",
            "from-T1w_to-template_mode-image_xfm",
            "from-template_to-T1w_mode-image_xfm",
            "space-template_desc-brain_T1w",
            "space-template_desc-preproc_T1w",
        )
    ],
    outputs={
        "space-template_desc-preproc_T1w": {"Template": "T1w-template"},
        "space-template_desc-head_T1w": {"Template": "T1w-template"},
        "space-template_desc-T1w_mask": {"Template": "T1w-template"},
        "from-T1w_to-template_mode-image_xfm": {"Template": "T1w-template"},
        "from-template_to-T1w_mode-image_xfm": {"Template": "T1w-template"},
    },
)
def overwrite_transform_anat_to_template(wf, cfg, strat_pool, pipe_num,
                                         opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-T1w_to-template_mode-image_xfm')

    reg_tool = check_prov_for_regtool(xfm_prov)

    if opt.lower() == 'fsl' and reg_tool.lower() == 'ants':

        # Apply head-to-head transforms on brain using ABCD-style registration
        # Convert ANTs warps to FSL warps to be consistent with the functional registration
        # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/scripts/AtlasRegistrationToMNI152_ANTsbased.sh#L134-L172

        # antsApplyTransforms -d 3 -i ${T1wRestore}.nii.gz -r ${Reference} \
        # -t ${WD}/xfms/T1w_to_MNI_3Warp.nii.gz \
        # -t ${WD}/xfms/T1w_to_MNI_2Affine.mat \
        # -t ${WD}/xfms/T1w_to_MNI_1Rigid.mat \
        # -t ${WD}/xfms/T1w_to_MNI_0DerivedInitialMovingTranslation.mat \
        # -o [${WD}/xfms/ANTs_CombinedWarp.nii.gz,1]
        ants_apply_warp_t1_to_template = pe.Node(interface=ants.ApplyTransforms(),
                                                name=f'ANTS-ABCD_T1_to_template_{pipe_num}')
        ants_apply_warp_t1_to_template.inputs.dimension = 3
        ants_apply_warp_t1_to_template.inputs.print_out_composite_warp_file = True
        ants_apply_warp_t1_to_template.inputs.output_image = 'ANTs_CombinedWarp.nii.gz'

        node, out = strat_pool.get_data(['desc-restore_T1w', 'desc-preproc_T1w'])
        wf.connect(node, out, ants_apply_warp_t1_to_template, 'input_image')

        node, out = strat_pool.get_data('T1w-template')
        wf.connect(node, out, ants_apply_warp_t1_to_template, 'reference_image')

        node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
        wf.connect(node, out, ants_apply_warp_t1_to_template, 'transforms')

        # antsApplyTransforms -d 3 -i ${T1wImage}.nii.gz -r ${Reference} \
        # -t [${WD}/xfms/T1w_to_MNI_0DerivedInitialMovingTranslation.mat,1] \
        # -t [${WD}/xfms/T1w_to_MNI_1Rigid.mat,1] \
        # -t [${WD}/xfms/T1w_to_MNI_2Affine.mat,1] \
        # -t ${WD}/xfms/T1w_to_MNI_3InverseWarp.nii.gz \
        # -o [${WD}/xfms/ANTs_CombinedInvWarp.nii.gz,1]

        # T1wImage is ACPC aligned head
        ants_apply_warp_template_to_t1 = pe.Node(interface=ants.ApplyTransforms(),
                                                name=f'ANTS-ABCD_template_to_T1_{pipe_num}')
        ants_apply_warp_template_to_t1.inputs.dimension = 3
        ants_apply_warp_template_to_t1.inputs.print_out_composite_warp_file = True
        ants_apply_warp_template_to_t1.inputs.output_image = 'ANTs_CombinedInvWarp.nii.gz'

        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, ants_apply_warp_template_to_t1, 'input_image')

        node, out = strat_pool.get_data('T1w-template')
        wf.connect(node, out, ants_apply_warp_template_to_t1, 'reference_image')

        node, out = strat_pool.get_data('from-template_to-T1w_mode-image_xfm')
        wf.connect(node, out, ants_apply_warp_template_to_t1, 'transforms')

        # c4d -mcs ${WD}/xfms/ANTs_CombinedWarp.nii.gz -oo ${WD}/xfms/e1.nii.gz ${WD}/xfms/e2.nii.gz ${WD}/xfms/e3.nii.gz
        # -mcs: -multicomponent-split, -oo: -output-multiple
        split_combined_warp = pe.Node(util.Function(input_names=['input',
                                                                'output_name'],
                                                    output_names=['output1',
                                                                'output2',
                                                                'output3'],
                                                    function=run_c4d),
                                    name=f'split_combined_warp_{pipe_num}')
        split_combined_warp.inputs.output_name = 'e'

        wf.connect(ants_apply_warp_t1_to_template, 'output_image',
            split_combined_warp, 'input')

        # c4d -mcs ${WD}/xfms/ANTs_CombinedInvWarp.nii.gz -oo ${WD}/xfms/e1inv.nii.gz ${WD}/xfms/e2inv.nii.gz ${WD}/xfms/e3inv.nii.gz
        split_combined_inv_warp = pe.Node(util.Function(input_names=['input',
                                                                    'output_name'],
                                                        output_names=['output1',
                                                                    'output2',
                                                                    'output3'],
                                                        function=run_c4d),
                                    name=f'split_combined_inv_warp_{pipe_num}')
        split_combined_inv_warp.inputs.output_name = 'einv'

        wf.connect(ants_apply_warp_template_to_t1, 'output_image',
            split_combined_inv_warp, 'input')

        # fslmaths ${WD}/xfms/e2.nii.gz -mul -1 ${WD}/xfms/e-2.nii.gz
        change_e2_sign = pe.Node(interface=fsl.maths.MathsCommand(),
                            name=f'change_e2_sign_{pipe_num}')
        change_e2_sign.inputs.args = '-mul -1'

        wf.connect(split_combined_warp, 'output2',
            change_e2_sign, 'in_file')

        # fslmaths ${WD}/xfms/e2inv.nii.gz -mul -1 ${WD}/xfms/e-2inv.nii.gz
        change_e2inv_sign = pe.Node(interface=fsl.maths.MathsCommand(),
                                name=f'change_e2inv_sign_{pipe_num}')
        change_e2inv_sign.inputs.args = '-mul -1'

        wf.connect(split_combined_inv_warp, 'output2',
            change_e2inv_sign, 'in_file')

        # fslmerge -t ${OutputTransform} ${WD}/xfms/e1.nii.gz ${WD}/xfms/e-2.nii.gz ${WD}/xfms/e3.nii.gz
        merge_xfms_to_list = pe.Node(util.Merge(3),
                                     name=f'merge_t1_to_template_xfms_to_list_{pipe_num}')

        wf.connect(split_combined_warp, 'output1',
            merge_xfms_to_list, 'in1')
        wf.connect(change_e2_sign, 'out_file',
            merge_xfms_to_list, 'in2')
        wf.connect(split_combined_warp, 'output3',
            merge_xfms_to_list, 'in3')

        merge_xfms = pe.Node(interface=fslMerge(),
                             name=f'merge_t1_to_template_xfms_{pipe_num}')
        merge_xfms.inputs.dimension = 't'

        wf.connect(merge_xfms_to_list, 'out',
            merge_xfms, 'in_files')

        # fslmerge -t ${OutputInvTransform} ${WD}/xfms/e1inv.nii.gz ${WD}/xfms/e-2inv.nii.gz ${WD}/xfms/e3inv.nii.gz
        merge_inv_xfms_to_list = pe.Node(util.Merge(3),
                                         name=f'merge_template_to_t1_xfms_to_list_{pipe_num}')

        wf.connect(split_combined_inv_warp, 'output1',
            merge_inv_xfms_to_list, 'in1')
        wf.connect(change_e2inv_sign, 'out_file',
            merge_inv_xfms_to_list, 'in2')
        wf.connect(split_combined_inv_warp, 'output3',
            merge_inv_xfms_to_list, 'in3')

        merge_inv_xfms = pe.Node(interface=fslMerge(),
                                 name=f'merge_template_to_t1_xfms_{pipe_num}')
        merge_inv_xfms.inputs.dimension = 't'

        wf.connect(merge_inv_xfms_to_list, 'out',
            merge_inv_xfms, 'in_files')

        # applywarp --rel --interp=spline -i ${T1wRestore} -r ${Reference} -w ${OutputTransform} -o ${OutputT1wImageRestore}
        fsl_apply_warp_t1_to_template = pe.Node(interface=fsl.ApplyWarp(),
                                                name=f'FSL-ABCD_T1_to_template_{pipe_num}')
        fsl_apply_warp_t1_to_template.inputs.relwarp = True
        fsl_apply_warp_t1_to_template.inputs.interp = 'spline'

        node, out = strat_pool.get_data(['desc-restore_T1w', 'desc-preproc_T1w'])
        wf.connect(node, out, fsl_apply_warp_t1_to_template, 'in_file')

        node, out = strat_pool.get_data('T1w-template')
        wf.connect(node, out, fsl_apply_warp_t1_to_template, 'ref_file')

        wf.connect(merge_xfms, 'merged_file',
            fsl_apply_warp_t1_to_template, 'field_file')

        # applywarp --rel --interp=nn -i ${T1wRestoreBrain} -r ${Reference} -w ${OutputTransform} -o ${OutputT1wImageRestoreBrain}
        fsl_apply_warp_t1_brain_to_template = pe.Node(interface=fsl.ApplyWarp(),
                                                      name=f'FSL-ABCD_T1_brain_to_template_{pipe_num}')
        fsl_apply_warp_t1_brain_to_template.inputs.relwarp = True
        fsl_apply_warp_t1_brain_to_template.inputs.interp = 'nn'

        # TODO connect T1wRestoreBrain, check T1wRestoreBrain quality
        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, fsl_apply_warp_t1_brain_to_template, 'in_file')

        node, out = strat_pool.get_data('T1w-template')
        wf.connect(node, out, fsl_apply_warp_t1_brain_to_template, 'ref_file')

        wf.connect(merge_xfms, 'merged_file',
            fsl_apply_warp_t1_brain_to_template, 'field_file')

        fsl_apply_warp_t1_brain_mask_to_template = pe.Node(interface=fsl.ApplyWarp(),
                                                        name=f'FSL-ABCD_T1_brain_mask_to_template_{pipe_num}')
        fsl_apply_warp_t1_brain_mask_to_template.inputs.relwarp = True
        fsl_apply_warp_t1_brain_mask_to_template.inputs.interp = 'nn'

        node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
        wf.connect(node, out, fsl_apply_warp_t1_brain_mask_to_template, 'in_file')

        node, out = strat_pool.get_data('T1w-template')
        wf.connect(node, out, fsl_apply_warp_t1_brain_mask_to_template, 'ref_file')

        wf.connect(merge_xfms, 'merged_file',
            fsl_apply_warp_t1_brain_mask_to_template, 'field_file')

        # fslmaths ${OutputT1wImageRestore} -mas ${OutputT1wImageRestoreBrain} ${OutputT1wImageRestoreBrain}
        apply_mask = pe.Node(interface=fsl.maths.ApplyMask(),
                            name=f'get_t1_brain_{pipe_num}')

        wf.connect(fsl_apply_warp_t1_to_template, 'out_file',
            apply_mask, 'in_file')

        wf.connect(fsl_apply_warp_t1_brain_to_template, 'out_file',
            apply_mask, 'mask_file')

        outputs = {
            'space-template_desc-preproc_T1w': (apply_mask, 'out_file'),
            'space-template_desc-head_T1w': (fsl_apply_warp_t1_to_template, 'out_file'),
            'space-template_desc-T1w_mask': (fsl_apply_warp_t1_brain_mask_to_template, 'out_file'),
            'from-T1w_to-template_mode-image_xfm': (merge_xfms, 'merged_file'),
            'from-template_to-T1w_mode-image_xfm': (merge_inv_xfms, 'merged_file')
        }

    return (wf, outputs)


@nodeblock(
    name="coregistration_prep_vol",
    switch=["functional_preproc", "run"],
    option_key=[
        "registration_workflows",
        "functional_registration",
        "coregistration",
        "func_input_prep",
        "input",
    ],
    option_val="Selected_Functional_Volume",
    inputs=[("desc-brain_bold", ["desc-motion_bold", "bold"], "sbref")],
    outputs=["sbref"],
)
def coregistration_prep_vol(wf, cfg, strat_pool, pipe_num, opt=None):

    get_func_volume = pe.Node(interface=afni.Calc(),
                              name=f'get_func_volume_{pipe_num}')

    get_func_volume.inputs.set(
        expr='a',
        single_idx=cfg.registration_workflows['functional_registration']['coregistration'][
            'func_input_prep']['Selected Functional Volume']['func_reg_input_volume'],
        outputtype='NIFTI_GZ'
    )

    if not cfg.registration_workflows['functional_registration'][
        'coregistration']['func_input_prep']['reg_with_skull']:
        node, out = strat_pool.get_data("desc-brain_bold")
    else:
        # TODO check which file is functional_skull_leaf
        # TODO add a function to choose brain or skull?
        node, out = strat_pool.get_data(["desc-motion_bold", "bold"])

    wf.connect(node, out, get_func_volume, 'in_file_a')

    coreg_input = (get_func_volume, 'out_file')

    outputs = {
        'sbref': coreg_input
    }

    return (wf, outputs)


@nodeblock(
    name="coregistration_prep_mean",
    switch=[["functional_preproc", "run"],
            ["functional_preproc", "coreg_prep", "run"]],
    option_key=[
        "registration_workflows",
        "functional_registration",
        "coregistration",
        "func_input_prep",
        "input",
    ],
    option_val="Mean_Functional",
    inputs=["desc-mean_bold"],
    outputs=["sbref"],
)
def coregistration_prep_mean(wf, cfg, strat_pool, pipe_num, opt=None):

    coreg_input = strat_pool.get_data("desc-mean_bold")

    # TODO add mean skull
    if cfg.registration_workflows['functional_registration'][
            'coregistration']['func_input_prep']['Mean Functional'][
            'n4_correct_func']:
        n4_correct_func = pe.Node(
            interface=
            ants.N4BiasFieldCorrection(dimension=3,
                                       copy_header=True,
                                       bspline_fitting_distance=200),
            shrink_factor=2,
            name=f'func_mean_n4_corrected_{pipe_num}')
        n4_correct_func.inputs.args = '-r True'

        node, out = coreg_input
        wf.connect(node, out, n4_correct_func, 'input_image')

        coreg_input = (n4_correct_func, 'output_image')

    outputs = {
        'sbref': coreg_input
    }

    return (wf, outputs)


@nodeblock(
    name="coregistration_prep_fmriprep",
    switch=[["functional_preproc", "run"],
            ["functional_preproc", "coreg_prep", "run"]],
    option_key=[
        "registration_workflows",
        "functional_registration",
        "coregistration",
        "func_input_prep",
        "input",
    ],
    option_val="fmriprep_reference",
    inputs=["desc-ref_bold"],
    outputs=["sbref"],
)
def coregistration_prep_fmriprep(wf, cfg, strat_pool, pipe_num, opt=None):

    coreg_input = strat_pool.get_data("desc-ref_bold")

    outputs = {
        'sbref': coreg_input
    }

    return (wf, outputs)


@nodeblock(
    name="coregistration",
    config=["registration_workflows", "functional_registration", "coregistration"],
    switch=["run"],
    inputs=[
        (
            "sbref",
            "desc-motion_bold",
            "space-bold_label-WM_mask",
            "despiked-fieldmap",
            "fieldmap-mask",
            "effectiveEchoSpacing",
            "pe-direction",
        ),
        (
            "desc-preproc_T1w",
            "desc-restore-brain_T1w",
            "desc-preproc_T2w",
            "desc-preproc_T2w",
            "T2w",
            ["label-WM_probseg", "label-WM_mask"],
            ["label-WM_pveseg", "label-WM_mask"],
            "desc-head_T1w",
            "desc-head_T2w",
        ),
    ],
    outputs=[
        "space-T1w_sbref",
        "from-bold_to-T1w_mode-image_desc-linear_xfm",
        "from-bold_to-T1w_mode-image_desc-linear_warp",
    ],
)
def coregistration(wf, cfg, strat_pool, pipe_num, opt=None):
    diff_complete = False
    if strat_pool.check_rpool("despiked-fieldmap") and \
            strat_pool.check_rpool("fieldmap-mask"):
        diff_complete = True

    if strat_pool.check_rpool('T2w') and cfg.anatomical_preproc['run_t2']:
        # monkey data
        func_to_anat = create_register_func_to_anat_use_T2(cfg,
                                                    f'func_to_anat_FLIRT_'
                                                    f'{pipe_num}')

        # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh#L177
        # fslmaths "$fMRIFolder"/"$NameOffMRI"_mc -Tmean "$fMRIFolder"/"$ScoutName"_gdc
        func_mc_mean = pe.Node(interface=afni_utils.TStat(),
                            name=f'func_motion_corrected_mean_{pipe_num}')

        func_mc_mean.inputs.options = '-mean'
        func_mc_mean.inputs.outputtype = 'NIFTI_GZ'

        node, out = strat_pool.get_data("desc-motion_bold")
        wf.connect(node, out, func_mc_mean, 'in_file')

        wf.connect(func_mc_mean, 'out_file', func_to_anat, 'inputspec.func')

        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, func_to_anat, 'inputspec.T1_brain')

        node, out = strat_pool.get_data('desc-head_T2w')
        wf.connect(node, out, func_to_anat, 'inputspec.T2_head')

        node, out = strat_pool.get_data('desc-preproc_T2w')
        wf.connect(node, out, func_to_anat, 'inputspec.T2_brain')

    else:
        # if field map-based distortion correction is on, but BBR is off,
        # send in the distortion correction files here
        func_to_anat = create_register_func_to_anat(cfg, diff_complete,
                                                    f'func_to_anat_FLIRT_'
                                                    f'{pipe_num}')

        func_to_anat.inputs.inputspec.dof = cfg.registration_workflows[
        'functional_registration']['coregistration']['dof']

        func_to_anat.inputs.inputspec.interp = cfg.registration_workflows[
        'functional_registration']['coregistration']['interpolation']

        node, out = strat_pool.get_data('sbref')
        wf.connect(node, out, func_to_anat, 'inputspec.func')

        if cfg.registration_workflows['functional_registration'][
            'coregistration']['reference'] == 'brain':
            # TODO: use JSON meta-data to confirm
            node, out = strat_pool.get_data('desc-preproc_T1w')
        elif cfg.registration_workflows['functional_registration'][
            'coregistration']['reference'] == 'restore-brain':
            node, out = strat_pool.get_data('desc-restore-brain_T1w')
        wf.connect(node, out, func_to_anat, 'inputspec.anat')

    if diff_complete:
        node, out = strat_pool.get_data('effectiveEchoSpacing')
        wf.connect(node, out, func_to_anat, 'echospacing_input.echospacing')

        node, out = strat_pool.get_data('pe-direction')
        wf.connect(node, out, func_to_anat, 'pedir_input.pedir')

        node, out = strat_pool.get_data("despiked-fieldmap")
        wf.connect(node, out, func_to_anat, 'inputspec.fieldmap')

        node, out = strat_pool.get_data("fieldmap-mask")
        wf.connect(node, out, func_to_anat, 'inputspec.fieldmapmask')

    if strat_pool.check_rpool('T2w') and cfg.anatomical_preproc['run_t2']:
        outputs = {
            'space-T1w_sbref':
                (func_to_anat, 'outputspec.anat_func_nobbreg'),
            'from-bold_to-T1w_mode-image_desc-linear_xfm':
                (func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg'),
            'from-bold_to-T1w_mode-image_desc-linear_warp':
                (func_to_anat, 'outputspec.func_to_anat_linear_warp_nobbreg')
        }
    else:
        outputs = {
            'space-T1w_sbref':
                (func_to_anat, 'outputspec.anat_func_nobbreg'),
            'from-bold_to-T1w_mode-image_desc-linear_xfm':
                (func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg')
        }

    if True in cfg.registration_workflows['functional_registration'][
        'coregistration']["boundary_based_registration"]["run"]:

        func_to_anat_bbreg = create_bbregister_func_to_anat(diff_complete,
                                                            f'func_to_anat_'
                                                            f'bbreg_'
                                                            f'{pipe_num}')
        func_to_anat_bbreg.inputs.inputspec.bbr_schedule = \
            cfg.registration_workflows['functional_registration'][
                'coregistration']['boundary_based_registration'][
                'bbr_schedule']

        func_to_anat_bbreg.inputs.inputspec.bbr_wm_mask_args = \
            cfg.registration_workflows['functional_registration'][
                'coregistration']['boundary_based_registration'][
                'bbr_wm_mask_args']

        node, out = strat_pool.get_data('sbref')
        wf.connect(node, out, func_to_anat_bbreg, 'inputspec.func')

        if cfg.registration_workflows['functional_registration'][
                'coregistration']['boundary_based_registration'][
                'reference'] == 'whole-head':
            node, out = strat_pool.get_data('desc-head_T1w')
            wf.connect(node, out, func_to_anat_bbreg, 'inputspec.anat')

        elif cfg.registration_workflows['functional_registration'][
                'coregistration']['boundary_based_registration'][
                'reference'] == 'brain':
            node, out = strat_pool.get_data('desc-preproc_T1w')
            wf.connect(node, out, func_to_anat_bbreg, 'inputspec.anat')

        wf.connect(func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg',
                   func_to_anat_bbreg, 'inputspec.linear_reg_matrix')

        if strat_pool.check_rpool('space-bold_label-WM_mask'):
            node, out = strat_pool.get_data(["space-bold_label-WM_mask"])
            wf.connect(node, out,
                       func_to_anat_bbreg, 'inputspec.anat_wm_segmentation')
        else:
            if cfg.registration_workflows['functional_registration'][
                'coregistration']['boundary_based_registration']['bbr_wm_map'] == 'probability_map':
                node, out = strat_pool.get_data(["label-WM_probseg",
                                                 "label-WM_mask"])
            elif cfg.registration_workflows['functional_registration'][
                'coregistration']['boundary_based_registration']['bbr_wm_map'] == 'partial_volume_map':
                node, out = strat_pool.get_data(["label-WM_pveseg",
                                                 "label-WM_mask"])
            wf.connect(node, out,
                       func_to_anat_bbreg, 'inputspec.anat_wm_segmentation')

        if diff_complete:
            node, out = strat_pool.get_data('effectiveEchoSpacing')
            wf.connect(node, out,
                       func_to_anat_bbreg, 'echospacing_input.echospacing')

            node, out = strat_pool.get_data('pe-direction')
            wf.connect(node, out, func_to_anat_bbreg, 'pedir_input.pedir')

            node, out = strat_pool.get_data("despiked-fieldmap")
            wf.connect(node, out, func_to_anat_bbreg, 'inputspec.fieldmap')

            node, out = strat_pool.get_data("fieldmap-mask")
            wf.connect(node, out,
                       func_to_anat_bbreg, 'inputspec.fieldmapmask')

        outputs = {
            'space-T1w_sbref':
                (func_to_anat_bbreg, 'outputspec.anat_func'),
            'from-bold_to-T1w_mode-image_desc-linear_xfm':
                (func_to_anat_bbreg, 'outputspec.func_to_anat_linear_xfm')
        }

    return (wf, outputs)


@nodeblock(
    name="create_func_to_T1template_xfm",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run"],
    option_key=["target_template", "using"],
    option_val="T1_template",
    inputs=[
        (
            "sbref",
            "from-bold_to-T1w_mode-image_desc-linear_xfm",
            "ants-blip-warp",
            "fsl-blip-warp",
        ),
        (
            "from-T1w_to-template_mode-image_xfm",
            "from-template_to-T1w_mode-image_xfm",
            "desc-brain_T1w",
        ),
        "T1w-brain-template-funcreg",
    ],
    outputs={
        "from-bold_to-template_mode-image_xfm": {
            "Template": "T1w-brain-template-funcreg"
        },
        "from-template_to-bold_mode-image_xfm": {
            "Template": "T1w-brain-template-funcreg"
        },
    },
)
def create_func_to_T1template_xfm(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Condense the BOLD-to-T1 coregistration transform and the T1-to-template
    transform into one transform matrix.
    '''
    xfm_prov = strat_pool.get_cpac_provenance(
        'from-T1w_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    xfm, outputs = bold_to_T1template_xfm_connector('create_func_to_T1w'
                                                    f'template_xfm_{pipe_num}',
                                                    cfg, reg_tool,
                                                    symmetric=False)

    node, out = strat_pool.get_data(
        'from-bold_to-T1w_mode-image_desc-linear_xfm')
    wf.connect(node, out, xfm, 'inputspec.coreg_xfm')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, xfm, 'inputspec.input_brain')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, xfm, 'inputspec.mean_bold')

    node, out = strat_pool.get_data('T1w-brain-template-funcreg')
    wf.connect(node, out, xfm, 'inputspec.T1w-brain-template_funcreg')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, xfm, 'inputspec.T1w_to_template_xfm')

    # FNIRT pipelines don't have an inverse nonlinear warp, make optional
    if strat_pool.check_rpool('from-template_to-T1w_mode-image_xfm'):
        node, out = strat_pool.get_data('from-template_to-T1w_mode-image_xfm')
        wf.connect(node, out, xfm, 'inputspec.template_to_T1w_xfm')

    if strat_pool.check_rpool('ants-blip-warp'):
        if reg_tool == 'ants':
            node, out = strat_pool.get_data('ants-blip-warp')
            wf.connect(node, out, xfm, 'inputspec.blip_warp')
        elif reg_tool == 'fsl':
            # apply the ants blip warp separately
            pass
    elif strat_pool.check_rpool('fsl-blip-warp'):
        if reg_tool == 'fsl':
            node, out = strat_pool.get_data('fsl-blip-warp')
            wf.connect(node, out, xfm, 'inputspec.blip_warp')
        elif reg_tool == 'ants':
            # apply the fsl blip warp separately
            pass

    return (wf, outputs)


@nodeblock(
    name="create_func_to_T1template_symmetric_xfm",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run"],
    option_key=["target_template", "using"],
    option_val="T1_template",
    inputs=[
        (
            "from-T1w_to-symtemplate_mode-image_xfm",
            "from-symtemplate_to-T1w_mode-image_xfm",
            "desc-brain_T1w",
        ),
        ("sbref", "from-bold_to-T1w_mode-image_desc-linear_xfm"),
        "T1w-brain-template-symmetric-deriv",
    ],
    outputs={
        "from-bold_to-symtemplate_mode-image_xfm": {
            "Template": "T1w-brain-template-symmetric-deriv"
        },
        "from-symtemplate_to-bold_mode-image_xfm": {
            "Template": "T1w-brain-template-symmetric-deriv"
        },
    },
)
def create_func_to_T1template_symmetric_xfm(wf, cfg, strat_pool, pipe_num,
                                            opt=None):
    '''Condense the BOLD-to-T1 coregistration transform and the T1-to-
    symmetric-template transform into one transform matrix.
    '''

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-T1w_to-symtemplate_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    xfm, outputs = bold_to_T1template_xfm_connector('create_func_to_T1wsymtem'
                                                    f'plate_xfm_{pipe_num}',
                                                    cfg, reg_tool,
                                                    symmetric=True)

    node, out = strat_pool.get_data(
        'from-bold_to-T1w_mode-image_desc-linear_xfm')
    wf.connect(node, out, xfm, 'inputspec.coreg_xfm')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, xfm, 'inputspec.input_brain')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, xfm, 'inputspec.mean_bold')

    node, out = strat_pool.get_data('T1w-brain-template-symmetric-deriv')
    wf.connect(node, out, xfm, 'inputspec.T1w-brain-template_funcreg')

    node, out = strat_pool.get_data('from-T1w_to-symtemplate_mode-image_xfm')
    wf.connect(node, out, xfm, 'inputspec.T1w_to_template_xfm')

    # FNIRT pipelines don't have an inverse nonlinear warp, make optional
    if strat_pool.check_rpool('from-symtemplate_to-T1w_mode-image_xfm'):
        node, out = \
            strat_pool.get_data('from-symtemplate_to-T1w_mode-image_xfm')
        wf.connect(node, out, xfm, 'inputspec.template_to_T1w_xfm')

    return (wf, outputs)


@nodeblock(
    name="apply_phasediff_to_timeseries_separately",
    switch=[
        [
            "registration_workflows",
            "functional_registration",
            "func_registration_to_template",
            "run",
        ],
        ["functional_preproc", "distortion_correction", "run"],
    ],
    option_key=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
        "apply_transform",
        "using",
    ],
    option_val=["default", "single_step_resampling_from_stc", "abcd"],
    inputs=[
        (
            "sbref",
            "desc-preproc_bold",
            "desc-stc_bold",
            "bold",
            "from-bold_to-T1w_mode-image_desc-linear_xfm",
        ),
        "despiked-fieldmap",
        "pe-direction",
        "effectiveEchoSpacing",
    ],
    outputs=["sbref", "desc-preproc_bold", "desc-stc_bold", "bold"],
)
def apply_phasediff_to_timeseries_separately(wf, cfg, strat_pool, pipe_num,
                                             opt=None):

    outputs = {'desc-preproc_bold': strat_pool.get_data("desc-preproc_bold")}
    if not strat_pool.check_rpool("despiked-fieldmap"):
        return (wf, outputs)

    invert_coreg_xfm = pe.Node(interface=fsl.ConvertXFM(),
        name=f'invert_coreg_xfm_{pipe_num}')
    invert_coreg_xfm.inputs.invert_xfm = True

    node, out = strat_pool.get_data("from-bold_to-T1w_mode-image_desc-linear_xfm")
    wf.connect(node, out, invert_coreg_xfm, 'in_file')

    warp_fmap = pe.Node(interface=fsl.ApplyWarp(),
        name=f'warp_fmap_{pipe_num}')

    node, out = strat_pool.get_data('despiked-fieldmap')
    wf.connect(node, out, warp_fmap, 'in_file')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, warp_fmap, 'ref_file')

    wf.connect(invert_coreg_xfm, 'out_file', warp_fmap, 'premat')

    mask_fmap = pe.Node(interface=fsl.maths.MathsCommand(),
                        name=f'mask_fmap_{pipe_num}')
    mask_fmap.inputs.args = '-abs -bin'

    wf.connect(warp_fmap, 'out_file', mask_fmap, 'in_file')

    conv_pedir = \
        pe.Node(interface=util.Function(input_names=['pedir',
                                                     'convert'],
                                        output_names=['pedir'],
                                        function=convert_pedir),
                name=f'apply_phasediff_convert_pedir_{pipe_num}')
    conv_pedir.inputs.convert = 'ijk_to_xyz'

    node, out = strat_pool.get_data('pe-direction')
    wf.connect(node, out, conv_pedir, 'pedir')

    fugue_saveshift = pe.Node(interface=fsl.FUGUE(),
                              name=f'fugue_saveshift_{pipe_num}')
    fugue_saveshift.inputs.save_shift = True

    wf.connect(warp_fmap, 'out_file', fugue_saveshift, 'fmap_in_file')
    wf.connect(mask_fmap, 'out_file', fugue_saveshift, 'mask_file')

    # FSL calls effective echo spacing = dwell time (not accurate)
    node, out = strat_pool.get_data('effectiveEchoSpacing')
    wf.connect(node, out, fugue_saveshift, 'dwell_time')

    wf.connect(conv_pedir, 'pedir', fugue_saveshift, 'unwarp_direction')

    shift_warp = pe.Node(interface=fsl.ConvertWarp(), 
                         name=f'shift_warp_{pipe_num}')
    shift_warp.inputs.out_relwarp = True

    wf.connect(fugue_saveshift, 'shift_out_file', shift_warp, 'shift_in_file')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, shift_warp, 'reference')

    wf.connect(conv_pedir, 'pedir', shift_warp, 'shift_direction')

    warp_bold = pe.Node(interface=fsl.ApplyWarp(),
                        name=f'warp_bold_phasediff_{pipe_num}')
    warp_bold.inputs.relwarp = True
    warp_bold.inputs.interp = 'spline'

    if opt == 'default':
        node, out = strat_pool.get_data('desc-preproc_bold')
        out_label = 'desc-preproc_bold'
    elif opt == 'single_step_resampling_from_stc':
        node, out = strat_pool.get_data('desc-stc_bold')
        out_label = 'desc-stc_bold'
    elif opt == 'abcd':
        node, out = strat_pool.get_data('bold')
        out_label = 'bold'

    wf.connect(node, out, warp_bold, 'in_file')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, warp_bold, 'ref_file')

    wf.connect(shift_warp, 'out_file', warp_bold, 'field_file')
    
    warp_sbref = pe.Node(interface=fsl.ApplyWarp(),
                        name=f'warp_sbref_phasediff_{pipe_num}')
    warp_sbref.inputs.relwarp = True
    warp_sbref.inputs.interp = 'spline'

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, warp_sbref, 'in_file')
    wf.connect(node, out, warp_sbref, 'ref_file')

    wf.connect(shift_warp, 'out_file', warp_sbref, 'field_file')

    outputs = {
        out_label: (warp_bold, 'out_file'),
        'sbref': (warp_sbref, 'out_file')
    }

    return (wf, outputs)


@nodeblock(
    name="apply_blip_to_timeseries_separately",
    switch=[
        [
            "registration_workflows",
            "functional_registration",
            "func_registration_to_template",
            "run",
        ],
        ["functional_preproc", "distortion_correction", "run"],
    ],
    option_key=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
        "apply_transform",
        "using",
    ],
    option_val=["default", "single_step_resampling_from_stc", "abcd"],
    inputs=[
        (
            "sbref",
            "desc-preproc_bold",
            "desc-stc_bold",
            "bold",
            "from-bold_to-template_mode-image_xfm",
            "ants-blip-warp",
            "fsl-blip-warp",
        )
    ],
    outputs=["desc-preproc_bold", "desc-stc_bold", "bold"],
)
def apply_blip_to_timeseries_separately(wf, cfg, strat_pool, pipe_num,
                                        opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-bold_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    outputs = {'desc-preproc_bold': strat_pool.get_data("desc-preproc_bold")}
    if strat_pool.check_rpool("ants-blip-warp"):
        if reg_tool == 'fsl':
            blip_node, blip_out = strat_pool.get_data("ants-blip-warp")
            reg_tool = 'ants'
        else:
            return (wf, outputs)
    elif strat_pool.check_rpool("fsl-blip-warp"):
        if reg_tool == 'ants':
            blip_node, blip_out = strat_pool.get_data("fsl-blip-warp")
            reg_tool = 'fsl'
        else:
            return (wf, outputs)
    else:
        return (wf, outputs)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    apply_xfm = apply_transform(f'warp_ts_to_blip_sep_{pipe_num}', reg_tool,
                                time_series=True, num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)

    if reg_tool == 'ants':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'ANTs_pipelines']['interpolation']
    elif reg_tool == 'fsl':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'FNIRT_pipelines']['interpolation']

    connect = strat_pool.get_data("desc-preproc_bold")

    if opt == 'default':
        node, out = strat_pool.get_data('desc-preproc_bold')
        out_label = 'desc-preproc_bold'
    elif opt == 'single_step_resampling_from_stc':
        node, out = strat_pool.get_data('desc-stc_bold')
        out_label = 'desc-stc_bold'
    elif opt == 'abcd':
        node, out = strat_pool.get_data('bold')
        out_label = 'bold'

    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("sbref")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    wf.connect(blip_node, blip_out, apply_xfm, 'inputspec.transform')

    outputs = {
        out_label: (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="transform_whole_head_T1w_to_T1template",
    config=["registration_workflows", "anatomical_registration"],
    switch=["run"],
    inputs=[
        (
            "desc-head_T1w",
            "from-T1w_to-template_mode-image_xfm",
            "space-template_desc-head_T1w",
        ),
        "T1w-template",
    ],
    outputs={"space-template_desc-head_T1w": {"Template": "T1w-template"}},
)
def warp_wholeheadT1_to_template(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-T1w_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    apply_xfm = apply_transform(f'warp_wholehead_T1w_to_T1template_{pipe_num}', 
                                reg_tool, time_series=False, num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)

    if reg_tool == 'ants':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'ANTs_pipelines']['interpolation']
    elif reg_tool == 'fsl':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'FNIRT_pipelines']['interpolation']

    connect = strat_pool.get_data("desc-head_T1w")
    node, out = connect
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("T1w-template")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-T1w_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        'space-template_desc-head_T1w': (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="transform_T1mask_to_T1template",
    switch=[
        ["registration_workflows", "anatomical_registration", "run"],
        ["anatomical_preproc", "run"],
        ["anatomical_preproc", "brain_extraction", "run"],
    ],
    inputs=[
        ("space-T1w_desc-brain_mask", "from-T1w_to-template_mode-image_xfm"),
        "T1w-template",
    ],
    outputs={"space-template_desc-brain_mask": {"Template": "T1w-template"}},
)
def warp_T1mask_to_template(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-T1w_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    apply_xfm = apply_transform(f'warp_T1mask_to_T1template_{pipe_num}', 
                                reg_tool, time_series=False, num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)

    apply_xfm.inputs.inputspec.interpolation = "NearestNeighbor"
    '''
    if reg_tool == 'ants':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'ANTs_pipelines']['interpolation']
    elif reg_tool == 'fsl':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'FNIRT_pipelines']['interpolation']
    '''
    connect = strat_pool.get_data("space-T1w_desc-brain_mask")
    node, out = connect
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("T1w-template")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-T1w_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        'space-template_desc-brain_mask': (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="transform_timeseries_to_T1template",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run"],
    option_key=["apply_transform", "using"],
    option_val="default",
    inputs=[
        ("desc-preproc_bold", "from-bold_to-template_mode-image_xfm"),
        "T1w-brain-template-funcreg",
    ],
    outputs={
        "space-template_desc-preproc_bold": {
            "Template": "T1w-brain-template-funcreg"}
    },
)
def warp_timeseries_to_T1template(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-bold_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    apply_xfm = apply_transform(f'warp_ts_to_T1template_{pipe_num}', reg_tool,
                                time_series=True, num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)

    if reg_tool == 'ants':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'ANTs_pipelines']['interpolation']
    elif reg_tool == 'fsl':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'FNIRT_pipelines']['interpolation']

    connect = strat_pool.get_data("desc-preproc_bold")
    node, out = connect
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("T1w-brain-template-funcreg")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-bold_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        'space-template_desc-preproc_bold': (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="transform_timeseries_to_T1template_deriv",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run"],
    option_key=["apply_transform", "using"],
    option_val="default",
    inputs=[
        ("desc-preproc_bold", "from-bold_to-template_mode-image_xfm"),
        "T1w-brain-template-funcreg",
    ],
    outputs={
        "space-template_res-derivative_desc-preproc_bold": {
            "Template": "T1w-brain-template-deriv"
        }
    },
)
def warp_timeseries_to_T1template_deriv(wf, cfg, strat_pool, pipe_num,
                                        opt=None):
    xfm_prov = strat_pool.get_cpac_provenance(
        'from-bold_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    apply_xfm = apply_transform(f'warp_ts_to_T1template_{pipe_num}', reg_tool,
                                time_series=True, num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)

    if reg_tool == 'ants':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'ANTs_pipelines']['interpolation']
    elif reg_tool == 'fsl':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'FNIRT_pipelines']['interpolation']

    connect = strat_pool.get_data("desc-preproc_bold")
    node, out = connect
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("T1w-brain-template-deriv")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-bold_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        'space-template_res-derivative_desc-preproc_bold': 
            (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="transform_timeseries_to_T1template_abcd",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run"],
    option_key=["apply_transform", "using"],
    option_val="abcd",
    inputs=[
        ("desc-preproc_bold", "bold", "motion-basefile",
         "coordinate-transformation"),
        "from-T1w_to-template_mode-image_xfm",
        "from-bold_to-T1w_mode-image_desc-linear_xfm",
        "from-bold_to-template_mode-image_xfm",
        "fsl-blip-warp",
        "desc-preproc_T1w",
        "space-template_res-bold_desc-brain_T1w",
        "space-template_desc-bold_mask",
        "T1w-brain-template-funcreg",
    ],
    outputs={
        "space-template_desc-preproc_bold": {
            "Template": "T1w-brain-template-funcreg"},
        "space-template_desc-scout_bold": {
            "Template": "T1w-brain-template-funcreg"},
        "space-template_desc-head_bold": {
            "Template": "T1w-brain-template-funcreg"},
    },
)
def warp_timeseries_to_T1template_abcd(wf, cfg, strat_pool, pipe_num, opt=None
                                       ):
    # Apply motion correction, coreg, anat-to-template transforms on raw functional timeseries using ABCD-style registration
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/scripts/OneStepResampling.sh#L168-L197

    # https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/scripts/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased.sh#L548
    # convertwarp --relout --rel -m ${WD}/fMRI2str.mat --ref=${T1wImage} --out=${WD}/fMRI2str.nii.gz
    convert_func_to_anat_linear_warp = pe.Node(interface=fsl.ConvertWarp(),
        name=f'convert_func_to_anat_linear_warp_{pipe_num}')

    convert_func_to_anat_linear_warp.inputs.out_relwarp = True
    convert_func_to_anat_linear_warp.inputs.relwarp = True
    
    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, convert_func_to_anat_linear_warp, 'reference')
    
    if strat_pool.check_rpool('fsl-blip-warp'):
        node, out = strat_pool.get_data('from-bold_to-T1w_mode-image_desc-linear_xfm')
        wf.connect(node, out, convert_func_to_anat_linear_warp, 'postmat')

        node, out = strat_pool.get_data('fsl-blip-warp')
        wf.connect(node, out, convert_func_to_anat_linear_warp, 'warp1')
    else:
        node, out = strat_pool.get_data('from-bold_to-T1w_mode-image_desc-linear_xfm')
        wf.connect(node, out, convert_func_to_anat_linear_warp, 'premat')

    # https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/scripts/OneStepResampling.sh#L140
    # convertwarp --relout --rel --warp1=${fMRIToStructuralInput} --warp2=${StructuralToStandard} --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${OutputTransform}
    convert_func_to_standard_warp = pe.Node(interface=fsl.ConvertWarp(),
        name=f'convert_func_to_standard_warp_{pipe_num}')

    convert_func_to_standard_warp.inputs.out_relwarp = True
    convert_func_to_standard_warp.inputs.relwarp = True

    wf.connect(convert_func_to_anat_linear_warp, 'out_file',
        convert_func_to_standard_warp, 'warp1')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, convert_func_to_standard_warp, 'warp2')

    node, out = strat_pool.get_data('space-template_res-bold_desc-brain_T1w')
    wf.connect(node, out, convert_func_to_standard_warp, 'reference')

    # TODO add condition: if no gradient distortion
    # https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh#L283-L284
    # fslroi "$fMRIFolder"/"$NameOffMRI"_gdc "$fMRIFolder"/"$NameOffMRI"_gdc_warp 0 3
    extract_func_roi = pe.Node(interface=fsl.ExtractROI(),
        name=f'extract_func_roi_{pipe_num}')

    extract_func_roi.inputs.t_min = 0
    extract_func_roi.inputs.t_size = 3

    node, out = strat_pool.get_data('bold')
    wf.connect(node, out, extract_func_roi, 'in_file')

    # fslmaths "$fMRIFolder"/"$NameOffMRI"_gdc_warp -mul 0 "$fMRIFolder"/"$NameOffMRI"_gdc_warp
    multiply_func_roi_by_zero = pe.Node(interface=fsl.maths.MathsCommand(),
                                        name=f'multiply_func_roi_by_zero_{pipe_num}')

    multiply_func_roi_by_zero.inputs.args = '-mul 0'

    wf.connect(extract_func_roi, 'roi_file',
        multiply_func_roi_by_zero, 'in_file')

    # https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/scripts/OneStepResampling.sh#L168-L193
    # fslsplit ${InputfMRI} ${WD}/prevols/vol -t
    split_func = pe.Node(interface=fsl.Split(),
        name=f'split_func_{pipe_num}')

    split_func.inputs.dimension = 't'

    node, out = strat_pool.get_data('bold')
    wf.connect(node, out, split_func, 'in_file')

    ### Loop starts! ###
    # convertwarp --relout --rel --ref=${WD}/prevols/vol${vnum}.nii.gz --warp1=${GradientDistortionField} --postmat=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum} --out=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_gdc_warp.nii.gz
    convert_motion_distortion_warp = pe.MapNode(interface=fsl.ConvertWarp(),
        name=f'convert_motion_distortion_warp_{pipe_num}',
        iterfield=['reference', 'postmat'])

    convert_motion_distortion_warp.inputs.out_relwarp = True
    convert_motion_distortion_warp.inputs.relwarp = True

    wf.connect(multiply_func_roi_by_zero, 'out_file',
        convert_motion_distortion_warp, 'warp1')

    wf.connect(split_func, 'out_files',
        convert_motion_distortion_warp, 'reference')

    node, out = strat_pool.get_data('coordinate-transformation')
    wf.connect(node, out, convert_motion_distortion_warp, 'postmat')

    # convertwarp --relout --rel --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --warp1=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_gdc_warp.nii.gz --warp2=${OutputTransform} --out=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz
    convert_registration_warp = pe.MapNode(interface=fsl.ConvertWarp(),
        name=f'convert_registration_warp_{pipe_num}',
        iterfield=['warp1'])

    convert_registration_warp.inputs.out_relwarp = True
    convert_registration_warp.inputs.relwarp = True

    node, out = strat_pool.get_data('space-template_res-bold_desc-brain_T1w')
    wf.connect(node, out, convert_registration_warp, 'reference')

    wf.connect(convert_motion_distortion_warp, 'out_file',
        convert_registration_warp, 'warp1')

    wf.connect(convert_func_to_standard_warp, 'out_file',
        convert_registration_warp, 'warp2')

    # fslmaths ${WD}/prevols/vol${vnum}.nii.gz -mul 0 -add 1 ${WD}/prevols/vol${vnum}_mask.nii.gz
    generate_vol_mask = pe.MapNode(interface=fsl.maths.MathsCommand(),
                        name=f'generate_mask_{pipe_num}',
                        iterfield=['in_file'])

    generate_vol_mask.inputs.args = '-mul 0 -add 1'

    wf.connect(split_func, 'out_files',
        generate_vol_mask, 'in_file')

    # applywarp --rel --interp=spline --in=${WD}/prevols/vol${vnum}.nii.gz --warp=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${WD}/postvols/vol${vnum}.nii.gz
    applywarp_func_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                                    name=f'applywarp_func_to_standard_{pipe_num}',
                                    iterfield=['in_file', 'field_file'])

    applywarp_func_to_standard.inputs.relwarp = True
    applywarp_func_to_standard.inputs.interp = 'spline'

    wf.connect(split_func, 'out_files',
        applywarp_func_to_standard, 'in_file')

    wf.connect(convert_registration_warp, 'out_file',
        applywarp_func_to_standard, 'field_file')

    node, out = strat_pool.get_data('space-template_res-bold_desc-brain_T1w')
    wf.connect(node, out,
        applywarp_func_to_standard, 'ref_file')

    # applywarp --rel --interp=nn --in=${WD}/prevols/vol${vnum}_mask.nii.gz --warp=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${WD}/postvols/vol${vnum}_mask.nii.gz
    applywarp_func_mask_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                                    name=f'applywarp_func_mask_to_standard_{pipe_num}',
                                    iterfield=['in_file', 'field_file'])

    applywarp_func_mask_to_standard.inputs.relwarp = True
    applywarp_func_mask_to_standard.inputs.interp = 'nn'

    wf.connect(generate_vol_mask, 'out_file',
        applywarp_func_mask_to_standard, 'in_file')

    wf.connect(convert_registration_warp, 'out_file',
        applywarp_func_mask_to_standard, 'field_file')

    node, out = strat_pool.get_data('space-template_res-bold_desc-brain_T1w')
    wf.connect(node, out,
        applywarp_func_mask_to_standard, 'ref_file')

    ### Loop ends! ###

    # fslmerge -tr ${OutputfMRI} $FrameMergeSTRING $TR_vol
    merge_func_to_standard = pe.Node(interface=fslMerge(),
                                     name=f'merge_func_to_standard_{pipe_num}')

    merge_func_to_standard.inputs.dimension = 't'

    wf.connect(applywarp_func_to_standard, 'out_file',
        merge_func_to_standard, 'in_files')

    # fslmerge -tr ${OutputfMRI}_mask $FrameMergeSTRINGII $TR_vol
    merge_func_mask_to_standard = pe.Node(interface=fslMerge(),
                                          name='merge_func_mask_to_'
                                               f'standard_{pipe_num}')

    merge_func_mask_to_standard.inputs.dimension = 't'

    wf.connect(applywarp_func_mask_to_standard, 'out_file',
        merge_func_mask_to_standard, 'in_files')

    # fslmaths ${OutputfMRI}_mask -Tmin ${OutputfMRI}_mask
    find_min_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                        name=f'find_min_mask_{pipe_num}')

    find_min_mask.inputs.args = '-Tmin'

    wf.connect(merge_func_mask_to_standard, 'merged_file',
        find_min_mask, 'in_file')

    # Combine transformations: gradient non-linearity distortion + fMRI_dc to standard
    # convertwarp --relout --rel --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --warp1=${GradientDistortionField} --warp2=${OutputTransform} --out=${WD}/Scout_gdc_MNI_warp.nii.gz
    convert_dc_warp = pe.Node(interface=fsl.ConvertWarp(),
        name=f'convert_dc_warp_{pipe_num}')

    convert_dc_warp.inputs.out_relwarp = True
    convert_dc_warp.inputs.relwarp = True

    node, out = strat_pool.get_data('space-template_res-bold_desc-brain_T1w')
    wf.connect(node, out, convert_dc_warp, 'reference')

    wf.connect(multiply_func_roi_by_zero, 'out_file',
        convert_dc_warp, 'warp1')

    wf.connect(convert_func_to_standard_warp, 'out_file',
        convert_dc_warp, 'warp2')

    # applywarp --rel --interp=spline --in=${ScoutInput} -w ${WD}/Scout_gdc_MNI_warp.nii.gz -r ${WD}/${T1wImageFile}.${FinalfMRIResolution} -o ${ScoutOutput}
    applywarp_scout = pe.Node(interface=fsl.ApplyWarp(),
        name=f'applywarp_scout_input_{pipe_num}')

    applywarp_scout.inputs.relwarp = True
    applywarp_scout.inputs.interp = 'spline'

    node, out = strat_pool.get_data('motion-basefile')
    wf.connect(node, out, applywarp_scout, 'in_file')

    node, out = strat_pool.get_data('space-template_res-bold_desc-brain_T1w')
    wf.connect(node, out, applywarp_scout, 'ref_file')

    wf.connect(convert_dc_warp, 'out_file', applywarp_scout, 'field_file')

    # https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/scripts/IntensityNormalization.sh#L124-L127
    # fslmaths ${InputfMRI} -mas ${BrainMask} -mas ${InputfMRI}_mask -thr 0 -ing 10000 ${OutputfMRI} -odt float
    merge_func_mask = pe.Node(util.Merge(2),
        name=f'merge_func_mask_{pipe_num}')

    node, out = strat_pool.get_data('space-template_desc-bold_mask')
    wf.connect(node, out, merge_func_mask, 'in1')

    wf.connect(find_min_mask, 'out_file', merge_func_mask, 'in2')

    extract_func_brain = pe.Node(interface=fsl.MultiImageMaths(),
                        name=f'extract_func_brain_{pipe_num}')

    extract_func_brain.inputs.op_string = '-mas %s -mas %s -thr 0 -ing 10000'
    extract_func_brain.inputs.output_datatype = 'float'

    wf.connect(merge_func_to_standard, 'merged_file',
        extract_func_brain, 'in_file')

    wf.connect(merge_func_mask, 'out',
        extract_func_brain, 'operand_files')

    # fslmaths ${ScoutInput} -mas ${BrainMask} -mas ${InputfMRI}_mask -thr 0 -ing 10000 ${ScoutOutput} -odt float
    extract_scout_brain = pe.Node(interface=fsl.MultiImageMaths(),
        name=f'extract_scout_brain_{pipe_num}')

    extract_scout_brain.inputs.op_string = '-mas %s -mas %s -thr 0 -ing 10000'
    extract_scout_brain.inputs.output_datatype = 'float'

    wf.connect(applywarp_scout, 'out_file',
        extract_scout_brain, 'in_file')

    wf.connect(merge_func_mask, 'out',
        extract_scout_brain, 'operand_files')

    outputs = {
        'space-template_desc-preproc_bold': (extract_func_brain, 'out_file'),
        'space-template_desc-scout_bold': (extract_scout_brain, 'out_file'),
        'space-template_desc-head_bold': (merge_func_to_standard, 'merged_file')
    }

    return (wf, outputs)


@nodeblock(
    name="transform_timeseries_to_T1template_dcan_nhp",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run"],
    option_key=["apply_transform", "using"],
    option_val="dcan_nhp",
    inputs=[
        (
            ["desc-reorient_bold", "bold"],
            "coordinate-transformation",
            "from-T1w_to-template_mode-image_warp",
            "from-bold_to-T1w_mode-image_desc-linear_warp",
            "T1w-template",
            "space-template_desc-head_T1w",
            "space-template_desc-T1w_mask",
            "space-template_desc-T1wT2w_biasfield",
        )
    ],
    outputs={
        "space-template_desc-preproc_bold": {"Template": "T1w-template"},
        "space-template_desc-bold_mask": {"Template": "T1w-template"},
    },
)
def warp_timeseries_to_T1template_dcan_nhp(wf, cfg, strat_pool, pipe_num,
                                           opt=None):
    # Apply motion correction, coreg, anat-to-template transforms on raw functional timeseries
    # Ref: https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/scripts/OneStepResampling.sh

    # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/scripts/OneStepResampling.sh#L131
    # ${FSLDIR}/bin/flirt -interp spline -in ${T1wImage} -ref ${T1wImage} -applyisoxfm $FinalfMRIResolution -out ${WD}/${T1wImageFile}.${FinalfMRIResolution}
    anat_resample = pe.Node(interface=fsl.FLIRT(),
                            name=f'anat_resample_func_res_{pipe_num}'
                            )
    anat_resample.inputs.apply_isoxfm = float(cfg.registration_workflows['functional_registration']['func_registration_to_template']['output_resolution']['func_preproc_outputs'].replace("mm", ""))
    anat_resample.inputs.interp = 'spline'

    node, out = strat_pool.get_data('space-template_desc-head_T1w')
    wf.connect(node, out, anat_resample, 'in_file')
    wf.connect(node, out, anat_resample, 'reference')

    # ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${T1wImage} -r ${ResampRefIm} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${WD}/${T1wImageFile}.${FinalfMRIResolution}
    applywarp_anat_res = pe.Node(interface=fsl.ApplyWarp(),
                                name=f'anat_func_res_{pipe_num}')

    applywarp_anat_res.inputs.relwarp = True
    applywarp_anat_res.inputs.interp = 'spline'
    applywarp_anat_res.inputs.premat = cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['identity_matrix']

    node, out = strat_pool.get_data('space-template_desc-head_T1w')
    wf.connect(node, out, applywarp_anat_res, 'in_file')
    wf.connect(anat_resample, 'out_file', applywarp_anat_res, 'ref_file')

    # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/scripts/OneStepResampling.sh#L136-L138
    # Create brain masks in this space (changing resolution)
    # ${FSLDIR}/bin/applywarp --rel --interp=nn -i ${FreeSurferBrainMask}.nii.gz -r ${WD}/${T1wImageFile}.${FinalfMRIResolution} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz
    applywarp_anat_mask_res = pe.Node(interface=fsl.ApplyWarp(),
                                name=f'anat_mask_func_res_{pipe_num}')
    applywarp_anat_mask_res.inputs.relwarp = True
    applywarp_anat_mask_res.inputs.interp = 'nn'
    applywarp_anat_mask_res.inputs.premat = cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['identity_matrix']

    node, out = strat_pool.get_data('space-template_desc-T1w_mask')
    wf.connect(node, out, applywarp_anat_mask_res, 'in_file')
    wf.connect(applywarp_anat_res, 'out_file', applywarp_anat_mask_res, 'ref_file')

    # ${FSLDIR}/bin/fslmaths ${WD}/${T1wImageFile}.${FinalfMRIResolution} -mas ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz
    T1_brain_res = pe.Node(interface=fsl.MultiImageMaths(),
                                  name=f't1_brain_func_res_{pipe_num}')
    T1_brain_res.inputs.op_string = "-mas %s "

    wf.connect(applywarp_anat_res, 'out_file', T1_brain_res, 'in_file')
    wf.connect(applywarp_anat_mask_res, 'out_file', T1_brain_res, 'operand_files')

    # Create versions of the biasfield (changing resolution)
    # ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${BiasField} -r ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${WD}/${BiasFieldFile}.${FinalfMRIResolution}
    applywarp_bias_field_res = pe.Node(interface=fsl.ApplyWarp(),
                                name=f'biasfiled_func_res_{pipe_num}')
    applywarp_bias_field_res.inputs.relwarp = True
    applywarp_bias_field_res.inputs.interp = 'spline'
    applywarp_bias_field_res.inputs.premat = cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['identity_matrix']

    node, out = strat_pool.get_data('space-template_desc-T1wT2w_biasfield')
    wf.connect(node, out, applywarp_bias_field_res, 'in_file')
    wf.connect(T1_brain_res, 'out_file', applywarp_bias_field_res, 'ref_file')

    # ${FSLDIR}/bin/fslmaths ${WD}/${BiasFieldFile}.${FinalfMRIResolution} -thr 0.1 ${WD}/${BiasFieldFile}.${FinalfMRIResolution}
    biasfield_thr = pe.Node(interface=fsl.MultiImageMaths(),
                                  name=f'biasfiedl_thr_{pipe_num}')
    biasfield_thr.inputs.op_string = "-thr 0.1"

    wf.connect(applywarp_bias_field_res, 'out_file', biasfield_thr, 'in_file')

    # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/scripts/OneStepResampling.sh#L144-L146
    # convertwarp --relout --rel --warp1=${fMRIToStructuralInput} --warp2=${StructuralToStandard} --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${OutputTransform}
    convert_func_to_standard_warp = pe.Node(interface=fsl.ConvertWarp(),
        name=f'convert_func_to_standard_warp_{pipe_num}')

    convert_func_to_standard_warp.inputs.out_relwarp = True
    convert_func_to_standard_warp.inputs.relwarp = True

    node, out = strat_pool.get_data('from-bold_to-T1w_mode-image_desc-linear_warp')
    wf.connect(node, out, convert_func_to_standard_warp, 'warp1')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_warp')
    wf.connect(node, out, convert_func_to_standard_warp, 'warp2')

    wf.connect(applywarp_anat_res, 'out_file', convert_func_to_standard_warp, 'reference')

    # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/GenericfMRIVolumeProcessingPipeline.sh#L157-L158
    # fslroi "$fMRIFolder"/"$NameOffMRI"_gdc "$fMRIFolder"/"$NameOffMRI"_gdc_warp 0 3
    extract_func_roi = pe.Node(interface=fsl.ExtractROI(),
        name=f'extract_func_roi_{pipe_num}')

    extract_func_roi.inputs.t_min = 0
    extract_func_roi.inputs.t_size = 3

    node, out = strat_pool.get_data(['desc-reorient_bold', 'bold'])
    wf.connect(node, out, extract_func_roi, 'in_file')

    # fslmaths "$fMRIFolder"/"$NameOffMRI"_gdc_warp -mul 0 "$fMRIFolder"/"$NameOffMRI"_gdc_warp
    multiply_func_roi_by_zero = pe.Node(interface=fsl.maths.MathsCommand(),
                                        name=f'multiply_func_roi_by_zero_{pipe_num}')

    multiply_func_roi_by_zero.inputs.args = '-mul 0'

    wf.connect(extract_func_roi, 'roi_file',
        multiply_func_roi_by_zero, 'in_file')

    # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/scripts/OneStepResampling.sh#L173
    # fslsplit ${InputfMRI} ${WD}/prevols/vol -t
    split_func = pe.Node(interface=fsl.Split(),
        name=f'split_func_{pipe_num}')

    split_func.inputs.dimension = 't'

    node, out = strat_pool.get_data(['desc-reorient_bold', 'bold'])
    wf.connect(node, out, split_func, 'in_file')

    ### Loop starts! ###
    # convertwarp --relout --rel --ref=${WD}/prevols/vol${vnum}.nii.gz --warp1=${GradientDistortionField} --postmat=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum} --out=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_gdc_warp.nii.gz
    convert_motion_distortion_warp = pe.MapNode(interface=fsl.ConvertWarp(),
        name=f'convert_motion_distortion_warp_{pipe_num}',
        iterfield=['reference', 'postmat'])

    convert_motion_distortion_warp.inputs.out_relwarp = True
    convert_motion_distortion_warp.inputs.relwarp = True

    wf.connect(multiply_func_roi_by_zero, 'out_file',
        convert_motion_distortion_warp, 'warp1')

    wf.connect(split_func, 'out_files',
        convert_motion_distortion_warp, 'reference')

    node, out = strat_pool.get_data('coordinate-transformation')
    wf.connect(node, out, convert_motion_distortion_warp, 'postmat')

    # convertwarp --relout --rel --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --warp1=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_gdc_warp.nii.gz --warp2=${OutputTransform} --out=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz
    convert_registration_warp = pe.MapNode(interface=fsl.ConvertWarp(),
        name=f'convert_registration_warp_{pipe_num}',
        iterfield=['warp1'])

    convert_registration_warp.inputs.out_relwarp = True
    convert_registration_warp.inputs.relwarp = True

    wf.connect(applywarp_anat_res, 'out_file', convert_registration_warp, 'reference')

    wf.connect(convert_motion_distortion_warp, 'out_file',
        convert_registration_warp, 'warp1')

    wf.connect(convert_func_to_standard_warp, 'out_file',
        convert_registration_warp, 'warp2')

    # fslmaths ${WD}/prevols/vol${vnum}.nii.gz -mul 0 -add 1 ${WD}/prevols/vol${vnum}_mask.nii.gz
    generate_vol_mask = pe.MapNode(interface=fsl.maths.MathsCommand(),
                        name=f'generate_mask_{pipe_num}',
                        iterfield=['in_file'])

    generate_vol_mask.inputs.args = '-mul 0 -add 1'

    wf.connect(split_func, 'out_files',
        generate_vol_mask, 'in_file')

    # applywarp --rel --interp=spline --in=${WD}/prevols/vol${vnum}.nii.gz --warp=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${WD}/postvols/vol${vnum}.nii.gz
    applywarp_func_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                                    name=f'applywarp_func_to_standard_{pipe_num}',
                                    iterfield=['in_file', 'field_file'])

    applywarp_func_to_standard.inputs.relwarp = True
    applywarp_func_to_standard.inputs.interp = 'spline'

    wf.connect(split_func, 'out_files',
        applywarp_func_to_standard, 'in_file')

    wf.connect(convert_registration_warp, 'out_file',
        applywarp_func_to_standard, 'field_file')

    wf.connect(applywarp_anat_res, 'out_file',
        applywarp_func_to_standard, 'ref_file')

    # applywarp --rel --interp=nn --in=${WD}/prevols/vol${vnum}_mask.nii.gz --warp=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_all_warp.nii.gz --ref=${WD}/${T1wImageFile}.${FinalfMRIResolution} --out=${WD}/postvols/vol${vnum}_mask.nii.gz
    applywarp_func_mask_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                                    name=f'applywarp_func_mask_to_standard_{pipe_num}',
                                    iterfield=['in_file', 'field_file'])

    applywarp_func_mask_to_standard.inputs.relwarp = True
    applywarp_func_mask_to_standard.inputs.interp = 'nn'

    wf.connect(generate_vol_mask, 'out_file',
        applywarp_func_mask_to_standard, 'in_file')

    wf.connect(convert_registration_warp, 'out_file',
        applywarp_func_mask_to_standard, 'field_file')

    wf.connect(applywarp_anat_res, 'out_file',
        applywarp_func_mask_to_standard, 'ref_file')

    ### Loop ends! ###

    # fslmerge -tr ${OutputfMRI} $FrameMergeSTRING $TR_vol
    merge_func_to_standard = pe.Node(interface=fslMerge(),
                                     name=f'merge_func_to_standard_{pipe_num}')

    merge_func_to_standard.inputs.dimension = 't'

    wf.connect(applywarp_func_to_standard, 'out_file',
        merge_func_to_standard, 'in_files')

    # fslmerge -tr ${OutputfMRI}_mask $FrameMergeSTRINGII $TR_vol
    merge_func_mask_to_standard = pe.Node(interface=fslMerge(),
                                          name='merge_func_mask_to_'
                                               f'standard_{pipe_num}')

    merge_func_mask_to_standard.inputs.dimension = 't'

    wf.connect(applywarp_func_mask_to_standard, 'out_file',
        merge_func_mask_to_standard, 'in_files')

    # fslmaths ${OutputfMRI}_mask -Tmin ${OutputfMRI}_mask
    find_min_mask = pe.Node(interface=fsl.maths.MathsCommand(),
        name=f'find_min_mask_{pipe_num}')

    find_min_mask.inputs.args = '-Tmin'

    wf.connect(merge_func_mask_to_standard, 'merged_file',
        find_min_mask, 'in_file')

    # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/fMRIVolume/scripts/IntensityNormalization.sh#L113-L119
    # fslmaths ${InputfMRI} -div ${BiasField} $jacobiancom -mas ${BrainMask} -mas ${InputfMRI}_mask -ing 10000 ${OutputfMRI} -odt float

    merge_func_mask = pe.Node(util.Merge(3),
                                name=f'merge_operand_files_{pipe_num}')

    wf.connect(biasfield_thr, 'out_file', merge_func_mask, 'in1')

    wf.connect(applywarp_anat_mask_res, 'out_file', merge_func_mask, 'in2')

    wf.connect(find_min_mask, 'out_file', merge_func_mask, 'in3')


    extract_func_brain = pe.Node(interface=fsl.MultiImageMaths(),
                        name=f'extract_func_brain_{pipe_num}')

    extract_func_brain.inputs.op_string = '-div %s -mas %s -mas %s -ing 10000'
    extract_func_brain.inputs.output_datatype = 'float'

    wf.connect(merge_func_to_standard, 'merged_file',
        extract_func_brain, 'in_file')

    wf.connect(merge_func_mask, 'out',
        extract_func_brain, 'operand_files')

    func_mask_final = pe.Node(interface=fsl.MultiImageMaths(),
                                name=f'func_mask_final_{pipe_num}')
    func_mask_final.inputs.op_string = "-mas %s "

    wf.connect(applywarp_anat_mask_res, 'out_file', func_mask_final, 'in_file')

    wf.connect(find_min_mask, 'out_file', func_mask_final, 'operand_files')

    outputs = {
        'space-template_desc-preproc_bold': (extract_func_brain, 'out_file'),
        'space-template_desc-bold_mask': (func_mask_final, 'out_file')
    }

    return (wf, outputs)


@nodeblock(
    name="single_step_resample_stc_timeseries_to_T1template",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run"],
    option_key=["apply_transform", "using"],
    option_val="single_step_resampling_from_stc",
    inputs=[
        (
            "sbref",
            "desc-stc_bold",
            "motion-basefile",
            "space-bold_desc-brain_mask",
            "coordinate-transformation",
            "from-T1w_to-template_mode-image_xfm",
            "from-bold_to-T1w_mode-image_desc-linear_xfm",
            "from-bold_to-template_mode-image_xfm",
            "ants-blip-warp",
            "fsl-blip-warp",
            "T1w",
            "desc-preproc_T1w",
            "T1w-brain-template-funcreg",
            "T1w-brain-template-deriv",
        )
    ],
    outputs={
        "space-template_desc-preproc_bold": {
            "Template": "T1w-brain-template-funcreg"},
        "space-template_desc-brain_bold": {
            "Template": "T1w-brain-template-funcreg"},
        "space-template_desc-bold_mask": {
            "Template": "T1w-brain-template-funcreg"},
        "space-template_desc-head_bold": {
            "Template": "T1w-brain-template-funcreg"},
        "space-template_res-derivative_desc-preproc_bold": {
            "Template": "T1w-brain-template-deriv"
        },
        "space-template_res-derivative_desc-bold_mask": {
            "Template": "T1w-brain-template-deriv"
        },
    },
)
def single_step_resample_timeseries_to_T1template(wf, cfg, strat_pool,
                                                  pipe_num, opt=None):
    '''
    Apply motion correction, coreg, anat-to-template transforms on
    slice-time corrected functional timeseries based on fMRIPrep
    pipeline

    Copyright (c) 2015-2018, the CRN developers team.
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

    * Neither the name of fmriprep nor the names of its contributors
    may be used to endorse or promote products derived from this
    software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
    OF THE POSSIBILITY OF SUCH DAMAGE.

    Ref: https://github.com/nipreps/fmriprep/blob/84a6005b/fmriprep/workflows/bold/resampling.py#L159-L419
    '''  # noqa: 501
    xfm_prov = strat_pool.get_cpac_provenance(
        'from-T1w_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    bbr2itk = pe.Node(util.Function(input_names=['reference_file',
                                                 'source_file',
                                                 'transform_file'],
                                    output_names=['itk_transform'],
                                    function=run_c3d),
                      name=f'convert_bbr2itk_{pipe_num}')

    if cfg.registration_workflows['functional_registration'][
            'coregistration']['boundary_based_registration'][
            'reference'] == 'whole-head':
        node, out = strat_pool.get_data('T1w')
        wf.connect(node, out, bbr2itk, 'reference_file')

    elif cfg.registration_workflows['functional_registration'][
            'coregistration']['boundary_based_registration'][
            'reference'] == 'brain':
        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, bbr2itk, 'reference_file')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, bbr2itk, 'source_file')

    node, out = strat_pool.get_data(
        'from-bold_to-T1w_mode-image_desc-linear_xfm')
    wf.connect(node, out, bbr2itk, 'transform_file')

    split_func = pe.Node(interface=fsl.Split(),
        name=f'split_func_{pipe_num}')

    split_func.inputs.dimension = 't'

    node, out = strat_pool.get_data('desc-stc_bold')
    wf.connect(node, out, split_func, 'in_file')

    ### Loop starts! ###
    motionxfm2itk = pe.MapNode(util.Function(
        input_names=['reference_file',
                     'source_file',
                     'transform_file'],
        output_names=['itk_transform'],
        function=run_c3d),
        name=f'convert_motionxfm2itk_{pipe_num}',
        iterfield=['transform_file'])

    node, out = strat_pool.get_data('motion-basefile')
    wf.connect(node, out, motionxfm2itk, 'reference_file')
    wf.connect(node, out, motionxfm2itk, 'source_file')

    node, out = strat_pool.get_data('coordinate-transformation')
    motion_correct_tool = check_prov_for_motion_tool(
        strat_pool.get_cpac_provenance('coordinate-transformation'))
    if motion_correct_tool == 'mcflirt':
        wf.connect(node, out, motionxfm2itk, 'transform_file')
    elif motion_correct_tool == '3dvolreg':
        convert_transform = pe.Node(util.Function(
            input_names=['one_d_filename'],
            output_names=['transform_directory'],
            function=one_d_to_mat,
            imports=['import os', 'import numpy as np']),
            name=f'convert_transform_{pipe_num}')
        wf.connect(node, out, convert_transform, 'one_d_filename')
        wf.connect(convert_transform, 'transform_directory',
                   motionxfm2itk, 'transform_file')

    merge_num = 4
    blip = False
    if strat_pool.check_rpool('ants-blip-warp') and reg_tool == 'ants':
        blip_node, blip_out = strat_pool.get_data('ants-blip-warp')
        merge_num = 5
        blip = True
    elif strat_pool.check_rpool('fsl-blip-warp') and reg_tool == 'fsl':
        blip_node, blip_out = strat_pool.get_data('fsl-blip-warp')
        merge_num = 5
        blip = True

    collectxfm = pe.MapNode(util.Merge(merge_num),
                            name=f'collectxfm_func_to_standard_{pipe_num}',
                            iterfield=[f'in{merge_num}'])

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, collectxfm, 'in1')

    wf.connect(bbr2itk, 'itk_transform',
               collectxfm, 'in2')

    collectxfm.inputs.in3 = 'identity'

    if blip:
        wf.connect(blip_node, blip_out, collectxfm, 'in4')

    wf.connect(motionxfm2itk, 'itk_transform',
               collectxfm, f'in{merge_num}')

    applyxfm_func_to_standard = pe.MapNode(
        interface=ants.ApplyTransforms(),
        name=f'applyxfm_func_to_standard_{pipe_num}',
        iterfield=['input_image', 'transforms'])
    applyxfm_func_to_standard.inputs.float = True
    applyxfm_func_to_standard.inputs.interpolation = 'LanczosWindowedSinc'

    applyxfm_derivfunc_to_standard = pe.MapNode(
        interface=ants.ApplyTransforms(),
        name=f'applyxfm_derivfunc_to_standard_{pipe_num}',
        iterfield=['input_image', 'transforms'])
    applyxfm_derivfunc_to_standard.inputs.float = True
    applyxfm_derivfunc_to_standard.inputs.interpolation = 'LanczosWindowedSinc'

    wf.connect(split_func, 'out_files',
               applyxfm_func_to_standard, 'input_image')
    wf.connect(split_func, 'out_files',
               applyxfm_derivfunc_to_standard, 'input_image')

    node, out = strat_pool.get_data('T1w-brain-template-funcreg')
    wf.connect(node, out, applyxfm_func_to_standard, 'reference_image')
    
    node, out = strat_pool.get_data('T1w-brain-template-deriv')
    wf.connect(node, out, applyxfm_derivfunc_to_standard, 'reference_image')

    wf.connect(collectxfm, 'out', applyxfm_func_to_standard, 'transforms')
    wf.connect(collectxfm, 'out', applyxfm_derivfunc_to_standard, 'transforms')

    ### Loop ends! ###

    merge_func_to_standard = pe.Node(interface=fslMerge(),
                                     name=f'merge_func_to_standard_{pipe_num}')
    merge_func_to_standard.inputs.dimension = 't'

    wf.connect(applyxfm_func_to_standard, 'output_image',
               merge_func_to_standard, 'in_files')

    merge_derivfunc_to_standard = pe.Node(
        interface=fslMerge(), name=f'merge_derivfunc_to_standard_{pipe_num}')
    merge_derivfunc_to_standard.inputs.dimension = 't'

    wf.connect(applyxfm_derivfunc_to_standard, 'output_image',
               merge_derivfunc_to_standard, 'in_files')

    applyxfm_func_mask_to_standard = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'applyxfm_func_mask_to_standard_{pipe_num}')
    applyxfm_func_mask_to_standard.inputs.interpolation = 'MultiLabel'

    node, out = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out, applyxfm_func_mask_to_standard, 'input_image')

    node, out = strat_pool.get_data('T1w-brain-template-funcreg')
    wf.connect(node, out, applyxfm_func_mask_to_standard, 'reference_image')

    collectxfm_mask = pe.Node(
        util.Merge(2), name=f'collectxfm_func_mask_to_standard_{pipe_num}')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, collectxfm_mask, 'in1')

    wf.connect(bbr2itk, 'itk_transform', collectxfm_mask, 'in2')

    wf.connect(collectxfm_mask, 'out',
               applyxfm_func_mask_to_standard, 'transforms')

    applyxfm_deriv_mask_to_standard = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'applyxfm_deriv_mask_to_standard_{pipe_num}')
    applyxfm_deriv_mask_to_standard.inputs.interpolation = 'MultiLabel'

    node, out = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out, applyxfm_deriv_mask_to_standard, 'input_image')

    node, out = strat_pool.get_data('T1w-brain-template-deriv')
    wf.connect(node, out, applyxfm_deriv_mask_to_standard, 'reference_image')

    collectxfm_deriv_mask = pe.Node(
        util.Merge(2), name=f'collectxfm_deriv_mask_to_standard_{pipe_num}')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, collectxfm_deriv_mask, 'in1')

    wf.connect(bbr2itk, 'itk_transform',
               collectxfm_deriv_mask, 'in2')

    wf.connect(collectxfm_deriv_mask, 'out',
               applyxfm_deriv_mask_to_standard, 'transforms')

    apply_mask = pe.Node(interface=fsl.maths.ApplyMask(),
                         name=f'get_func_brain_to_standard_{pipe_num}')

    wf.connect(merge_func_to_standard, 'merged_file',
               apply_mask, 'in_file')

    wf.connect(applyxfm_func_mask_to_standard, 'output_image',
               apply_mask, 'mask_file')

    outputs = {
        'space-template_desc-head_bold': (merge_func_to_standard,
                                          'merged_file'),
        'space-template_desc-brain_bold': (apply_mask, 'out_file'),
        'space-template_desc-preproc_bold': (apply_mask, 'out_file'),
        'space-template_desc-bold_mask': (applyxfm_func_mask_to_standard,
                                          'output_image'),
        'space-template_res-derivative_desc-preproc_bold':
            (merge_derivfunc_to_standard, 'merged_file'),
        'space-template_res-derivative_desc-bold_mask':
            (applyxfm_deriv_mask_to_standard, 'output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="transform_sbref_to_T1template",
    switch=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
        "run",
    ],
    inputs=[
        ("sbref", "from-bold_to-template_mode-image_xfm"),
        "T1w-brain-template-funcreg",
    ],
    outputs={
        "space-template_sbref": {
            "Description": "Single-volume sbref of the BOLD time-series "
                           "transformed to template space.",
            "Template": "T1w-brain-template-funcreg",
        }
    },
)
def warp_sbref_to_T1template(wf, cfg, strat_pool, pipe_num, opt=None):
    xfm = 'from-bold_to-template_mode-image_xfm'
    wf, apply_xfm = warp_resource_to_template(
        wf, cfg, strat_pool, pipe_num, 'sbref', xfm,
        reference='T1w-brain-template-funcreg', time_series=False)[:2]
    outputs = {'space-template_sbref':
               (apply_xfm, 'outputspec.output_image')}
    return _warp_return(wf, apply_xfm, outputs)


@nodeblock(
    name="transform_bold_mask_to_T1template",
    switch=[
        [
            "registration_workflows",
            "functional_registration",
            "func_registration_to_template",
            "run",
        ],
        ["registration_workflows", "anatomical_registration", "run"],
    ],
    option_key=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
        "apply_transform",
        "using",
    ],
    option_val=["default", "abcd", "dcan_nhp"],
    inputs=[
        ("space-bold_desc-brain_mask", "from-bold_to-template_mode-image_xfm"),
        "T1w-brain-template-funcreg",
    ],
    outputs={
        "space-template_desc-bold_mask": {
            "Template": "T1w-brain-template-funcreg"}})
def warp_bold_mask_to_T1template(wf, cfg, strat_pool, pipe_num, opt=None):
    xfm = 'from-bold_to-template_mode-image_xfm'
    wf, apply_xfm = warp_resource_to_template(
        wf, cfg, strat_pool, pipe_num, 'space-bold_desc-brain_mask', xfm,
        reference='T1w-brain-template-funcreg', time_series=False)[:2]
    outputs = {'space-template_desc-bold_mask':
               (apply_xfm, 'outputspec.output_image')}
    return _warp_return(wf, apply_xfm, outputs)


@nodeblock(
    name="transform_deriv_mask_to_T1template",
    switch=[
        [
            "registration_workflows",
            "functional_registration",
            "func_registration_to_template",
            "run",
        ],
        ["registration_workflows", "anatomical_registration", "run"],
    ],
    option_key=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
        "apply_transform",
        "using",
    ],
    option_val=["default", "abcd", "dcan_nhp"],
    inputs=[
        ("space-bold_desc-brain_mask", "from-bold_to-template_mode-image_xfm"),
        "T1w-brain-template-deriv",
    ],
    outputs={
        "space-template_res-derivative_desc-bold_mask": {
            "Template": "T1w-brain-template-deriv"
        }
    },
)
def warp_deriv_mask_to_T1template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Transform the BOLD mask to template space and to the resolution set for
    the derivative outputs.
    '''
    xfm = 'from-bold_to-template_mode-image_xfm'
    wf, apply_xfm = warp_resource_to_template(
        wf, cfg, strat_pool, pipe_num, 'space-bold_desc-brain_mask', xfm,
        reference='T1w-brain-template-deriv', time_series=False)[:2]
    outputs = {'space-template_res-derivative_desc-bold_mask':
               (apply_xfm, 'outputspec.output_image')}
    return _warp_return(wf, apply_xfm, outputs)


@nodeblock(
    name="transform_timeseries_to_EPItemplate",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run_EPI"],
    inputs=[
        ("desc-preproc_bold", "from-bold_to-EPItemplate_mode-image_xfm"),
        "EPI-template",
    ],
    outputs={"space-template_desc-preproc_bold": {"Template": "EPI-template"}},
)
def warp_timeseries_to_EPItemplate(wf, cfg, strat_pool, pipe_num, opt=None):
    xfm = 'from-bold_to-EPItemplate_mode-image_xfm'
    wf, apply_xfm, resource = warp_resource_to_template(
        wf, cfg, strat_pool, pipe_num, 'desc-preproc_bold', xfm,
        time_series=True)
    outputs = {f'space-template_{resource}':
               (apply_xfm, 'outputspec.output_image')}
    return _warp_return(wf, apply_xfm, outputs)


@nodeblock(
    name="transform_bold_mean_to_EPItemplate",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run_EPI"],
    inputs=[
        ("desc-mean_bold", "from-bold_to-EPItemplate_mode-image_xfm"),
        "EPI-template",
    ],
    outputs={"space-template_desc-mean_bold": {"Template": "EPI-template"}},
)
def warp_bold_mean_to_EPItemplate(wf, cfg, strat_pool, pipe_num, opt=None):
    xfm = 'from-bold_to-EPItemplate_mode-image_xfm'
    wf, apply_xfm = warp_resource_to_template(
        wf, cfg, strat_pool, pipe_num, 'desc-mean_bold', xfm,
        time_series=False)[:2]
    outputs = {'space-template_desc-mean_bold':
               (apply_xfm, 'outputspec.output_image')}
    return _warp_return(wf, apply_xfm, outputs)


@nodeblock(
    name="transform_bold_mask_to_EPItemplate",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run_EPI"],
    inputs=[
        ("space-bold_desc-brain_mask",
         "from-bold_to-EPItemplate_mode-image_xfm"),
        "EPI-template",
    ],
    outputs={"space-template_desc-bold_mask": {"Template": "EPI-template"}},
)
def warp_bold_mask_to_EPItemplate(wf, cfg, strat_pool, pipe_num, opt=None):
    xfm = 'from-bold_to-EPItemplate_mode-image_xfm'
    wf, apply_xfm = warp_resource_to_template(
        wf, cfg, strat_pool, pipe_num, 'space-bold_desc-brain_mask', xfm,
        time_series=False)[:2]
    outputs = {'space-template_desc-bold_mask':
               (apply_xfm, 'outputspec.output_image')}
    return _warp_return(wf, apply_xfm, outputs)


@nodeblock(
    name="transform_deriv_mask_to_EPItemplate",
    config=[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
    ],
    switch=["run_EPI"],
    inputs=[
        ("space-bold_desc-brain_mask",
         "from-bold_to-EPItemplate_mode-image_xfm"),
        "EPI-template",
    ],
    outputs={
        "space-template_res-derivative_desc-bold_mask": {
            "Template": "EPI-template"}
    },
)
def warp_deriv_mask_to_EPItemplate(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Transform the BOLD mask to template space and to the resolution set for
    the derivative outputs.
    '''
    xfm = 'from-bold_to-EPItemplate_mode-image_xfm'
    wf, apply_xfm = warp_resource_to_template(
        wf, cfg, strat_pool, pipe_num, 'space-bold_desc-brain_mask', xfm,
        time_series=False)[:2]
    outputs = {'space-template_res-derivative_desc-bold_mask':
               (apply_xfm, 'outputspec.output_image')}
    return _warp_return(wf, apply_xfm, outputs)


@nodeblock(
    name="warp_tissuemask_to_T1template",
    switch=["registration_workflows", "anatomical_registration", "run"],
    inputs=[
        (
            "label-CSF_mask",
            "label-WM_mask",
            "label-GM_mask",
            "from-T1w_to-template_mode-image_xfm",
        ),
        "T1w-template",
    ],
    outputs={
        "space-template_label-CSF_mask": {"Template": "T1w-template"},
        "space-template_label-WM_mask": {"Template": "T1w-template"},
        "space-template_label-GM_mask": {"Template": "T1w-template"},
    },
)
def warp_tissuemask_to_T1template(wf, cfg, strat_pool, pipe_num, opt=None):
    return warp_tissuemask_to_template(wf, cfg, strat_pool, pipe_num,
                                       xfm='from-T1w_to-template_mode-image_'
                                           'xfm', template_space='T1')


@nodeblock(
    name="warp_tissuemask_to_EPItemplate",
    switch=[
        "registration_workflows",
        "functional_registration",
        "EPI_registration",
        "run",
    ],
    inputs=[
        (
            "label-CSF_mask",
            "label-WM_mask",
            "label-GM_mask",
            "from-bold_to-EPItemplate_mode-image_xfm",
        ),
        "EPI-template",
    ],
    outputs={
        "space-template_label-CSF_mask": {"Template": "EPI-template"},
        "space-template_label-WM_mask": {"Template": "EPI-template"},
        "space-template_label-GM_mask": {"Template": "EPI-template"},
    },
)
def warp_tissuemask_to_EPItemplate(wf, cfg, strat_pool, pipe_num, opt=None):
    return warp_tissuemask_to_template(wf, cfg, strat_pool, pipe_num,
                                       xfm='from-bold_to-EPItemplate_'
                                           'mode-image_xfm',
                                       template_space='EPI')


def warp_tissuemask_to_template(wf, cfg, strat_pool, pipe_num, xfm,
                                template_space):
    '''Function to apply transforms to tissue masks

    Parameters
    ----------
    wf, cfg, strat_pool, pipe_num
        passed through from Node Block

    xfm : str
        transform

    template_space : str
        T1 or EPI

    Returns
    -------
    wf : nipype.pipeline.engine.workflows.Workflow

    outputs : dict
    '''
    tissue_types = ['CSF', 'WM', 'GM']
    apply_xfm = {}
    for tissue in tissue_types:
        wf, apply_xfm[tissue] = warp_resource_to_template(
            wf, cfg, strat_pool, pipe_num, f'label-{tissue}_mask', xfm,
            time_series=False)[:2]
    if template_space == 'T1':
        template_space = ''
    outputs = {f'space-{template_space}template_label-{tissue}_mask': (
        apply_xfm[tissue], 'outputspec.output_image') for
               tissue in tissue_types}
    return _warp_return(wf, apply_xfm, outputs)


def warp_resource_to_template(wf: pe.Workflow, cfg, strat_pool, pipe_num: int,
                              input_resource: LIST_OR_STR, xfm: str,
                              reference: Optional[str] = None,
                              time_series: Optional[bool] = False
                              ) -> TUPLE[pe.Workflow, pe.Workflow, str]:
    '''Function to warp a resource into a template space

    Parameters
    ----------
    wf : pe.Workflow

    cfg : CPAC.utils.configuration.Configuration

    strat_pool : CPAC.pipeline.engine.ResourcePool

    pipe_num : int

    input_resource : str or list
        key for the resource to warp to template

    xfm : str
        key for the transform to apply

    reference : str, optional
        key for reference if not using f'{template_space}-template'

    time_series : boolean, optional
        resource to transform is 4D?

    Returns
    -------
    wf : pe.Workflow
        original workflow with subworkflow to warp resource to template
        connected

    apply_xfm : pe.Workflow
        subworkflow added to warp resource to template

    resource : str
        key of input resource in strat_pool
    '''
    # determine space we're warping to
    template_space = xfm.split('_to-', 1)[1].split('template')[0]
    if template_space == '':
        template_space = 'T1w'
    # determine tool used for registration
    xfm_prov = strat_pool.get_cpac_provenance(xfm)
    reg_tool = check_prov_for_regtool(xfm_prov)
    # set 'resource'
    if strat_pool.check_rpool(input_resource):
        resource, input_resource = strat_pool.get_data(input_resource,
                                                       report_fetched=True)
    else:
        return wf, None, input_resource
    # set 'reference' if not passed and determine subworkflow name
    if reference is None:
        subwf_input_name = input_resource
        reference = f'{template_space}-template'
    else:
        subwf_input_name = '-'.join([
            reference.split('-')[-1].split('_')[-1],
            input_resource.split('-')[-1].split('_')[-1]])
    # set up 'apply_transform' subworkflow
    apply_xfm = apply_transform(f'warp_{subwf_input_name}_to_'
                                f'{template_space}template_{pipe_num}',
                                reg_tool, time_series=time_series,
                                num_cpus=cfg.pipeline_setup['system_config'][
                                    'max_cores_per_participant'],
                                num_ants_cores=cfg.pipeline_setup[
                                    'system_config']['num_ants_threads'])
    # set appropriate 'interpolation' input based on registration tool
    if reg_tool == 'ants':
        apply_xfm.inputs.inputspec.interpolation = 'NearestNeighbor'
    elif reg_tool == 'fsl':
        apply_xfm.inputs.inputspec.interpolation = 'nn'
    # connect nodes to subworkflow
    node, out = resource
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')
    node, out = strat_pool.get_data(reference)
    wf.connect(node, out, apply_xfm, 'inputspec.reference')
    node, out = strat_pool.get_data(xfm)
    wf.connect(node, out, apply_xfm, 'inputspec.transform')
    return wf, apply_xfm, input_resource


def _warp_return(wf: pe.Workflow, apply_xfm: Optional[pe.Workflow],
                 outputs: dict) -> TUPLE[pe.Workflow, dict]:
    """Check if we have a transform to apply; if not, don't add the outputs"""
    if apply_xfm is None:
        return wf, {}
    return wf, outputs
