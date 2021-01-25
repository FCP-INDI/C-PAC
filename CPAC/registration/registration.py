import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants

from nipype.interfaces.afni import utils as afni_utils

import nipype.interfaces.c3 as c3


from CPAC.anat_preproc.lesion_preproc import create_lesion_preproc

from CPAC.func_preproc.utils import chunk_ts, split_ts_chunks

from CPAC.registration.utils import seperate_warps_list, \
                                    check_transforms, \
                                    generate_inverse_transform_flags, \
                                    single_ants_xfm_to_list, \
                                    interpolation_string, \
                                    change_itk_transform_type, \
                                    hardcoded_reg

from CPAC.utils.utils import check_prov_for_regtool


def apply_transform(wf_name, reg_tool, time_series=False, multi_input=False,
                    num_cpus=1, num_ants_cores=1):

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

        input_image_type = 1
        if time_series:
            input_image_type = 3

        if multi_input:
            apply_warp = pe.MapNode(interface=ants.ApplyTransforms(),
                                    name=f'apply_warp_{wf_name}',
                                    iterfield=['input_image'])
        else:
            apply_warp = pe.Node(interface=ants.ApplyTransforms(),
                                 name=f'apply_warp_{wf_name}')

        apply_warp.inputs.dimension = 3
        apply_warp.inputs.input_image_type = input_image_type
        apply_warp.interface.num_threads = int(num_ants_cores)

        wf.connect(inputNode, 'reference', apply_warp, 'reference_image')

        interp_string = pe.Node(util.Function(input_names=['interpolation',
                                                           'reg_tool'],
                                              output_names=['interpolation'],
                                              function=interpolation_string),
                                name=f'interp_string')
        interp_string.inputs.reg_tool = reg_tool

        wf.connect(inputNode, 'interpolation', interp_string, 'interpolation')
        wf.connect(interp_string, 'interpolation',
                   apply_warp, 'interpolation')

        ants_xfm_list = \
            pe.Node(util.Function(input_names=['transform'],
                                  output_names=['transform_list'],
                                  function=single_ants_xfm_to_list),
                    name=f'single_ants_xfm_to_list')

        wf.connect(inputNode, 'transform', ants_xfm_list, 'transform')
        wf.connect(ants_xfm_list, 'transform_list', apply_warp, 'transforms')

        # parallelize the apply warp, if multiple CPUs, and it's a time
        # series!
        if int(num_cpus) > 1 and time_series:

            chunk_imports = ['import nibabel as nb']
            chunk = pe.Node(util.Function(input_names=['func_file',
                                                       'n_cpus'],
                                     output_names=['TR_ranges'],
                                     function=chunk_ts,
                                     imports=chunk_imports),
                            name=f'chunk_{wf_name}')

            chunk.inputs.n_cpus = int(num_cpus)
            wf.connect(inputNode, 'input_image', chunk, 'func_file')

            split_imports = ['import os', 'import subprocess']
            split = pe.Node(util.Function(input_names=['func_file',
                                                       'tr_ranges'],
                                     output_names=['split_funcs'],
                                     function=split_ts_chunks,
                                     imports=split_imports),
                            name=f'split_{wf_name}')

            wf.connect(inputNode, 'input_image', split, 'func_file')
            wf.connect(chunk, 'TR_ranges', split, 'tr_ranges')

            wf.connect(split, 'split_funcs', apply_warp, 'input_image')

            func_concat = pe.Node(interface=afni_utils.TCat(),
                                  name=f'func_concat_{wf_name}')
            func_concat.inputs.outputtype = 'NIFTI_GZ'

            wf.connect(apply_warp, 'output_image', func_concat, 'in_files')

            wf.connect(func_concat, 'out_file', outputNode, 'output_image')

        else:
            wf.connect(inputNode, 'input_image', apply_warp, 'input_image')
            wf.connect(apply_warp, 'output_image', outputNode, 'output_image')

    elif reg_tool == 'fsl':

        if multi_input:
            apply_warp = pe.MapNode(interface=fsl.FLIRT(),
                                    name=f'fsl_apply_warp',
                                    iterfield=['in_file'])
        else:
            apply_warp = pe.Node(interface=fsl.FLIRT(),
                                 name='fsl_apply_warp')

        apply_warp.inputs.apply_xfm = True

        interp_string = pe.Node(util.Function(input_names=['interpolation',
                                                           'reg_tool'],
                                              output_names=['interpolation'],
                                              function=interpolation_string),
                                name=f'interp_string')
        interp_string.inputs.reg_tool = reg_tool

        wf.connect(inputNode, 'interpolation', interp_string, 'interpolation')
        wf.connect(interp_string, 'interpolation', apply_warp, 'interp')

        # mni to t1
        wf.connect(inputNode, 'reference', apply_warp, 'reference')
        wf.connect(inputNode, 'transform', apply_warp, 'in_matrix_file')

        # parallelize the apply warp, if multiple CPUs, and it's a time
        # series!
        if int(num_cpus) > 1 and time_series:

            chunk_imports = ['import nibabel as nb']
            chunk = pe.Node(util.Function(input_names=['func_file',
                                                       'n_cpus'],
                                     output_names=['TR_ranges'],
                                     function=chunk_ts,
                                     imports=chunk_imports),
                            name=f'chunk_{wf_name}')

            chunk.inputs.n_cpus = int(num_cpus)

            wf.connect(inputNode, 'input_image', chunk, 'func_file')

            split_imports = ['import os', 'import subprocess']
            split = pe.Node(util.Function(input_names=['func_file',
                                                       'tr_ranges'],
                                     output_names=['split_funcs'],
                                     function=split_ts_chunks,
                                     imports=split_imports),
                            name=f'split_{wf_name}')

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



    apply_xfm = apply_transform(f'warp_{label}_to_template', reg_tool,
                                time_series=True, num_cpus=num_cpus,
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


def create_register_func_to_anat(phase_diff_distcor=False,
                                 name='register_func_to_anat'):

    """
    Registers a functional scan in native space to anatomical space using a
    linear transform and does not include bbregister.

    Parameters
    ----------
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
    linear_reg.inputs.cost = 'corratio'
    linear_reg.inputs.dof = 6

    # if fieldmap_distortion:

    def convert_pedir(pedir):
        # FSL Flirt requires pedir input encoded as an int
        conv_dct = {'x': 1, 'y': 2, 'z': 3, 'x-': -1, 'y-': -2, 'z-': -3,
                    'i': 1, 'j': 2, 'k': 3, 'i-': -1, 'j-': -2, 'k-': -3,
                    '-x': -1, '-i': -1, '-y': -2,
                    '-j': -2, '-z': -3, '-k': -3}
        if not isinstance(pedir, str):
            raise Exception("\n\nPhase-encoding direction must be a "
                            "string value.\n\nValue: {0}"
                            "\n\n".format(pedir))
        if pedir not in conv_dct.keys():
            raise Exception("\n\nInvalid phase-encoding direction "
                            "entered: {0}\n\n".format(pedir))
        return conv_dct[pedir]

    if phase_diff_distcor:
        register_func_to_anat.connect(
            inputNode_pedir, ('pedir', convert_pedir),
            linear_reg, 'pedir')
        register_func_to_anat.connect(inputspec, 'fieldmap',
                                      linear_reg, 'fieldmap')
        register_func_to_anat.connect(inputspec, 'fieldmapmask',
                                      linear_reg, 'fieldmapmask')
        register_func_to_anat.connect(inputNode_echospacing, 'echospacing',
                                      linear_reg, 'echospacing')

    register_func_to_anat.connect(inputspec, 'func', linear_reg, 'in_file')

    register_func_to_anat.connect(inputspec, 'anat', linear_reg, 'reference')

    register_func_to_anat.connect(inputspec, 'interp', linear_reg, 'interp')

    register_func_to_anat.connect(linear_reg, 'out_matrix_file',
                                  outputspec,
                                  'func_to_anat_linear_xfm_nobbreg')

    register_func_to_anat.connect(linear_reg, 'out_file',
                                  outputspec, 'anat_func_nobbreg')

    return register_func_to_anat


def create_bbregister_func_to_anat(phase_diff_distcor=False,
                                   name='bbregister_func_to_anat'):

    """
    Registers a functional scan in native space to structural.  This is
    meant to be used after create_nonlinear_register() has been run and
    relies on some of it's outputs.

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
        inputspec.anat_skull : string (nifti file)
            Corresponding full-head scan of subject
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
                                                       'anat_skull',
                                                       'linear_reg_matrix',
                                                       'anat_wm_segmentation',
                                                       'bbr_schedule',
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
    wm_bb_mask.inputs.op_string = '-thr 0.5 -bin'

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
        inputspec, 'anat_skull',
        bbreg_func_to_anat, 'reference')

    register_bbregister_func_to_anat.connect(
        inputspec, 'linear_reg_matrix',
        bbreg_func_to_anat, 'in_matrix_file')

    # if fieldmap_distortion:

    def convert_pedir(pedir):
        # FSL Flirt requires pedir input encoded as an int
        conv_dct = {'x': 1, 'y': 2, 'z': 3, 'x-': -1, 'y-': -2, 'z-': -3,
                    'i': 1, 'j': 2, 'k': 3, 'i-': -1, 'j-': -2, 'k-': -3,
                    '-x': -1, '-i': -1, '-y': -2,
                    '-j': -2, '-z': -3, '-k': -3}
        if not isinstance(pedir, str):
            raise Exception("\n\nPhase-encoding direction must be a "
                            "string value.\n\nValue: {0}"
                            "\n\n".format(pedir))
        if pedir not in conv_dct.keys():
            raise Exception("\n\nInvalid phase-encoding direction "
                            "entered: {0}\n\n".format(pedir))
        return conv_dct[pedir]

    if phase_diff_distcor:
        register_bbregister_func_to_anat.connect(
            inputNode_pedir, ('pedir', convert_pedir),
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
                                                     'reference_mask',
                                                     'moving_mask',
                                                     'ants_para',
                                                     'fixed_image_mask',
                                                     'interp'],
                                        output_names=['warp_list',
                                                      'warped_image'],
                                        function=hardcoded_reg,
                                        imports=reg_imports),
                name='calc_ants_warp')

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

    if reg_ants_skull == 1 and not reg_ants_skull == 0:
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
                               symmetric=False):

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

    if opt == 'FSL' or opt == 'FSL-linear':

        flirt_reg_anat_mni = create_fsl_flirt_linear_reg(
            f'anat_mni_flirt_register{symm}'
        )

        # Input registration parameters
        wf.connect(inputNode, 'interpolation', flirt_reg_anat_mni, 'interp')

        wf.connect(inputNode, 'input_brain',
                   flirt_reg_anat_mni, 'inputspec.input_brain')

        wf.connect(inputNode, 'reference_brain', flirt_reg_anat_mni,
                   'inputspec.reference_brain')

        outputs = {
            f'space-{sym}template_desc-brain_{orig}': (
                flirt_reg_anat_mni, 'outputspec.output_brain'),
            f'from-{orig}_to-{sym}template_mode-image_desc-linear_xfm': (
                flirt_reg_anat_mni, 'outputspec.linear_xfm'),
            f'from-{sym}template_to-{orig}_mode-image_desc-linear_xfm': (
                flirt_reg_anat_mni, 'outputspec.invlinear_xfm'),
            f'from-{orig}_to-{sym}template_mode-image_xfm': (
                flirt_reg_anat_mni, 'outputspec.linear_xfm')
        }

    if opt == 'FSL':
        fnirt_reg_anat_mni = create_fsl_fnirt_nonlinear_reg(
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

        write_composite_xfm = pe.Node(interface=fsl.ConvertWarp(),
                                      name=f'combine_fsl_warps{symm}')

        wf.connect(inputNode, 'reference_head',
                   write_composite_xfm, 'reference')

        wf.connect(flirt_reg_anat_mni, 'outputspec.linear_xfm',
                   write_composite_xfm, 'premat')

        wf.connect(fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm',
                   write_composite_xfm, 'warp1')

        # NOTE: this is an UPDATE because of the opt block above
        added_outputs = {
            f'space-{sym}template_desc-brain_{orig}': (
                fnirt_reg_anat_mni, 'outputspec.output_brain'),
            f'from-{orig}_to-{sym}template_mode-image_desc-nonlinear_xfm': (
                fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm'),
            f'from-{orig}_to-{sym}template_mode-image_xfm': (
                write_composite_xfm, 'out_file')
        }
        outputs.update(added_outputs)

    return (wf, outputs)


def ANTs_registration_connector(wf_name, cfg, params, orig='T1w',
                                symmetric=False):

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
        mem_gb=1.5)
    write_composite_linear_xfm.inputs.print_out_composite_warp_file = True
    write_composite_linear_xfm.inputs.output_image = \
        "from-T1w_to-template_mode-image_desc-linear_xfm.nii.gz"

    wf.connect(inputNode, 'input_brain',
               write_composite_linear_xfm, 'input_image')

    wf.connect(inputNode, 'reference_brain',
               write_composite_linear_xfm, 'reference_image')

    wf.connect(inputNode, 'interpolation',
               write_composite_linear_xfm, 'interpolation')

    write_composite_linear_xfm.inputs.input_image_type = 0
    write_composite_linear_xfm.inputs.dimension = 3

    collect_transforms = pe.Node(util.Merge(3),
                                 name=f'collect_transforms{symm}')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
               collect_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_transforms, 'in3')

    # check transform list to exclude Nonetype (missing) init/rig/affine
    check_transform = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['checked_transform_list',
                                    'list_length'],
                      function=check_transforms),
        name=f'check_transforms')

    wf.connect(collect_transforms, 'out', check_transform, 'transform_list')

    wf.connect(check_transform, 'checked_transform_list',
               write_composite_linear_xfm, 'transforms')

    # combine the linear xfm's into one - makes it easier downstream
    write_composite_invlinear_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_invlinear{symm}_xfm',
        mem_gb=1.5)
    write_composite_invlinear_xfm.inputs.print_out_composite_warp_file = True
    write_composite_invlinear_xfm.inputs.output_image = \
        "from-template_to-T1w_mode-image_desc-linear_xfm.nii.gz"

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

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_inv_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_inv_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
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
        "from-T1w_to-template_mode-image_xfm.nii.gz"

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

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
               collect_all_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_all_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_all_transforms, 'in3')

    wf.connect(ants_reg_anat_mni, 'outputspec.warp_field',
               collect_all_transforms, 'in4')

    wf.connect(collect_all_transforms, 'out',
               write_composite_xfm, 'transforms')

    # combine ALL xfm's into one - makes it easier downstream
    write_composite_inv_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_inv_{symm}xfm',
        mem_gb=1.5)
    write_composite_inv_xfm.inputs.print_out_composite_warp_file = True
    write_composite_inv_xfm.inputs.output_image = \
        "from-template_to-T1w_mode-image_xfm.nii.gz"

    wf.connect(inputNode, 'reference_brain',
               write_composite_inv_xfm, 'input_image')

    wf.connect(inputNode, 'input_brain',
               write_composite_inv_xfm, 'reference_image')

    wf.connect(inputNode, 'interpolation',
               write_composite_inv_xfm, 'interpolation')

    write_composite_inv_xfm.inputs.input_image_type = 0
    write_composite_inv_xfm.inputs.dimension = 3
    write_composite_inv_xfm.inputs.invert_transform_flags = [False, True,
                                                             True, True]

    collect_all_inv_transforms = pe.Node(util.Merge(4),
                                         name=f'collect_all_inv_transforms'
                                         f'{symm}')

    wf.connect(ants_reg_anat_mni, 'outputspec.inverse_warp_field',
               collect_all_inv_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_all_inv_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_all_inv_transforms, 'in3')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
               collect_all_inv_transforms, 'in4')

    wf.connect(collect_all_inv_transforms, 'out',
               write_composite_inv_xfm, 'transforms')

    outputs = {
        f'space-{sym}template_desc-brain_{orig}': (
            ants_reg_anat_mni, 'outputspec.normalized_output_brain'),
        f'from-{orig}_to-{sym}template_mode-image_xfm': (
            write_composite_xfm, 'output_image'),
        f'from-{sym}template_to-{orig}_mode-image_xfm': (
            write_composite_inv_xfm, 'output_image'),
        f'from-{orig}_to-{sym}template_mode-image_desc-linear_xfm': (
            write_composite_linear_xfm, 'output_image'),
        f'from-{sym}template_to-{orig}_mode-image_desc-linear_xfm': (
            write_composite_invlinear_xfm, 'output_image'),
        f'from-{orig}_to-{sym}template_mode-image_desc-nonlinear_xfm': (
            ants_reg_anat_mni, 'outputspec.warp_field'),
        f'from-{sym}template_to-{orig}_mode-image_desc-nonlinear_xfm': (
            ants_reg_anat_mni, 'outputspec.inverse_warp_field')
    }

    return (wf, outputs)


def bold_to_T1template_xfm_connector(wf_name, cfg, reg_tool, symmetric=False):

    wf = pe.Workflow(name=wf_name)

    inputNode = pe.Node(
        util.IdentityInterface(fields=['input_brain',
                                       'mean_bold',
                                       'coreg_xfm',
                                       'T1w_brain_template_funcreg',
                                       'T1w_to_template_xfm',
                                       'template_to_T1w_xfm']),
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

        wf.connect(inputNode, 'input_brain',
                   write_composite_xfm, 'input_image')

        wf.connect(inputNode, 'T1w_brain_template_funcreg',
                   write_composite_xfm, 'reference_image')

        write_composite_xfm.inputs.input_image_type = 0
        write_composite_xfm.inputs.dimension = 3
        write_composite_xfm.inputs.interpolation = \
            cfg.registration_workflows['anatomical_registration'][
                'registration']['ANTs']['interpolation']

        collect_all_transforms = pe.Node(util.Merge(2),
                                         name=f'collect_all_transforms')

        wf.connect(change_transform, 'updated_affine_file',
                   collect_all_transforms, 'in1')

        wf.connect(inputNode, 'T1w_to_template_xfm',
                   collect_all_transforms, 'in2')

        wf.connect(collect_all_transforms, 'out',
                   write_composite_xfm, 'transforms')

        write_composite_inv_xfm = pe.Node(
            interface=ants.ApplyTransforms(),
            name=f'write_composite_inv_xfm',
            mem_gb=1.5)
        write_composite_inv_xfm.inputs.print_out_composite_warp_file = True
        write_composite_inv_xfm.inputs.invert_transform_flags = [False, True]
        write_composite_inv_xfm.inputs.output_image = \
            f"from-{sym}template_to-bold_mode-image_xfm.nii.gz"

        wf.connect(inputNode, 'T1w_brain_template_funcreg',
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

        wf.connect(inputNode, 'template_to_T1w_xfm',
                   collect_inv_transforms, 'in1')

        wf.connect(change_transform, 'updated_affine_file',
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

        wf.connect(inputNode, 'T1w_brain_template',
                   write_composite_xfm, 'reference')

        wf.connect(inputNode, 'coreg_xfm', write_composite_xfm, 'premat')

        wf.connect(inputNode, 'T1w_to_template_xfm',
                   write_composite_xfm, 'warp1')

        outputs = {
            f'from-bold_to-{sym}template_mode-image_xfm':
                (write_composite_xfm, 'out_file'),
        }

    return (wf, outputs)


def register_FSL_anat_to_template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "register_FSL_anat_to_template",
     "config": ["registration_workflows", "anatomical_registration",
                "registration"],
     "switch": "run",
     "option_key": "using",
     "option_val": ["FSL", "FSL-linear"],
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "desc-brain_T1w",
                "T1w_template",
                "T1w_brain_template",
                "template_ref_mask"],
     "outputs": ["space-template_desc-brain_T1w",
                 "from-T1w_to-template_mode-image_desc-linear_xfm",
                 "from-template_to-T1w_mode-image_desc-linear_xfm",
                 "from-T1w_to-template_mode-image_desc-nonlinear_xfm",
                 "from-T1w_to-template_mode-image_xfm"]}
    '''

    fsl, outputs = FSL_registration_connector('register_FSL_anat_to_'
                                              f'template_{pipe_idx}', cfg,
                                              'T1w', opt)

    fsl.inputs.inputspec.interpolation = cfg['registration_workflows'][
        'anatomical_registration']['FSL-FNIRT']['interpolation']

    fsl.inputs.inputspec.fnirt_config = cfg['registration_workflows'][
        'registration']['FSL-FNIRT']['fnirt_config']

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, fsl, 'inputspec.input_brain')

    node, out = strat_pool.get_data('T1w_brain_template')
    wf.connect(node, out, fsl, 'inputspec.reference_brain')

    node, out = strat_pool.get_data(["desc-preproc_T1w",
                                     "desc-reorient_T1w", "T1w"])
    wf.connect(node, out, fsl, 'inputspec.input_head')

    node, out = strat_pool.get_data('T1w_template')
    wf.connect(node, out, fsl, 'inputspec.reference_head')

    node, out = strat_pool.get_data('ref_mask')
    wf.connect(node, out, fsl, 'inputspec.reference_mask')

    return (wf, outputs)


def register_symmetric_FSL_anat_to_template(wf, cfg, strat_pool, pipe_num,
                                            opt=None):
    '''
    {"name": "register_symmetric_FSL_anat_to_template",
     "config": ["registration_workflows", "anatomical_registration",
                "registration"],
     "switch": "run",
     "option_key": "using",
     "option_val": ["FSL", "FSL-linear"],
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "desc-brain_T1w",
                "T1w_template",
                "T1w_brain_template",
                "template_ref_mask"],
     "outputs": ["space-symtemplate_desc-brain_T1w",
                 "from-T1w_to-symtemplate_mode-image_desc-linear_xfm",
                 "from-symtemplate_to-T1w_mode-image_desc-linear_xfm",
                 "from-T1w_to-symtemplate_mode-image_desc-nonlinear_xfm",
                 "from-T1w_to-symtemplate_mode-image_xfm"]}
    '''

    fsl, outputs = FSL_registration_connector('register_FSL_anat_to_'
                                              f'template_symmetric_'
                                              f'{pipe_idx}', cfg, 'T1w', opt,
                                              symmetric=True)

    fsl.inputs.inputspec.interpolation = cfg['registration_workflows'][
        'anatomical_registration']['FSL-FNIRT']['interpolation']

    fsl.inputs.inputspec.fnirt_config = cfg['registration_workflows'][
        'registration']['FSL-FNIRT']['fnirt_config']

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, fsl, 'inputspec.input_brain')

    node, out = strat_pool.get_data('T1w_brain_template_symmetric')
    wf.connect(node, out, fsl, 'inputspec.reference_brain')

    node, out = strat_pool.get_data(["desc-preproc_T1w",
                                     "desc-reorient_T1w", "T1w"])
    wf.connect(node, out, fsl, 'inputspec.input_head')

    node, out = strat_pool.get_data('T1w_template_symmetric')
    wf.connect(node, out, fsl, 'inputspec.reference_head')

    node, out = strat_pool.get_data('dilated_symmetric_brain_mask')
    wf.connect(node, out, fsl, 'inputspec.reference_mask')

    return (wf, outputs)


def register_FSL_EPI_to_template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Directly register the mean functional to an EPI template. No T1w
    involved.

    Node Block:
    {"name": "register_FSL_EPI_to_template",
     "config": ["registration_workflows", "functional_registration",
                "EPI_registration"],
     "switch": ["run"],
     "option_key": "using",
     "option_val": ["FSL", "FSL-linear"],
     "inputs": ["bold_coreg_input",
                "space-bold_desc-brain_mask",
                "EPI_template",
                "EPI_template_mask"],
     "outputs": ["space-template_desc-brain_bold",
                 "from-bold_to-template_mode-image_desc-linear_xfm",
                 "from-template_to-bold_mode-image_desc-linear_xfm",
                 "from-bold_to-template_mode-image_desc-nonlinear_xfm",
                 "from-template_to-bold_mode-image_desc-nonlinear_xfm",
                 "from-bold_to-template_mode-image_xfm"]}
    '''

    fsl, outputs = FSL_registration_connector('register_FSL_EPI_to_'
                                              f'template_{pipe_idx}', cfg,
                                              opt, 'bold')

    fsl.inputs.inputspec.interpolation = cfg['registration_workflows'][
        'functional_registration']['EPI_registration']['FSL-FNIRT'][
        'interpolation']

    fsl.inputs.inputspec.fnirt_config = cfg['registration_workflows'][
        'functional_registration']['EPI_registration']['FSL-FNIRT'][
        'fnirt_config']

    node, out = strat_pool.get_data('bold_coreg_input')
    wf.connect(node, out, fsl, 'inputspec.input_brain')

    node, out = strat_pool.get_data('EPI_template')
    wf.connect(node, out, fsl, 'inputspec.reference_brain')

    node, out = strat_pool.get_data('bold_coreg_input')
    wf.connect(node, out, fsl, 'inputspec.input_head')

    node, out = strat_pool.get_data('EPI_template')
    wf.connect(node, out, fsl, 'inputspec.reference_head')

    node, out = strat_pool.get_data('EPI_template_mask')
    wf.connect(node, out, fsl, 'inputspec.reference_mask')

    return (wf, outputs)


def register_ANTs_anat_to_template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "register_ANTs_anat_to_template",
     "config": ["registration_workflows", "anatomical_registration"],
     "switch": ["run"],
     "option_key": ["registration", "using"],
     "option_val": "ANTS",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "desc-brain_T1w",
                "space-T1w_desc-brain_mask",
                "T1w_template",
                "T1w_brain_template",
                "T1w_brain_template_mask",
                "label-lesion_mask"],
     "outputs": ["space-template_desc-brain_T1w",
                 "from-T1w_to-template_mode-image_desc-linear_xfm",
                 "from-template_to-T1w_mode-image_desc-linear_xfm",
                 "from-T1w_to-template_mode-image_desc-nonlinear_xfm",
                 "from-template_to-T1w_mode-image_desc-nonlinear_xfm",
                 "from-T1w_to-template_mode-image_xfm",
                 "from-template_to-T1w_mode-image_xfm"]}
    '''

    params = cfg.registration_workflows['anatomical_registration'][
        'registration']['ANTs']['T1_registration']

    ants, outputs = ANTs_registration_connector('ANTS_T1_to_template', cfg,
                                                params, 'T1w')

    ants.inputs.inputspec.interpolation = cfg.registration_workflows[
        'anatomical_registration']['registration']['ANTs']['interpolation']

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, ants, 'inputspec.input_brain')

    node, out = strat_pool.get_data('T1w_brain_template')
    wf.connect(node, out, ants, 'inputspec.reference_brain')

    node, out = strat_pool.get_data(["desc-preproc_T1w", "desc-reorient_T1w",
                                     "T1w"])
    wf.connect(node, out, ants, 'inputspec.input_head')

    node, out = strat_pool.get_data(f'T1w_template')
    wf.connect(node, out, ants, 'inputspec.reference_head')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, ants, 'inputspec.input_mask')

    node, out = strat_pool.get_data('T1w_brain_template_mask')
    wf.connect(node, out, ants, 'inputspec.reference_mask')

    if strat_pool.check_rpool('label-lesion_mask'):
        node, out = strat_pool.get_data('label-lesion_mask')
        wf.connect(node, out, ants, 'inputspec.lesion_mask')

    return (wf, outputs)


def register_symmetric_ANTs_anat_to_template(wf, cfg, strat_pool, pipe_num,
                                             opt=None):
    '''
    {"name": "register_symmetric_ANTs_anat_to_template",
     "config": ["registration_workflows", "anatomical_registration"],
     "switch": ["run"],
     "option_key": ["registration", "using"],
     "option_val": "ANTS",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "desc-brain_T1w",
                "space-T1w_desc-brain_mask",
                "T1w_template_symmetric",
                "T1w_brain_template_symmetric",
                "dilated_symmetric_brain_mask",
                "label-lesion_mask"],
     "outputs": ["space-symtemplate_desc-brain_T1w",
                 "from-T1w_to-symtemplate_mode-image_desc-linear_xfm",
                 "from-symtemplate_to-T1w_mode-image_desc-linear_xfm",
                 "from-T1w_to-symtemplate_mode-image_desc-nonlinear_xfm",
                 "from-symtemplate_to-T1w_mode-image_desc-nonlinear_xfm",
                 "from-T1w_to-symtemplate_mode-image_xfm",
                 "from-symtemplate_to-T1w_mode-image_xfm"]}
    '''

    params = cfg.registration_workflows['anatomical_registration'][
        'registration']['ANTs']['T1_registration']

    ants, outputs = ANTs_registration_connector('ANTS_T1_to_template_'
                                                'symmetric', cfg, params,
                                                'T1w', True)

    ants.inputs.inputspec.interpolation = cfg.registration_workflows[
        'anatomical_registration']['registration']['ANTs']['interpolation']

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, ants, 'inputspec.input_brain')

    node, out = strat_pool.get_data('T1w_brain_template_symmetric')
    wf.connect(node, out, ants, 'inputspec.reference_brain')

    node, out = strat_pool.get_data(["desc-preproc_T1w", "desc-reorient_T1w",
                                     "T1w"])
    wf.connect(node, out, ants, 'inputspec.input_head')

    node, out = strat_pool.get_data('T1w_template_symmetric')
    wf.connect(node, out, ants, 'inputspec.reference_head')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, ants, 'inputspec.input_mask')

    node, out = strat_pool.get_data('dilated_symmetric_brain_mask')
    wf.connect(node, out, ants, 'inputspec.reference_mask')

    if strat_pool.check_rpool('label-lesion_mask'):
        node, out = strat_pool.get_data('label-lesion_mask')
        wf.connect(node, out, ants, 'inputspec.lesion_mask')

    return (wf, outputs)


def register_ANTs_EPI_to_template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Directly register the mean functional to an EPI template. No T1w
    involved.

    Node Block:
    {"name": "register_ANTs_EPI_to_template",
     "config": ["registration_workflows", "functional_registration",
                "EPI_registration"],
     "switch": ["run"],
     "option_key": "using",
     "option_val": "ANTS",
     "inputs": ["bold_coreg_input",
                "space-bold_desc-brain_mask",
                "EPI_template",
                "EPI_template_mask"],
     "outputs": ["space-template_desc-brain_bold",
                 "from-bold_to-template_mode-image_desc-linear_xfm",
                 "from-template_to-bold_mode-image_desc-linear_xfm",
                 "from-bold_to-template_mode-image_desc-nonlinear_xfm",
                 "from-template_to-bold_mode-image_desc-nonlinear_xfm",
                 "from-bold_to-template_mode-image_xfm",
                 "from-template_to-bold_mode-image_xfm"]}
    '''

    params = cfg.registration_workflows['functional_registration'][
        'EPI_registration']['ANTs']['parameters']

    ants, outputs = ANTs_registration_connector('ANTS_T1_to_EPI_template',
                                                cfg, params, 'bold')

    ants.inputs.inputspec.interpolation = cfg.registration_workflows[
        'functional_registration']['EPI_registration']['ANTs'][
        'interpolation']

    node, out = strat_pool.get_data('bold_coreg_input')
    wf.connect(node, out, ants, 'inputspec.input_brain')

    node, out = strat_pool.get_data('EPI_template')
    wf.connect(node, out, ants, 'inputspec.reference_brain')

    node, out = strat_pool.get_data('bold_coreg_input')
    wf.connect(node, out, ants, 'inputspec.input_head')

    node, out = strat_pool.get_data('EPI_template')
    wf.connect(node, out, ants, 'inputspec.reference_head')

    node, out = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out, ants, 'inputspec.input_mask')

    node, out = strat_pool.get_data('EPI_template_mask')
    wf.connect(node, out, ants, 'inputspec.reference_mask')

    return (wf, outputs)


def coregistration_prep_vol(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "coregistration_prep_vol",
     "config": ["registration_workflows", "functional_registration",
                "coregistration", "func_input_prep"],
     "switch": "None",
     "option_key": "input",
     "option_val": "Selected_Functional_Volume",
     "inputs": ["desc-brain_bold"],
     "outputs": ["bold_coreg_input"]}
    '''

    get_func_volume = pe.Node(interface=afni.Calc(),
                              name=f'get_func_volume_{pipe_num}')

    get_func_volume.inputs.set(
        expr='a',
        single_idx=cfg.func_reg_input_volume,
        outputtype='NIFTI_GZ'
    )
    node, out = strat_pool.get_data("desc-brain_bold")
    wf.connect(node, out, get_func_volume, 'in_file_a')

    coreg_input = (get_func_volume, 'out_file')

    outputs = {
        'bold_coreg_input': coreg_input
    }

    return (wf, outputs)


def coregistration_prep_mean(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "coregistration_prep_mean",
     "config": ["registration_workflows", "functional_registration",
                "coregistration", "func_input_prep"],
     "switch": "None",
     "option_key": "input",
     "option_val": "Mean_Functional",
     "inputs": ["desc-mean_bold"],
     "outputs": ["bold_coreg_input"]}
    '''

    coreg_input = strat_pool.get_data("desc-mean_bold")

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
        'bold_coreg_input': coreg_input
    }

    return (wf, outputs)


def coregistration(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "coregistration",
     "config": ["registration_workflows", "functional_registration",
                "coregistration"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["bold_coreg_input",
                "desc-brain_T1w",
                "diff_phase_dwell",
                "diff_phase_pedir",
                "despiked_fieldmap",
                "fieldmap_mask"],
     "outputs": ["space-T1w_desc-mean_bold",
                 "from-bold_to-T1w_mode-image_desc-linear_xfm"]}
    '''

    diff_complete = False
    if strat_pool.check_rpool("despiked_fieldmap") and \
            strat_pool.check_rpool("fieldmap_mask"):
        diff_complete = True

    # if field map-based distortion correction is on, but BBR is off,
    # send in the distortion correction files here
    func_to_anat = create_register_func_to_anat(diff_complete,
                                                f'func_to_anat_FLIRT_'
                                                f'{pipe_num}')
    func_to_anat.inputs.inputspec.interp = 'trilinear'

    node, out = strat_pool.get_data('bold_coreg_input')
    wf.connect(node, out, func_to_anat, 'inputspec.func')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, func_to_anat, 'inputspec.anat')

    if diff_complete:
        node, out = strat_pool.get_data('diff_phase_dwell')
        wf.connect(node, out, func_to_anat, 'echospacing_input.echospacing')

        node, out = strat_pool.get_data('diff_phase_pedir')
        wf.connect(node, out, func_to_anat, 'pedir_input.pedir')

        node, out = strat_pool.get_data("despiked_fieldmap")
        wf.connect(node, out, func_to_anat, 'inputspec.fieldmap')

        node, out = strat_pool.get_data("fieldmap_mask")
        wf.connect(node, out, func_to_anat, 'inputspec.fieldmapmask')

    outputs = {
        'space-T1w_desc-mean_bold':
            (func_to_anat, 'outputspec.anat_func_nobbreg'),
        'from-bold_to-T1w_mode-image_desc-linear_xfm':
            (func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg')
    }

    return (wf, outputs)


def bbr_coregistration(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "bbr_coregistration",
     "config": ["registration_workflows", "functional_registration",
                "coregistration", "boundary_based_registration"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["bold_coreg_input",
                "T1w",
                "from-bold_to-T1w_mode-image_desc-linear_xfm",
                "label-CSF_probseg",
                "diff_phase_dwell",
                "diff_phase_pedir",
                "despiked_fieldmap",
                "fieldmap_mask"],
     "outputs": ["space-T1w_desc-mean_bold",
                 "from-bold_to-T1w_mode-image_desc-linear_xfm"]}
    '''

    diff_complete = False
    if strat_pool.check_rpool("despiked_fieldmap") and \
            strat_pool.check_rpool("fieldmap_mask"):
        diff_complete = True

    func_to_anat_bbreg = create_bbregister_func_to_anat(diff_complete,
                                                        f'func_to_anat_'
                                                        f'bbreg_'
                                                        f'{pipe_num}')
    func_to_anat_bbreg.inputs.inputspec.bbr_schedule = \
        cfg.registration_workflows['functional_registration'][
            'coregistration']['boundary_based_registration'][
            'bbr_schedule']

    node, out = strat_pool.get_data('bold_coreg_input')
    wf.connect(node, out, func_to_anat_bbreg, 'inputspec.func')

    node, out = strat_pool.get_data('T1w')
    wf.connect(node, out, func_to_anat_bbreg, 'inputspec.anat_skull')

    node, out = strat_pool.get_data(
        'from-bold_to-T1w_mode-image_desc-linear_xfm')
    wf.connect(node, out,
               func_to_anat_bbreg, 'inputspec.linear_reg_matrix')

    node, out = strat_pool.get_data('label-CSF_probseg')
    wf.connect(node, out,
               func_to_anat_bbreg, 'inputspec.anat_wm_segmentation')

    if diff_complete:
        node, out = strat_pool.get_data('diff_phase_dwell')
        wf.connect(node, out,
                   func_to_anat_bbreg, 'echospacing_input.echospacing')

        node, out = strat_pool.get_data('diff_phase_pedir')
        wf.connect(node, out, func_to_anat_bbreg, 'pedir_input.pedir')

        node, out = strat_pool.get_data("despiked_fieldmap")
        wf.connect(node, out, func_to_anat_bbreg, 'inputspec.fieldmap')

        node, out = strat_pool.get_data("fieldmap_mask")
        wf.connect(node, out,
                   func_to_anat_bbreg, 'inputspec.fieldmapmask')

    outputs = {
        'space-T1w_desc-mean_bold':
            (func_to_anat_bbreg, 'outputspec.anat_func'),
        'from-bold_to-T1w_mode-image_desc-linear_xfm':
            (func_to_anat_bbreg, 'outputspec.func_to_anat_linear_xfm')
    }

    return (wf, outputs)


def create_func_to_T1template_xfm(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Condense the BOLD-to-T1 coregistration transform and the T1-to-template
    transform into one transform matrix.

    Node Block:
    {"name": "create_func_to_T1template_xfm",
     "config": ["registration_workflows", "functional_registration",
                "func_registration_to_template"],
     "switch": ["run"],
     "option_key": ["target_template", "using"],
     "option_val": "T1_template",
     "inputs": ["from-bold_to-T1w_mode-image_desc-linear_xfm",
                "from-T1w_to-template_mode-image_xfm",
                "from-template_to-T1w_mode-image_xfm",
                "desc-mean_bold",
                "desc-brain_T1w",
                "T1w_brain_template_funcreg"],
     "outputs": ["from-bold_to-template_mode-image_xfm",
                 "from-template_to-bold_mode-image_xfm"]}
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

    node, out = strat_pool.get_data('desc-mean_bold')
    wf.connect(node, out, xfm, 'inputspec.mean_bold')

    node, out = strat_pool.get_data('T1w_brain_template_funcreg')
    wf.connect(node, out, xfm, 'inputspec.T1w_brain_template_funcreg')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, xfm, 'inputspec.T1w_to_template_xfm')

    # FNIRT pipelines don't have an inverse nonlinear warp, make optional
    if strat_pool.check_rpool('from-template_to-T1w_mode-image_xfm'):
        node, out = strat_pool.get_data('from-template_to-T1w_mode-image_xfm')
        wf.connect(node, out, xfm, 'inputspec.template_to_T1w_xfm')

    return (wf, outputs)


def create_func_to_T1template_symmetric_xfm(wf, cfg, strat_pool, pipe_num,
                                            opt=None):
    '''Condense the BOLD-to-T1 coregistration transform and the T1-to-
    symmetric-template transform into one transform matrix.

    Node Block:
    {"name": "create_func_to_T1template_symmetric_xfm",
     "config": ["registration_workflows", "functional_registration",
                "func_registration_to_template"],
     "switch": ["run"],
     "option_key": ["target_template", "using"],
     "option_val": "T1_template",
     "inputs": ["from-bold_to-T1w_mode-image_desc-linear_xfm",
                "from-T1w_to-symtemplate_mode-image_xfm",
                "from-symtemplate_to-T1w_mode-image_xfm",
                "desc-mean_bold",
                "desc-brain_T1w",
                "T1w_brain_template_symmetric_deriv"],
     "outputs": ["from-bold_to-symtemplate_mode-image_xfm",
                 "from-symtemplate_to-bold_mode-image_xfm"]}
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

    node, out = strat_pool.get_data('desc-mean_bold')
    wf.connect(node, out, xfm, 'inputspec.mean_bold')

    node, out = strat_pool.get_data('T1w_brain_template_symmetric_deriv')
    wf.connect(node, out, xfm, 'inputspec.T1w_brain_template_funcreg')

    node, out = strat_pool.get_data('from-T1w_to-symtemplate_mode-image_xfm')
    wf.connect(node, out, xfm, 'inputspec.T1w_to_template_xfm')

    # FNIRT pipelines don't have an inverse nonlinear warp, make optional
    if strat_pool.check_rpool('from-symtemplate_to-T1w_mode-image_xfm'):
        node, out = \
            strat_pool.get_data('from-symtemplate_to-T1w_mode-image_xfm')
        wf.connect(node, out, xfm, 'inputspec.template_to_T1w_xfm')

    return (wf, outputs)


def warp_timeseries_to_T1template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    Node Block:
    {"name": "transform_timeseries_to_T1template",
     "config": ["registration_workflows", "functional_registration",
                "func_registration_to_template"],
     "switch": ["run"],
     "option_key": ["target_template", "using"],
     "option_val": "T1_template",
     "inputs": [["desc-cleaned_bold", "desc-preproc_bold",
                 "desc-reorient_bold", "bold"],
                "from-bold_to-template_mode-image_xfm",
                "T1w_brain_template_funcreg"],
     "outputs": ["space-template_desc-cleaned_bold",
                 "space-template_desc-preproc_bold",
                 "space-template_desc-reorient_bold",
                 "space-template_bold"]}
    '''

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

    connect, resource = strat_pool.get_data(["desc-cleaned_bold",
                                             "desc-preproc_bold",
                                             "desc-reorient_bold",
                                             "bold"],
                                            report_fetched=True)

    node, out = connect
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("T1w_brain_template_funcreg")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-bold_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        f'space-template_{resource}': (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


def warp_bold_mask_to_T1template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    Node Block:
    {"name": "transform_bold_mask_to_T1template",
     "config": ["registration_workflows", "functional_registration",
                "func_registration_to_template"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["space-bold_desc-brain_mask",
                "from-bold_to-template_mode-image_xfm",
                "T1w_brain_template_funcreg"],
     "outputs": ["space-template_desc-bold_mask"]}
    '''

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

    node, out = strat_pool.get_data("space-bold_desc-brain_mask")
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("T1w_brain_template_funcreg")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-bold_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        'space-template_desc-bold_mask':
            (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


def warp_deriv_mask_to_T1template(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Transform the BOLD mask to template space and to the resolution set for
    the derivative outputs.

    Node Block:
    {"name": "transform_bold_mask_to_T1template",
     "config": ["registration_workflows", "functional_registration",
                "func_registration_to_template"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["space-bold_desc-brain_mask",
                "from-bold_to-template_mode-image_xfm",
                "T1w_brain_template_deriv"],
     "outputs": ["space-template_res-derivative_desc-bold_mask"]}
    '''

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

    node, out = strat_pool.get_data("space-bold_desc-brain_mask")
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("T1w_brain_template_deriv")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-bold_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        f'space-template_res-derivative_desc-bold_mask':
            (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


def warp_timeseries_to_EPItemplate(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    Node Block:
    {"name": "transform_timeseries_to_EPItemplate",
     "config": ["registration_workflows", "functional_registration",
                "func_registration_to_template"],
     "switch": ["run"],
     "option_key": ["target_template", "using"],
     "option_val": "EPI_template",
     "inputs": [["desc-cleaned_bold", "desc-preproc_bold",
                 "desc-reorient_bold", "bold"],
                "from-bold_to-template_mode-image_xfm",
                "EPI_template"],
     "outputs": ["space-template_desc-cleaned_bold",
                 "space-template_desc-preproc_bold",
                 "space-template_desc-reorient_bold",
                 "space-template_bold"]}
    '''

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-bold_to-template_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    apply_xfm = apply_transform('warp_ts_to_EPItemplate', reg_tool,
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

    connect, resource = strat_pool.get_data(["desc-cleaned_bold",
                                             "desc-preproc_bold",
                                             "desc-reorient_bold",
                                             "bold"],
                                            report_fetched=True)
    node, out = connect
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("EPI_template")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-bold_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        f'space-template_{resource}': (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)
