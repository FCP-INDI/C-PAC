import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants

from CPAC.anat_preproc.lesion_preproc import create_lesion_preproc

from CPAC.registration.utils import seperate_warps_list, \
    check_transforms, \
    generate_inverse_transform_flags, \
    single_ants_xfm_to_list, \
    interpolation_string, \
    hardcoded_reg


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


def create_register_func_to_mni(name='register_func_to_mni'):
    """
    Registers a functional scan in native space to MNI standard space.
    This is meant to be used after create_nonlinear_register() has been
    run and relies on some of it's outputs.

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    register_func_to_mni : nipype.pipeline.engine.Workflow

    Notes
    -----

    Workflow Inputs::

        inputspec.func : string (nifti file)
            Input functional scan to be registered to MNI space
        inputspec.mni : string (nifti file)
            Reference MNI file
        inputspec.anat : string (nifti file)
            Corresponding anatomical scan of subject
        inputspec.interp : string
            Type of interpolation to use
            ('trilinear' or 'nearestneighbour' or 'sinc')
        inputspec.anat_to_mni_nonlinear_xfm : string (warp file)
            Corresponding anatomical native space to MNI warp file
        inputspec.anat_to_mni_linear_xfm : string (mat file)
            Corresponding anatomical native space to MNI mat file

    Workflow Outputs::

        outputspec.func_to_anat_linear_xfm : string (mat file)
            Affine transformation from functional to anatomical native space
        outputspec.func_to_mni_linear_xfm : string (mat file)
            Affine transformation from functional to MNI space
        outputspec.mni_to_func_linear_xfm : string (mat file)
            Affine transformation from MNI to functional space
        outputspec.mni_func : string (nifti file)
            Functional scan registered to MNI standard space

    .. exec::
        from CPAC.registration import create_register_func_to_mni
        wf = create_register_func_to_mni()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/register_func_to_mni.dot'
        )

    Workflow Graph:

    .. image:: ../../images/generated/register_func_to_mni.png
        :width: 500

    Detailed Workflow Graph:
    .. image:: ../images/register_func_to_mni_detailed.dot.png
        :width: 500
    """
    register_func_to_mni = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(
        fields=['func', 'mni', 'anat', 'interp', 'anat_to_mni_nonlinear_xfm',
                'anat_to_mni_linear_xfm']),
        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(
        fields=['func_to_anat_linear_xfm', 'func_to_mni_linear_xfm',
                'mni_to_func_linear_xfm', 'mni_func']),
        name='outputspec')

    linear_reg = pe.Node(interface=fsl.FLIRT(),
                         name='linear_func_to_anat')
    linear_reg.inputs.cost = 'corratio'
    linear_reg.inputs.dof = 6

    mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                       name='mni_warp')

    mni_affine = pe.Node(interface=fsl.ConvertXFM(),
                         name='mni_affine')
    mni_affine.inputs.concat_xfm = True
    register_func_to_mni.connect(linear_reg, 'out_matrix_file',
                                 mni_affine, 'in_file2')
    register_func_to_mni.connect(inputspec, 'anat_to_mni_linear_xfm',
                                 mni_affine, 'in_file')
    register_func_to_mni.connect(mni_affine, 'out_file',
                                 outputspec, 'func_to_mni_linear_xfm')

    inv_mni_affine = pe.Node(interface=fsl.ConvertXFM(),
                             name='inv_mni_affine')
    inv_mni_affine.inputs.invert_xfm = True
    register_func_to_mni.connect(mni_affine, 'out_file',
                                 inv_mni_affine, 'in_file')
    register_func_to_mni.connect(inv_mni_affine, 'out_file',
                                 outputspec, 'mni_to_func_linear_xfm')

    register_func_to_mni.connect(inputspec, 'func',
                                 linear_reg, 'in_file')
    register_func_to_mni.connect(inputspec, 'anat',
                                 linear_reg, 'reference')
    register_func_to_mni.connect(inputspec, 'interp',
                                 linear_reg, 'interp')

    register_func_to_mni.connect(inputspec, 'func',
                                 mni_warp, 'in_file')
    register_func_to_mni.connect(inputspec, 'mni',
                                 mni_warp, 'ref_file')
    register_func_to_mni.connect(inputspec, 'anat_to_mni_nonlinear_xfm',
                                 mni_warp, 'field_file')

    register_func_to_mni.connect(linear_reg, 'out_matrix_file',
                                 mni_warp, 'premat')

    register_func_to_mni.connect(linear_reg, 'out_matrix_file',
                                 outputspec, 'func_to_anat_linear_xfm')
    register_func_to_mni.connect(mni_warp, 'out_file',
                                 outputspec, 'mni_func')

    return register_func_to_mni


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
                                                       'fieldmapmask'
                                                       ]),
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


def create_register_func_to_epi(
        name='register_func_to_epi', reg_option='ANTS', reg_ants_skull=1
):
    register_func_to_epi = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['func_4d',
                                                       'func_3d',
                                                       'func_3d_mask',
                                                       'epi',
                                                       'interp',
                                                       'ants_para']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['ants_initial_xfm',
                                                        'ants_rigid_xfm',
                                                        'ants_affine_xfm',
                                                        'warp_field',
                                                        'inverse_warp_field',
                                                        'fsl_flirt_xfm',
                                                        'fsl_fnirt_xfm',
                                                        'invlinear_xfm',
                                                        'func_in_epi']),
                         # 'func_mask_in_epi']),
                         name='outputspec')

    if reg_option == 'ANTS':
        # linear + non-linear registration
        func_to_epi_ants = \
            create_wf_calculate_ants_warp(
                name='func_to_epi_ants',
                num_threads=1,
                reg_ants_skull=1)

        register_func_to_epi.connect([
            (inputspec, func_to_epi_ants, [
                ('func_3d', 'inputspec.moving_brain'),
                ('epi', 'inputspec.reference_brain'),
                ('func_3d', 'inputspec.moving_skull'),
                ('epi', 'inputspec.reference_skull'),
                ('interp', 'inputspec.interp'),
                ('ants_para', 'inputspec.ants_para')
            ]),
        ])

        register_func_to_epi.connect([
            (func_to_epi_ants, outputspec, [
                ('outputspec.ants_initial_xfm', 'ants_initial_xfm'),
                ('outputspec.ants_rigid_xfm', 'ants_rigid_xfm'),
                ('outputspec.ants_affine_xfm', 'ants_affine_xfm'),
                ('outputspec.warp_field', 'warp_field'),
                ('outputspec.inverse_warp_field', 'inverse_warp_field'),
            ]),
        ])

        # combine transforms
        collect_transforms = pe.Node(util.Merge(4),
                                     name='collect_transforms_ants')
        register_func_to_epi.connect([
            (func_to_epi_ants, collect_transforms, [
                ('outputspec.ants_initial_xfm', 'in1'),
                ('outputspec.ants_rigid_xfm', 'in2'),
                ('outputspec.ants_affine_xfm', 'in3'),
                ('outputspec.warp_field', 'in4'),
            ]),
        ])

        # check transform list to exclude Nonetype (missing) init/rig/affine
        check_transform = pe.Node(util.Function(
            input_names=['transform_list'],
            output_names=['checked_transform_list', 'list_length'],
            function=check_transforms), name='{0}_check_transforms'.format(
            name))

        register_func_to_epi.connect(
            collect_transforms, 'out', check_transform, 'transform_list')

        # apply transform to func
        func_in_epi = pe.Node(
            interface=ants.ApplyTransforms(), name='func_in_epi_ants')
        func_in_epi.inputs.dimension = 3
        func_in_epi.inputs.input_image_type = 3
        register_func_to_epi.connect(
            inputspec, 'func_4d', func_in_epi, 'input_image')
        register_func_to_epi.connect(
            inputspec, 'epi', func_in_epi, 'reference_image')
        register_func_to_epi.connect(
            check_transform, 'checked_transform_list',
            func_in_epi, 'transforms')
        register_func_to_epi.connect(
            func_in_epi, 'output_image', outputspec, 'func_in_epi')

        # # apply transform to functional mask
        # func_mask_in_epi = pe.Node(
        #    interface=ants.ApplyTransforms(), name='func_mask_in_epi_ants')
        # func_mask_in_epi.inputs.dimension = 3
        # func_mask_in_epi.inputs.input_image_type = 0
        # register_func_to_epi.connect(
        #    inputspec, 'func_3d_mask', func_mask_in_epi, 'input_image')
        # register_func_to_epi.connect(
        #    inputspec, 'epi', func_mask_in_epi, 'reference_image')
        # register_func_to_epi.connect(
        #     check_transform, 'checked_transform_list',
        #    func_mask_in_epi, 'transforms')
        # register_func_to_epi.connect(
        #     func_mask_in_epi, 'output_image', outputspec, 'func_mask_in_epi')

    elif reg_option == 'FSL':
        # flirt linear registration
        func_to_epi_linear = pe.Node(
            interface=fsl.FLIRT(), name='func_to_epi_linear_fsl')
        func_to_epi_linear.inputs.dof = 6

        register_func_to_epi.connect(
            inputspec, 'func_3d', func_to_epi_linear, 'in_file')
        register_func_to_epi.connect(
            inputspec, 'epi', func_to_epi_linear, 'reference')
        register_func_to_epi.connect(
            func_to_epi_linear, 'out_matrix_file', outputspec,
            'fsl_flirt_xfm')

        inv_flirt_xfm = pe.Node(
            interface=fsl.utils.ConvertXFM(), name='inv_linear_reg0_xfm')
        inv_flirt_xfm.inputs.invert_xfm = True

        # fnirt non-linear registration
        func_to_epi_nonlinear = pe.Node(
            interface=fsl.FNIRT(), name='func_to_epi_nonlinear_fsl')
        func_to_epi_nonlinear.inputs.fieldcoeff_file = True

        register_func_to_epi.connect(
            inputspec, 'func_3d', func_to_epi_nonlinear, 'in_file')
        register_func_to_epi.connect(
            inputspec, 'epi', func_to_epi_nonlinear, 'ref_file')
        register_func_to_epi.connect(
            func_to_epi_linear, 'out_matrix_file',
            func_to_epi_nonlinear, 'affine_file')
        register_func_to_epi.connect(
            func_to_epi_nonlinear, 'fieldcoeff_file',
            outputspec, 'fsl_fnirt_xfm')

        register_func_to_epi.connect(
            func_to_epi_linear, 'out_matrix_file', inv_flirt_xfm, 'in_file')
        register_func_to_epi.connect(
            inv_flirt_xfm, 'out_file', outputspec, 'invlinear_xfm')

        # apply warp
        func_in_epi = pe.Node(
            interface=fsl.ApplyWarp(), name='func_in_epi_fsl')
        func_in_epi.inputs.interp = 'sinc'

        register_func_to_epi.connect(
            inputspec, 'func_4d', func_in_epi, 'in_file')
        register_func_to_epi.connect(inputspec, 'epi', func_in_epi,
                                     'ref_file')

        # --premat input disabled because it was throwing off the transform
        # application quality --- why, though?
        # register_func_to_epi.connect(
        # func_to_epi_linear, 'out_matrix_file', func_in_epi, 'premat')

        register_func_to_epi.connect(
            func_to_epi_nonlinear, 'fieldcoeff_file',
            func_in_epi, 'field_file')
        register_func_to_epi.connect(
            func_in_epi, 'out_file', outputspec, 'func_in_epi')

    return register_func_to_epi


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


def connect_func_to_anat_init_reg(workflow, strat_list, c):
    new_strat_list = []

    diff_complete = False
    if True in c.functional_registration['1-coregistration']['run']:

        for num_strat, strat in enumerate(strat_list):

            diff_complete = False
            if "despiked_fieldmap" in strat and "fieldmap_mask" in strat:
                diff_complete = True

            # if field map-based distortion correction is on, but BBR is off,
            # send in the distortion correction files here
            # TODO: is this robust to the possibility of forking both
            # TODO: distortion correction and BBR at the same time?
            # TODO: (note if you are forking with BBR on/off, at this point
            # TODO:  there is still only one strat, so you would have to fork
            # TODO:  here instead to have a func->anat with fieldmap and
            # TODO:  without, and send the without-fieldmap to the BBR fork)

            # TODO: if we're moving the distortion correction warp
            #       application, then the below is unnecessary
            '''
            dist_corr = False
            if (
                'diff_distcor' in nodes and
                True not in c.functional_registration[
                    '1-coregistration'
                ]['boundary_based_registration']['run']
            ):
                dist_corr = True
                # TODO: for now, disabling dist corr when BBR is disabled
                err = "\n\n[!] Field map distortion correction is enabled, " \
                    "but Boundary-Based Registration is off- BBR is " \
                    "required for distortion correction.\n\n"
                raise Exception(err)
            '''

            func_to_anat = create_register_func_to_anat(
                diff_complete, f'func_to_anat_FLIRT_{num_strat}')

            # Input registration parameters
            func_to_anat.inputs.inputspec.interp = 'trilinear'

            if 'Mean Functional' in c.functional_registration[
                '1-coregistration'
            ]['func_input_prep']['input']:
                # Input functional image (mean functional)
                node, out_file = strat['mean_functional']
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.func')

            elif 'Selected Functional Volume' in c.functional_registration[
                '1-coregistration'
            ]['func_input_prep']['input']:
                # Input functional image (specific volume)
                node, out_file = strat['selected_func_volume']
                workflow.connect(node, out_file,
                                 func_to_anat, 'inputspec.func')

            # Input skull-stripped anatomical
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

            if False in c.functional_registration['1-coregistration']['run']:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(func_to_anat.name)

            strat.update_resource_pool({
                'mean_functional_in_anat': (
                    func_to_anat, 'outputspec.anat_func_nobbreg'),
                'functional_to_anat_linear_xfm': (
                    func_to_anat,
                    'outputspec.func_to_anat_linear_xfm_nobbreg')
            })

    strat_list += new_strat_list

    return workflow, strat_list, diff_complete


def connect_func_to_anat_bbreg(workflow, strat_list, c, diff_complete):
    from CPAC.utils.utils import pick_wm

    new_strat_list = []

    if (
                    True in c.functional_registration['1-coregistration'][
                    'run'] and
                    True in c.functional_registration[
                    '1-coregistration'
                ]['boundary_based_registration']['run']
    ):

        template_based_segmentation = c.anatomical_preproc[
            'segmentation_workflow'
        ]['1-segmentation']['Template_Based']['template_for_segmentation']

        for num_strat, strat in enumerate(strat_list):

            # this is needed here in case tissue segmentation is set on/off
            # and you have bbreg enabled- this will ensure bbreg will run for
            # the strat that has segmentation but will not run (thus avoiding
            # a crash) on the strat without segmentation
            if 'anatomical_wm_mask' in strat:

                func_to_anat_bbreg = create_bbregister_func_to_anat(
                    diff_complete,
                    f'func_to_anat_bbreg_{num_strat}'
                )

                # Input registration parameters
                func_to_anat_bbreg.inputs.inputspec.bbr_schedule = \
                    c.functional_registration[
                        '1-coregistration'
                    ]['boundary_based_registration']['bbr_schedule']

                if 'Mean Functional' in c.functional_registration[
                    '1-coregistration'
                ]['func_input_prep']['input']:
                    # Input functional image (mean functional)
                    node, out_file = strat['mean_functional']
                    workflow.connect(node, out_file,
                                     func_to_anat_bbreg, 'inputspec.func')

                elif 'Selected Functional Volume' in \
                        c.functional_registration[
                            '1-coregistration'
                        ]['func_input_prep']['input']:
                    # Input functional image (specific volume)
                    node, out_file = strat['selected_func_volume']
                    workflow.connect(node, out_file,
                                     func_to_anat_bbreg, 'inputspec.func')

                # Input anatomical whole-head image (reoriented)
                node, out_file = strat['anatomical_skull_leaf']
                workflow.connect(node, out_file,
                                 func_to_anat_bbreg,
                                 'inputspec.anat_skull')

                node, out_file = strat['functional_to_anat_linear_xfm']
                workflow.connect(node, out_file,
                                 func_to_anat_bbreg,
                                 'inputspec.linear_reg_matrix')

                if (
                                True in c.anatomical_preproc[
                                'segmentation_workflow'
                            ]['1-segmentation']['Template_Based']['run'] or
                                True in c.anatomical_preproc[
                                'segmentation_workflow'
                            ]['1-segmentation']['ANTs_Prior_Based']['run']
                ):
                    # Input segmentation mask since template-based
                    # segmentation and ANTs-prior-based segmentation do not
                    # generate probability maps
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

                if False in c.functional_registration[
                    '1-coregistration'
                ]['boundary_based_registration']['run']:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(func_to_anat_bbreg.name)

                strat.update_resource_pool({
                    'mean_functional_in_anat': (
                        func_to_anat_bbreg, 'outputspec.anat_func'),
                    'functional_to_anat_linear_xfm': (
                        func_to_anat_bbreg,
                        'outputspec.func_to_anat_linear_xfm')
                }, override=True)

            else:
                # anatomical segmentation is not being run in this particular
                # strategy/fork - we don't want this to stop workflow building
                # unless there is only one strategy
                if len(strat_list) > 1:
                    pass
                else:
                    err = "\n\n[!] Boundary-based registration (BBR) " \
                          "for functional-to-anatomical registration is " \
                          "enabled, but anatomical segmentation is not. " \
                          "BBR requires the outputs of segmentation. " \
                          "Please modify your pipeline configuration and " \
                          "run again.\n\n"
                    raise Exception(err)

    strat_list += new_strat_list

    return workflow, strat_list


def apply_transform(wf_name, reg_tool):
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

    if reg_tool == 'ants':
        apply_warp = pe.Node(interface=ants.ApplyTransforms(),
                             name=f'apply_warp_{wf_name}')

        wf.connect(inputNode, 'input_image', apply_warp, 'input_image')
        wf.connect(inputNode, 'reference', apply_warp, 'reference_image')

        interp_string = pe.Node(util.Function(input_names=['interpolation',
                                                           'reg_tool'],
                                              output_names=['interpolation'],
                                              function=interpolation_string),
                                name=f'interp_string')
        interp_string.inputs.reg_tool = reg_tool

        wf.connect(inputNode, 'interpolation', interp_string, 'interpolation')
        wf.connect(interp_string, 'interpolation', apply_warp,
                   'interpolation')

        ants_xfm_list = pe.Node(util.Function(input_names=['transform'],
                                              output_names=['transform_list'],
                                              function=single_ants_xfm_to_list),
                                name=f'single_ants_xfm_to_list')

        wf.connect(inputNode, 'transform', ants_xfm_list, 'transform')
        wf.connect(ants_xfm_list, 'transform_list', apply_warp, 'transforms')

        wf.connect(apply_warp, 'output_image', outputNode, 'output_image')

    elif reg_tool == 'fsl':
        tissueprior_mni_to_t1 = pe.Node(interface=fsl.FLIRT(),
                                        name='{0}_prior_mni_to_t1'.format(
                                            wf_name))
        tissueprior_mni_to_t1.inputs.apply_xfm = True
        tissueprior_mni_to_t1.inputs.interp = 'nearestneighbour'

        # mni to t1
        preproc.connect(inputNode, 'tissue_prior',
                        tissueprior_mni_to_t1, 'in_file')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1,
                        'reference')

        preproc.connect(inputNode, 'standard2highres_mat',
                        tissueprior_mni_to_t1, 'in_matrix_file')

        preproc.connect(tissueprior_mni_to_t1, 'out_file',
                        overlap_segmentmap_with_prior, 'operand_files')

    return wf


def FSL_registration_connector(wf, cfg, strat_pool, pipe_num, opt=None,
                               symmetric=False):
    sym = ''
    symm = ''
    if symmetric:
        sym = 'sym'
        symm = '_symmetric'

    if opt == 'FSL' or opt == 'FSL-linear':
        # this is to prevent the user from running FNIRT if they are
        # providing already-skullstripped inputs. this is because
        # FNIRT requires an input with the skull still on

        flirt_reg_anat_mni = create_fsl_flirt_linear_reg(
            f'anat_mni_flirt_register{symm}_{pipe_num}'
        )

        # Input registration parameters
        flirt_reg_anat_mni.inputs.inputspec.interp = \
            cfg['registration_workflows']['anatomical_registration'][
                'FSL-FNIRT']['interpolation']

        node, out = strat_pool.get_data('desc-brain_T1w')
        wf.connect(node, out, flirt_reg_anat_mni, 'inputspec.input_brain')

        node, out = strat_pool.get_data(f'T1w_brain_template{symm}')
        wf.connect(node, out, flirt_reg_anat_mni,
                   'inputspec.reference_brain')

        outputs = {
            f'space-{sym}template_desc-brain_T1w': (
                flirt_reg_anat_mni, 'outputspec.output_brain'),
            f'from-T1w_to-{sym}template_mode-image_desc-linear_xfm': (
                flirt_reg_anat_mni, 'outputspec.linear_xfm'),
            f'from-{sym}template_to-T1w_mode-image_desc-linear_xfm': (
                flirt_reg_anat_mni, 'outputspec.invlinear_xfm'),
            f'from-T1w_to-{sym}template_mode-image_xfm': (
                flirt_reg_anat_mni, 'outputspec.linear_xfm')
        }

    if opt == 'FSL':
        fnirt_reg_anat_mni = create_fsl_fnirt_nonlinear_reg(
            f'anat_mni_fnirt_register{symm}_{pipe_num}'
        )

        node, out = strat_pool.get_data('desc-brain_T1w')
        wf.connect(node, out, fnirt_reg_anat_mni, 'inputspec.input_brain')

        node, out = strat_pool.get_data(f'T1w_brain_template{sym}')
        wf.connect(node, out, fnirt_reg_anat_mni, 'inputspec.reference_brain')

        node, out = strat_pool.get_data(["desc-preproc_T1w",
                                         "desc-reorient_T1w", "T1w"])
        wf.connect(node, out, fnirt_reg_anat_mni, 'inputspec.input_skull')

        # NOTE: crossover from above opt block
        wf.connect(flirt_reg_anat_mni, 'outputspec.linear_xfm',
                   fnirt_reg_anat_mni, 'inputspec.linear_aff')

        node, out = strat_pool.get_data(f'T1w_template{sym}')
        wf.connect(node, out, fnirt_reg_anat_mni, 'inputspec.reference_skull')

        ref_mask = 'ref_mask'
        if symmetric:
            ref_mask = 'dilated_symmetric_brain_mask'
        node, out = strat_pool.get_data(ref_mask)
        wf.connect(node, out, fnirt_reg_anat_mni, 'inputspec.ref_mask')

        # assign the FSL FNIRT config file specified in pipeline config.yml
        fnirt_reg_anat_mni.inputs.inputspec.fnirt_config = \
            cfg['registration_workflow'][
                'registration']['FSL-FNIRT']['fnirt_config']

        write_composite_xfm = pe.Node(interface=fsl.ConvertWarp(),
                                      name=f'combine_fsl_warps{symm}_'
                                           f'{pipe_num}')

        node, out = strat_pool.get_data(f'T1w_template{sym}')
        wf.connect(node, out, write_composite_xfm, 'reference')

        wf.connect(flirt_reg_anat_mni, 'outputspec.linear_xfm',
                   write_composite_xfm, 'premat')

        wf.connect(fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm',
                   write_composite_xfm, 'warp1')

        # NOTE: this is an UPDATE because of the opt block above
        added_outputs = {
            f'space-{sym}template_desc-brain_T1w': (
                fnirt_reg_anat_mni, 'outputspec.output_brain'),
            f'from-T1w_to-{sym}template_mode-image_desc-nonlinear_xfm': (
                fnirt_reg_anat_mni, 'outputspec.nonlinear_xfm'),
            f'from-T1w_to-{sym}template_mode-image_xfm': (
                write_composite_xfm, 'out_file')
        }
        outputs.update(added_outputs)

    return (wf, outputs)


def ANTs_registration_connector(wf, cfg, strat_pool, pipe_num, opt=None,
                                symmetric=False):
    sym = ''
    symm = ''
    if symmetric:
        sym = 'sym'
        symm = '_symmetric'

    if cfg.registration_workflows['anatomical_registration']['registration'][
        'ANTs']['T1_registration'] is None:
        err_msg = '\n\n[!] C-PAC says: \nYou have selected ANTs as your ' \
                  'anatomical registration method.\n' \
                  'However, no ANTs parameters were specified.\n' \
                  'Please specify ANTs parameters properly and try again.'
        raise Exception(err_msg)

    ants_reg_anat_mni = \
        create_wf_calculate_ants_warp(
            f'anat_mni_ants_register{symm}_{pipe_num}',
            num_threads=cfg.pipeline_setup['system_config'][
                'num_ants_threads'],
            reg_ants_skull=cfg['registration_workflows'][
                'anatomical_registration']['reg_with_skull']
        )
    ants_reg_anat_mni.inputs.inputspec.ants_para = \
        cfg.registration_workflows['anatomical_registration']['registration'][
            'ANTs']['T1_registration']
    ants_reg_anat_mni.inputs.inputspec.interp = \
        cfg.registration_workflows['anatomical_registration']['registration'][
            'ANTs']['interpolation']

    # calculating the transform with the skullstripped is
    # reported to be better, but it requires very high
    # quality skullstripping. If skullstripping is imprecise
    # registration with skull is preferred

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, ants_reg_anat_mni, 'inputspec.moving_brain')

    node, out = strat_pool.get_data(f'T1w_brain_template{symm}')
    wf.connect(node, out, ants_reg_anat_mni, 'inputspec.reference_brain')

    node, out = strat_pool.get_data(["desc-preproc_T1w", "desc-reorient_T1w",
                                     "T1w"])
    wf.connect(node, out, ants_reg_anat_mni, 'inputspec.moving_skull')

    node, out = strat_pool.get_data(f'T1w_template{symm}')
    wf.connect(node, out, ants_reg_anat_mni, 'inputspec.reference_skull')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, ants_reg_anat_mni, 'inputspec.moving_mask')

    if not symmetric:
        node, out = strat_pool.get_data('T1w_brain_template_mask')
    if symmetric:
        node, out = strat_pool.get_data('dilated_symmetric_brain_mask')
    wf.connect(node, out, ants_reg_anat_mni, 'inputspec.reference_mask')

    ants_reg_anat_mni.inputs.inputspec.fixed_image_mask = None

    if strat_pool.check_rpool('label-lesion_mask') and cfg[
        'registration_workflow']['registration']['ANTs'][
        'use_lesion_mask']:
        # Create lesion preproc node to apply afni Refit and Resample
        lesion_preproc = create_lesion_preproc(
            wf_name=f'lesion_preproc{symm}_{pipe_num}'
        )

        node, out = strat_pool.get_data('label-lesion_mask')

        wf.connect(node, out, lesion_preproc, 'inputspec.lesion')
        wf.connect(lesion_preproc, 'outputspec.reorient',
                   ants_reg_anat_mni, 'inputspec.fixed_image_mask')

    # combine the linear xfm's into one - makes it easier downstream
    write_composite_linear_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_linear{symm}_xfm_{pipe_num}',
        mem_gb=1.5)
    write_composite_linear_xfm.inputs.print_out_composite_warp_file = True
    write_composite_linear_xfm.inputs.output_image = \
        "from-T1w_to-template_mode-image_desc-linear_xfm.nii.gz"

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, write_composite_linear_xfm, 'input_image')

    node, out = strat_pool.get_data(f'T1w_brain_template{symm}')
    wf.connect(node, out, write_composite_linear_xfm, 'reference_image')

    write_composite_linear_xfm.inputs.input_image_type = 0
    write_composite_linear_xfm.inputs.dimension = 3
    write_composite_linear_xfm.inputs.interpolation = \
        cfg.registration_workflows['anatomical_registration']['registration'][
            'ANTs']['interpolation']

    collect_transforms = pe.Node(util.Merge(3),
                                 name=f'collect_transforms{symm}_{pipe_num}')

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
        name=f'check_transforms_{pipe_num}')

    wf.connect(collect_transforms, 'out', check_transform, 'transform_list')

    wf.connect(check_transform, 'checked_transform_list',
               write_composite_linear_xfm, 'transforms')

    # combine the linear xfm's into one - makes it easier downstream
    write_composite_invlinear_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_invlinear{symm}_xfm_{pipe_num}',
        mem_gb=1.5)
    write_composite_invlinear_xfm.inputs.print_out_composite_warp_file = True
    write_composite_invlinear_xfm.inputs.output_image = \
        "from-template_to-T1w_mode-image_desc-linear_xfm.nii.gz"

    node, out = strat_pool.get_data(f'T1w_brain_template{symm}')
    wf.connect(node, out, write_composite_invlinear_xfm, 'input_image')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, write_composite_invlinear_xfm, 'reference_image')

    write_composite_invlinear_xfm.inputs.input_image_type = 0
    write_composite_invlinear_xfm.inputs.dimension = 3
    write_composite_invlinear_xfm.inputs.interpolation = \
        cfg.registration_workflows['anatomical_registration']['registration'][
            'ANTs']['interpolation']

    collect_inv_transforms = pe.Node(util.Merge(3),
                                     name='collect_inv_transforms'
                                          f'{symm}_{pipe_num}')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_affine_xfm',
               collect_inv_transforms, 'in1')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_rigid_xfm',
               collect_inv_transforms, 'in2')

    wf.connect(ants_reg_anat_mni, 'outputspec.ants_initial_xfm',
               collect_inv_transforms, 'in3')

    # check transform list to exclude Nonetype (missing) init/rig/affine
    check_inv_transform = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['checked_transform_list',
                                    'list_length'],
                      function=check_transforms),
        name=f'check_inv_transforms_{pipe_num}')

    wf.connect(collect_transforms, 'out',
               check_inv_transform, 'transform_list')

    wf.connect(check_inv_transform, 'checked_transform_list',
               write_composite_invlinear_xfm, 'transforms')

    # generate inverse transform flags, which depends on the
    # number of transforms
    inverse_transform_flags = pe.Node(
        util.Function(input_names=['transform_list'],
                      output_names=['inverse_transform_flags'],
                      function=generate_inverse_transform_flags),
        name=f'inverse_transform_flags_{pipe_num}')

    wf.connect(check_inv_transform, 'checked_transform_list',
               inverse_transform_flags, 'transform_list')

    wf.connect(inverse_transform_flags, 'inverse_transform_flags',
               write_composite_invlinear_xfm, 'invert_transform_flags')

    # combine ALL xfm's into one - makes it easier downstream
    write_composite_xfm = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'write_composite_{symm}_xfm_{pipe_num}',
        mem_gb=1.5)
    write_composite_xfm.inputs.print_out_composite_warp_file = True
    write_composite_xfm.inputs.output_image = \
        "from-T1w_to-template_mode-image_xfm.nii.gz"

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, write_composite_xfm, 'input_image')

    node, out = strat_pool.get_data(f'T1w_brain_template{symm}')
    wf.connect(node, out, write_composite_xfm, 'reference_image')

    write_composite_xfm.inputs.input_image_type = 0
    write_composite_xfm.inputs.dimension = 3
    write_composite_xfm.inputs.interpolation = \
        cfg.registration_workflows['anatomical_registration']['registration'][
            'ANTs']['interpolation']

    collect_all_transforms = pe.Node(util.Merge(4),
                                     name=f'collect_all_transforms'
                                          f'{symm}_{pipe_num}')

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

    outputs = {
        f'space-{sym}template_desc-brain_T1w': (
            ants_reg_anat_mni, 'outputspec.normalized_output_brain'),
        f'from-T1w_to-{sym}template_mode-image_xfm': (
            write_composite_xfm, 'output_image'),
        f'from-T1w_to-{sym}template_mode-image_desc-linear_xfm': (
            write_composite_linear_xfm, 'output_image'),
        f'from-{sym}template_to-T1w_mode-image_desc-linear_xfm': (
            write_composite_invlinear_xfm, 'output_image'),
        f'from-T1w_to-{sym}template_mode-image_desc-nonlinear_xfm': (
            ants_reg_anat_mni, 'outputspec.warp_field'),
        f'from-{sym}template_to-T1w_mode-image_desc-nonlinear_xfm': (
            ants_reg_anat_mni, 'outputspec.inverse_warp_field')
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

    wf, outputs = FSL_registration_connector(wf, cfg, strat_pool, pipe_num,
                                             opt)

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

    wf, outputs = FSL_registration_connector(wf, cfg, strat_pool, pipe_num,
                                             opt, symmetric=True)

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
                 "from-T1w_to-template_mode-image_xfm"]}
    '''
    wf, outputs = ANTs_registration_connector(wf, cfg, strat_pool, pipe_num,
                                              opt)
    return (wf, outputs)


def register_symmetric_ANTs_anat_to_template(wf, cfg, strat_pool, pipe_num,
                                             opt=None):
    '''
    {"name": "register_symmetric_ANTs_anat_to_template",
     "config": ["registration_workflows", "anatomical_registration",
                "registration"],
     "switch": "run",
     "option_key": "using",
     "option_val": "ANTs",
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
                 "from-T1w_to-symtemplate_mode-image_xfm"]}
    '''

    wf, outputs = ANTs_registration_connector(wf, cfg, strat_pool, pipe_num,
                                              opt, symmetric=True)

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
            interface=ants.N4BiasFieldCorrection(dimension=3,
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
                 "from-bold_to-T1_mode-image_desc-linear_xfm"]}
    '''

    diff_complete = False
    if strat_pool.check_rpool("despiked_fieldmap") and \
            strat_pool.check_rpool("fieldmap_mask"):
        diff_complete = True

    # if field map-based distortion correction is on, but BBR is off,
    # send in the distortion correction files here
    # TODO: is this robust to the possibility of forking both
    # TODO: distortion correction and BBR at the same time?
    # TODO: (note if you are forking with BBR on/off, at this point
    # TODO:  there is still only one strat, so you would have to fork
    # TODO:  here instead to have a func->anat with fieldmap and
    # TODO:  without, and send the without-fieldmap to the BBR fork)

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
        'from-bold_to-T1_mode-image_desc-linear_xfm':
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
                "from-bold_to-T1_mode-image_desc-linear_xfm",
                "label-CSF_probseg",
                "diff_phase_dwell",
                "diff_phase_pedir",
                "despiked_fieldmap",
                "fieldmap_mask"],
     "outputs": ["space-T1w_desc-mean_bold",
                 "from-bold_to-T1_mode-image_desc-linear_xfm"]}
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
        'from-bold_to-T1_mode-image_desc-linear_xfm')
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
        'space-T1w_desc-mean_bold': (
        func_to_anat_bbreg, 'outputspec.anat_func'),
        'from-bold_to-T1_mode-image_desc-linear_xfm':
            (func_to_anat_bbreg, 'outputspec.func_to_anat_linear_xfm')
    }

    return (wf, outputs)
