import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.c3 as c3
import nipype.interfaces.ants as ants

from CPAC.utils.interfaces.function import Function
from CPAC.registration.utils import seperate_warps_list, \
                                    check_transforms, \
                                    generate_inverse_transform_flags, \
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
            Configuration file containing parameters that can be specified in fnirt
    Workflow Outputs::
    
        outputspec.output_brain : string (nifti file)
            Normalizion of input brain file
        outputspec.nonlinear_xfm : string
            Nonlinear field coefficients file of nonlinear transformation
            
    Registration Procedure:

    1. Perform a nonlinear registration on an input file to the reference file utilizing affine
       transformation from the previous step as a starting point.
    2. Invert the affine transformation to provide the user a transformation (affine only) from the
       space of the reference file to the input file.
       
    .. exec::
        from CPAC.registration import create_fsl_fnirt_nonlinear_reg
        wf = create_fsl_fnirt_nonlinear_reg()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/nonlinear_register.dot'
        )

    Workflow Graph:
    
    .. image:: ../../images/generated/nonlinear_register.png
        :width: 500
    
    Detailed Workflow Graph:
    
    .. image:: ../../images/generated/nonlinear_register_detailed.png
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
    Registers a functional scan in native space to MNI standard space.  This is meant to be used 
    after create_nonlinear_register() has been run and relies on some of it's outputs.

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
            Type of interpolation to use ('trilinear' or 'nearestneighbour' or 'sinc')
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
    
    .. image:: ../../images/generated/register_func_to_mni_detailed.png
        :width: 500
    """
    register_func_to_mni = pe.Workflow(name=name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['func',
                                                       'mni',
                                                       'anat',
                                                       'interp',
                                                       'anat_to_mni_nonlinear_xfm',
                                                       'anat_to_mni_linear_xfm']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['func_to_anat_linear_xfm',
                                                        'func_to_mni_linear_xfm',
                                                        'mni_to_func_linear_xfm',
                                                        'mni_func']),
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
            Type of interpolation to use ('trilinear' or 'nearestneighbour' or 'sinc')
            
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

    outputspec = pe.Node(util.IdentityInterface(fields=['func_to_anat_linear_xfm_nobbreg',
                                                        'anat_func_nobbreg']),
                         name='outputspec')
    
    linear_reg = pe.Node(interface=fsl.FLIRT(),
                         name='linear_func_to_anat')
    linear_reg.inputs.cost = 'corratio'
    linear_reg.inputs.dof = 6
    
    #if fieldmap_distortion:

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
        register_func_to_anat.connect(inputNode_pedir, ('pedir', convert_pedir),
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
    Registers a functional scan in native space to structural.  This is meant to be used 
    after create_nonlinear_register() has been run and relies on some of it's outputs.

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

    outputspec = pe.Node(util.IdentityInterface(fields=['func_to_anat_linear_xfm',
                                                        'anat_func']),
                         name='outputspec')

    wm_bb_mask = pe.Node(interface=fsl.ImageMaths(),
                         name='wm_bb_mask')
    wm_bb_mask.inputs.op_string = '-thr 0.5 -bin'

    register_bbregister_func_to_anat.connect(inputspec, 'anat_wm_segmentation',
                                             wm_bb_mask, 'in_file')

    def bbreg_args(bbreg_target):
        return '-cost bbr -wmseg ' + bbreg_target

    bbreg_func_to_anat = pe.Node(interface=fsl.FLIRT(),
                                 name='bbreg_func_to_anat')
    bbreg_func_to_anat.inputs.dof = 6    
 
    register_bbregister_func_to_anat.connect(inputspec, 'bbr_schedule',
                                 bbreg_func_to_anat, 'schedule')
 
    register_bbregister_func_to_anat.connect(wm_bb_mask, ('out_file', bbreg_args),
                                 bbreg_func_to_anat, 'args')
 
    register_bbregister_func_to_anat.connect(inputspec, 'func',
                                 bbreg_func_to_anat, 'in_file')
 
    register_bbregister_func_to_anat.connect(inputspec, 'anat_skull',
                                 bbreg_func_to_anat, 'reference')
 
    register_bbregister_func_to_anat.connect(inputspec, 'linear_reg_matrix',
                                 bbreg_func_to_anat, 'in_matrix_file')

    #if fieldmap_distortion:

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
        register_bbregister_func_to_anat.connect(inputNode_pedir, ('pedir', convert_pedir),
                                                 bbreg_func_to_anat, 'pedir')
        register_bbregister_func_to_anat.connect(inputspec, 'fieldmap',
                                                 bbreg_func_to_anat, 'fieldmap')
        register_bbregister_func_to_anat.connect(inputspec, 'fieldmapmask',
                                                 bbreg_func_to_anat, 'fieldmapmask')
        register_bbregister_func_to_anat.connect(inputNode_echospacing, 'echospacing',
                                                 bbreg_func_to_anat, 'echospacing')

    register_bbregister_func_to_anat.connect(bbreg_func_to_anat, 'out_matrix_file',
                                 outputspec, 'func_to_anat_linear_xfm')
    
    register_bbregister_func_to_anat.connect(bbreg_func_to_anat, 'out_file',
                                 outputspec, 'anat_func')
    
    return register_bbregister_func_to_anat
    

def create_register_func_to_epi(name='register_func_to_epi', reg_option='ANTS', reg_ants_skull=1):

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
        collect_transforms = pe.Node(util.Merge(4), name='collect_transforms_ants')
        register_func_to_epi.connect([
            (func_to_epi_ants, collect_transforms, [
                ('outputspec.ants_initial_xfm', 'in1'),
                ('outputspec.ants_rigid_xfm', 'in2'),
                ('outputspec.ants_affine_xfm', 'in3'),
                ('outputspec.warp_field', 'in4'),
            ]),
        ])

        # check transform list to exclude Nonetype (missing) init/rig/affine
        check_transform = pe.Node(util.Function(input_names=['transform_list'], 
                                                output_names=['checked_transform_list', 'list_length'],
                                                function=check_transforms), name='{0}_check_transforms'.format(name))
        
        register_func_to_epi.connect(collect_transforms, 'out', check_transform, 'transform_list')


        # apply transform to func 
        func_in_epi = pe.Node(interface=ants.ApplyTransforms(), name='func_in_epi_ants')
        func_in_epi.inputs.dimension = 3
        func_in_epi.inputs.input_image_type = 3
        register_func_to_epi.connect(inputspec, 'func_4d', func_in_epi, 'input_image')
        register_func_to_epi.connect(inputspec, 'epi', func_in_epi, 'reference_image')
        register_func_to_epi.connect(check_transform, 'checked_transform_list', func_in_epi, 'transforms')
        register_func_to_epi.connect(func_in_epi, 'output_image', outputspec, 'func_in_epi')

        # # apply transform to functional mask
        # func_mask_in_epi = pe.Node(interface=ants.ApplyTransforms(), name='func_mask_in_epi_ants')
        # func_mask_in_epi.inputs.dimension = 3
        # func_mask_in_epi.inputs.input_image_type = 0
        # register_func_to_epi.connect(inputspec, 'func_3d_mask', func_mask_in_epi, 'input_image')
        # register_func_to_epi.connect(inputspec, 'epi', func_mask_in_epi, 'reference_image')
        # register_func_to_epi.connect(check_transform, 'checked_transform_list', func_mask_in_epi, 'transforms')
        # register_func_to_epi.connect(func_mask_in_epi, 'output_image', outputspec, 'func_mask_in_epi')

    elif reg_option == 'FSL':
        # flirt linear registration 
        func_to_epi_linear = pe.Node(interface=fsl.FLIRT(), name='func_to_epi_linear_fsl')
        func_to_epi_linear.inputs.dof = 6

        register_func_to_epi.connect(inputspec, 'func_3d', func_to_epi_linear, 'in_file')
        register_func_to_epi.connect(inputspec, 'epi', func_to_epi_linear, 'reference')
        register_func_to_epi.connect(func_to_epi_linear, 'out_matrix_file', outputspec, 'fsl_flirt_xfm')
        
        inv_flirt_xfm = pe.Node(interface=fsl.utils.ConvertXFM(), name='inv_linear_reg0_xfm')
        inv_flirt_xfm.inputs.invert_xfm = True

        # fnirt non-linear registration
        func_to_epi_nonlinear = pe.Node(interface=fsl.FNIRT(), name='func_to_epi_nonlinear_fsl')
        func_to_epi_nonlinear.inputs.fieldcoeff_file = True

        register_func_to_epi.connect(inputspec, 'func_3d', func_to_epi_nonlinear, 'in_file')
        register_func_to_epi.connect(inputspec, 'epi', func_to_epi_nonlinear, 'ref_file')
        register_func_to_epi.connect(func_to_epi_linear, 'out_matrix_file', func_to_epi_nonlinear, 'affine_file')
        register_func_to_epi.connect(func_to_epi_nonlinear, 'fieldcoeff_file', outputspec, 'fsl_fnirt_xfm')

        register_func_to_epi.connect(func_to_epi_linear, 'out_matrix_file', inv_flirt_xfm, 'in_file')
        register_func_to_epi.connect(inv_flirt_xfm, 'out_file', outputspec, 'invlinear_xfm')

        # apply warp
        func_in_epi = pe.Node(interface=fsl.ApplyWarp(), name='func_in_epi_fsl')
        func_in_epi.inputs.interp = 'sinc'

        register_func_to_epi.connect(inputspec, 'func_4d', func_in_epi, 'in_file')
        register_func_to_epi.connect(inputspec, 'epi', func_in_epi, 'ref_file')

        # --premat input disabled because it was throwing off the transform
        # application quality --- why, though?
        #register_func_to_epi.connect(func_to_epi_linear, 'out_matrix_file', func_in_epi, 'premat')

        register_func_to_epi.connect(func_to_epi_nonlinear, 'fieldcoeff_file', func_in_epi, 'field_file')
        register_func_to_epi.connect(func_in_epi, 'out_file', outputspec, 'func_in_epi')

    return register_func_to_epi


def create_wf_calculate_ants_warp(name='create_wf_calculate_ants_warp', num_threads=1, reg_ants_skull=1):

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
            Type of interpolation to use ('Linear' or 'BSpline' or 'LanczosWindowedSinc')

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
    
    .. image:: ../../images/generated/calculate_ants_warp.png
        :width: 500

    Detailed Workflow Graph:
    
    .. image:: ../../images/generated/calculate_ants_warp_detailed.png
        :width: 500      
    '''

    calc_ants_warp_wf = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(
        fields=['moving_brain',
                'reference_brain',
                'moving_skull',
                'reference_skull',
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
                                                     'fixed_image_mask',
                                                     'interp'],
                                        output_names=['warp_list',
                                                      'warped_image'],
                                        function=hardcoded_reg,
                                        imports=reg_imports),
                name='calc_ants_warp')

    calculate_ants_warp.interface.num_threads = num_threads

    select_forward_initial = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=seperate_warps_list), name='select_forward_initial')

    select_forward_initial.inputs.selection = "Initial"

    select_forward_rigid = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=seperate_warps_list), name='select_forward_rigid')

    select_forward_rigid.inputs.selection = "Rigid"

    select_forward_affine = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=seperate_warps_list), name='select_forward_affine')

    select_forward_affine.inputs.selection = "Affine"

    select_forward_warp = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=seperate_warps_list), name='select_forward_warp')

    select_forward_warp.inputs.selection = "Warp"

    select_inverse_warp = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=seperate_warps_list), name='select_inverse_warp')

    select_inverse_warp.inputs.selection = "Inverse"

    calc_ants_warp_wf.connect(inputspec, 'moving_brain',
            calculate_ants_warp, 'moving_brain')

    calc_ants_warp_wf.connect(inputspec, 'reference_brain',
            calculate_ants_warp, 'reference_brain')

    if reg_ants_skull == 1 and not reg_ants_skull == 0:
        calc_ants_warp_wf.connect(inputspec, 'moving_skull',
                calculate_ants_warp, 'moving_skull')

        calc_ants_warp_wf.connect(inputspec, 'reference_skull',
                calculate_ants_warp, 'reference_skull')

    else:
        calc_ants_warp_wf.connect(inputspec, 'moving_brain',
                calculate_ants_warp, 'moving_skull')

        calc_ants_warp_wf.connect(inputspec, 'reference_brain',
                calculate_ants_warp, 'reference_skull')

    calc_ants_warp_wf.connect(inputspec, 'fixed_image_mask',
            calculate_ants_warp, 'fixed_image_mask')

    calc_ants_warp_wf.connect(inputspec, 'ants_para',
            calculate_ants_warp, 'ants_para')

    calc_ants_warp_wf.connect(inputspec, 'interp',
            calculate_ants_warp, 'interp')

    # inter-workflow connections

    calc_ants_warp_wf.connect(calculate_ants_warp, 'warp_list',
            select_forward_initial, 'warp_list')

    calc_ants_warp_wf.connect(calculate_ants_warp, 'warp_list',
            select_forward_rigid, 'warp_list')

    calc_ants_warp_wf.connect(calculate_ants_warp, 'warp_list',
            select_forward_affine, 'warp_list')

    calc_ants_warp_wf.connect(calculate_ants_warp, 'warp_list',
            select_forward_warp, 'warp_list')

    calc_ants_warp_wf.connect(calculate_ants_warp, 'warp_list',
            select_inverse_warp, 'warp_list')

    # connections to outputspec

    calc_ants_warp_wf.connect(select_forward_initial, 'selected_warp',
            outputspec, 'ants_initial_xfm')

    calc_ants_warp_wf.connect(select_forward_rigid, 'selected_warp',
            outputspec, 'ants_rigid_xfm')

    calc_ants_warp_wf.connect(select_forward_affine, 'selected_warp',
            outputspec, 'ants_affine_xfm')

    calc_ants_warp_wf.connect(select_forward_warp, 'selected_warp',
            outputspec, 'warp_field')

    calc_ants_warp_wf.connect(select_inverse_warp, 'selected_warp',
            outputspec, 'inverse_warp_field')

    calc_ants_warp_wf.connect(calculate_ants_warp, 'warped_image',
            outputspec, 'normalized_output_brain')

    return calc_ants_warp_wf


def connect_func_to_anat_init_reg(workflow, strat_list, c):

    new_strat_list = []

    diff_complete = False
    
    if 1 in c.runRegisterFuncToAnat:

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
            if 'diff_distcor' in nodes and 1 not in c.runBBReg:
                dist_corr = True
                # TODO: for now, disabling dist corr when BBR is disabled
                err = "\n\n[!] Field map distortion correction is enabled, " \
                    "but Boundary-Based Registration is off- BBR is " \
                    "required for distortion correction.\n\n"
                raise Exception(err)
            '''

            func_to_anat = create_register_func_to_anat(diff_complete,
                                                        f'func_to_anat_FLIRT_{num_strat}')

            # Input registration parameters
            func_to_anat.inputs.inputspec.interp = 'trilinear'

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

            if 0 in c.runRegisterFuncToAnat:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(func_to_anat.name)

            strat.update_resource_pool({
                'mean_functional_in_anat': (func_to_anat, 'outputspec.anat_func_nobbreg'),
                'functional_to_anat_linear_xfm': (func_to_anat, 'outputspec.func_to_anat_linear_xfm_nobbreg')
            })

    strat_list += new_strat_list

    return workflow, strat_list, diff_complete


def connect_func_to_anat_bbreg(workflow, strat_list, c, diff_complete):

    from CPAC.utils.utils import pick_wm

    new_strat_list = []

    if 1 in c.runRegisterFuncToAnat and 1 in c.runBBReg:

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
                    c.boundaryBasedRegistrationSchedule

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
                node, out_file = strat['anatomical_skull_leaf']
                workflow.connect(node, out_file,
                                    func_to_anat_bbreg,
                                    'inputspec.anat_skull')

                node, out_file = strat['functional_to_anat_linear_xfm']
                workflow.connect(node, out_file,
                                    func_to_anat_bbreg,
                                    'inputspec.linear_reg_matrix')

                if 'T1_template' in c.template_based_segmentation or \
                        'EPI_template' in c.template_based_segmentation or \
                            1 in c.ANTs_prior_based_segmentation :
                    # Input segmentation mask,
                    # since template_based_segmentation or ANTs_prior_based_segmentation cannot generate
                    # probability maps
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

                if 0 in c.runBBReg:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(func_to_anat_bbreg.name)

                strat.update_resource_pool({
                    'mean_functional_in_anat': (func_to_anat_bbreg, 'outputspec.anat_func'),
                    'functional_to_anat_linear_xfm': (func_to_anat_bbreg, 'outputspec.func_to_anat_linear_xfm')
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


def connect_func_to_template_reg(workflow, strat_list, c):
    
    from CPAC.registration import output_func_to_standard

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
                if reg.lower() == 'ants' and c.ANTs_para_EPI_registration is None:
                    err_msg = '\n\n[!] C-PAC says: \n'\
                        "You have selected \'regOption: [{0}]\' and \'runRegisterFuncToTemplate :  ['{1}']\'. \n"\
                                'However, no EPI-to-template ANTs parameters were specified. ' \
                                'Please specify ANTs parameters properly and try again'.format(str(c.regOption),
                                                                                               str(c.runRegisterFuncToTemplate))
                    raise Exception(err_msg)
                elif reg.lower() == 'ants':
                    func_to_epi.inputs.inputspec.ants_para = c.ANTs_para_EPI_registration
                    func_to_epi.inputs.inputspec.interp = c.funcRegANTSinterpolation
                else:
                    func_to_epi.inputs.inputspec.interp = c.funcRegFSLinterpolation

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

                if reg.lower() == 'fsl':
                    strat.update_resource_pool({
                        'epi_registration_method': 'FSL',
                        'func_to_epi_linear_xfm': (func_to_epi, 'outputspec.fsl_flirt_xfm'),  
                        'func_to_epi_nonlinear_xfm': (func_to_epi, 'outputspec.fsl_fnirt_xfm'),
                        'epi_to_func_linear_xfm': (func_to_epi, 'outputspec.invlinear_xfm'),
                    })

                elif reg.lower() == 'ants':
                    strat.update_resource_pool({
                        'epi_registration_method': 'ANTS',
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

        if 'T1_template' in c.runRegisterFuncToTemplate and \
                'functional_to_epi-standard' not in strat:

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

    return workflow, strat_list
