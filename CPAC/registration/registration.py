import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl

def create_nonlinear_register(name='nonlinear_register'):
    """
    Performs non-linear registration of an input file to a reference file.

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
    
        inputspec.input_brain : string (nifti file)
            File of brain to be normalized (registered)
        inputspec.input_skull : string (nifti file)
            File of input brain with skull
        inputspec.reference_brain : string (nifti file)
            Target brain file to normalize to
        inputspec.reference_skull : string (nifti file)
            Target brain with skull to normalize to
        inputspec.fnirt_config : string (fsl fnirt config file)
            Configuration file containing parameters that can be specified in fnirt
            
    Workflow Outputs::
    
        outputspec.output_brain : string (nifti file)
            Normalizion of input brain file
        outputspec.linear_xfm : string (.mat file)
            Affine matrix of linear transformation of brain file
        outputspec.invlinear_xfm : string
            Inverse of affine matrix of linear transformation of brain file
        outputspec.nonlinear_xfm : string
            Nonlinear field coefficients file of nonlinear transformation
            
    Registration Procedure:
    
    1. Perform a linear registration to get affine transformation matrix.
    2. Perform a nonlinear registration on an input file to the reference file utilizing affine
       transformation from the previous step as a starting point.
    3. Invert the affine transformation to provide the user a transformation (affine only) from the
       space of the reference file to the input file.
       
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
                                                       'fnirt_config']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['output_brain',
                                                       'linear_xfm',
                                                       'invlinear_xfm',
                                                       'nonlinear_xfm']),
                         name='outputspec')
    
    linear_reg = pe.Node(interface =fsl.FLIRT(),
                         name='linear_reg_0')
    linear_reg.inputs.cost = 'corratio'
    
    nonlinear_reg = pe.Node(interface=fsl.FNIRT(),
                            name='nonlinear_reg_1')
    nonlinear_reg.inputs.fieldcoeff_file = True
    nonlinear_reg.inputs.jacobian_file = True
    nonlinear_reg.inputs.warp_resolution = (10,10,10)
    
    brain_warp = pe.Node(interface=fsl.ApplyWarp(),
                         name='brain_warp')
    
    
    inv_flirt_xfm = pe.Node(interface=fsl.utils.ConvertXFM(),
                            name='inv_linear_reg0_xfm')
    inv_flirt_xfm.inputs.invert_xfm = True

    nonlinear_register.connect(inputspec, 'input_brain',
                               linear_reg, 'in_file')
    nonlinear_register.connect(inputspec, 'reference_brain',
                               linear_reg, 'reference')
        
    nonlinear_register.connect(inputspec, 'input_skull',
                               nonlinear_reg, 'in_file')
    nonlinear_register.connect(inputspec, 'reference_skull',
                               nonlinear_reg, 'ref_file')
    nonlinear_register.connect(inputspec, 'fnirt_config',
                               nonlinear_reg, 'config_file')
    nonlinear_register.connect(linear_reg, 'out_matrix_file',
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

    nonlinear_register.connect(linear_reg, 'out_matrix_file',
                               inv_flirt_xfm, 'in_file')
    nonlinear_register.connect(inv_flirt_xfm, 'out_file',
                               outputspec, 'invlinear_xfm')

    nonlinear_register.connect(linear_reg, 'out_matrix_file',
                               outputspec, 'linear_xfm')
    
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
            
    Workflow Graph:
    
    .. image:: ../images/register_func_to_mni.dot.png
        :width: 500
        
    Detailed Workflow Graph:
    
    .. image:: ../images/register_func_to_mni_detailed.dot.png
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

def create_bbregister_func_to_mni(name='bbregister_func_to_mni'):
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
        inputspec.anat_wm_segmentation : string (nifti file)
            White matter segmentation probability mask in anatomical space
        inputspec.bbr_schedule : string (.sch file)
            Boundary based registration schedule file for flirt command
        
    Workflow Outputs::
    
        outputspec.func_to_anat_linear_xfm : string (mat file)
            Affine transformation from functional to anatomical native space
        outputspec.func_to_mni_linear_xfm : string (mat file)
            Affine transformation from functional to MNI space
        outputspec.mni_to_func_linear_xfm : string (mat file)
            Affine transformation from MNI to functional space
        outputspec.anat_wm_edge : string (nifti file)
            White matter edge mask in anatomical space
        outputspec.anat_func : string (nifti file)
            Functional data in anatomical space
        outputspec.mni_func : string (nifti file)
            Functional scan registered to MNI standard space
            
    Workflow Graph:
    
    .. image:: ../images/register_func_to_mni.dot.png
        :width: 500
        
    Detailed Workflow Graph:
    
    .. image:: ../images/register_func_to_mni_detailed.dot.png
        :width: 500
    """
    register_func_to_mni = pe.Workflow(name=name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['func',
                                                       'mni',
                                                       'anat',
                                                       'anat_skull',
                                                       'interp',
                                                       'anat_to_mni_nonlinear_xfm',
                                                       'anat_to_mni_linear_xfm',
                                                       'anat_wm_segmentation',
                                                       'bbr_schedule']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['func_to_anat_linear_xfm',
                                                        'func_to_mni_linear_xfm',
                                                        'mni_to_func_linear_xfm',
                                                        'anat_wm_edge',
                                                        'anat_func',
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

    wm_bb_mask = pe.Node(interface=fsl.ImageMaths(),
                         name='wm_bb_mask')
    wm_bb_mask.inputs.op_string = '-thr 0.5 -bin'
    register_func_to_mni.connect(inputspec, 'anat_wm_segmentation',
                                 wm_bb_mask, 'in_file')

    def wm_bb_edge_args(mas_file):
        return '-edge -bin -mas ' + mas_file

    wm_bb_edge = pe.Node(interface=fsl.ImageMaths(),
                         name='wm_bb_edge')
    register_func_to_mni.connect(wm_bb_mask, 'out_file',
                                 wm_bb_edge, 'in_file')
    register_func_to_mni.connect(wm_bb_mask, ('out_file', wm_bb_edge_args),
                                 wm_bb_edge, 'op_string')

    def bbreg_args(bbreg_target):
        return '-cost bbr -wmseg ' + bbreg_target

    bbreg_func_to_anat = pe.Node(interface=fsl.FLIRT(),
                                 name='bbreg_func_to_anat')
    bbreg_func_to_anat.inputs.dof = 6    
    register_func_to_mni.connect(inputspec, 'bbr_schedule',
                                 bbreg_func_to_anat, 'schedule')
    register_func_to_mni.connect(wm_bb_mask, ('out_file', bbreg_args),
                                 bbreg_func_to_anat, 'args')
    register_func_to_mni.connect(inputspec, 'func',
                                 bbreg_func_to_anat, 'in_file')
    register_func_to_mni.connect(inputspec, 'anat_skull',
                                 bbreg_func_to_anat, 'reference')
    register_func_to_mni.connect(linear_reg, 'out_matrix_file',
                                 bbreg_func_to_anat, 'in_matrix_file')
 
    register_func_to_mni.connect(inputspec, 'anat_to_mni_linear_xfm',
                                 mni_affine, 'in_file')    
    register_func_to_mni.connect(bbreg_func_to_anat, 'out_matrix_file',
                                 mni_affine, 'in_file2')
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
    
    register_func_to_mni.connect(bbreg_func_to_anat, 'out_matrix_file',
                                 mni_warp, 'premat')

    register_func_to_mni.connect(bbreg_func_to_anat, 'out_matrix_file',
                                 outputspec, 'func_to_anat_linear_xfm')
    register_func_to_mni.connect(bbreg_func_to_anat, 'out_file',
                                 outputspec, 'anat_func')
    register_func_to_mni.connect(mni_warp, 'out_file',
                                 outputspec, 'mni_func')
    register_func_to_mni.connect(wm_bb_edge, 'out_file',
                                 outputspec, 'anat_wm_edge')
    
    return register_func_to_mni