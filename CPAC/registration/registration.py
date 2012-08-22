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
    
        inputspec.input : string (nifti file)
            File to be normalized (registered)
        inputspec.reference : string (nifti file)
            Target file to normalize to
        inputspec.fnirt_config : string (fsl fnirt config file)
            Configuration file containing parameters that can be specified in fnirt
            
    Workflow Outputs::
    
        outputspec.output : string (nifti file)
            Normalizion of input file
        outputspec.linear_xfm : string (.mat file)
            Affine matrix of linear transformation
        outputspec.linear_invxfm : string
            Inverse of affine matrix of linear transformation
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
    
    inputspec = pe.Node(util.IdentityInterface(fields=['input',
                                                       'reference',
                                                       'fnirt_config']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['output',
                                                       'linear_xfm',
                                                       'linear_invxfm',
                                                       'nonlinear_xfm']),
                         name='outputspec')
    
    linear_reg = pe.Node(interface =fsl.FLIRT(),
                         name='linear_reg_0')
    linear_reg.inputs.cost = 'corratio'
    linear_reg.inputs.dof = 6
    linear_reg.inputs.interp = 'nearestneighbour'
    
    nonlinear_reg = pe.Node(interface=fsl.FNIRT(),
                            name='nonlinear_reg_1')
    nonlinear_reg.inputs.fieldcoeff_file = True
    nonlinear_reg.inputs.jacobian_file = True
    nonlinear_reg.inputs.warp_resolution = (10,10,10)
    
    inv_flirt_xfm = pe.Node(interface=fsl.utils.ConvertXFM(),
                            name='inv_linear_reg0_xfm')
    inv_flirt_xfm.inputs.invert_xfm = True

    nonlinear_register.connect(inputspec, 'input',
                               linear_reg, 'in_file')
    nonlinear_register.connect(inputspec, 'reference',
                               linear_reg, 'reference')
        
    nonlinear_register.connect(inputspec, 'input',
                               nonlinear_reg, 'in_file')
    nonlinear_register.connect(inputspec, 'reference',
                               nonlinear_reg, 'ref_file')
    nonlinear_register.connect(inputspec, 'fnirt_config',
                               nonlinear_reg, 'config_file')
    nonlinear_register.connect(linear_reg, 'out_matrix_file',
                               nonlinear_reg, 'affine_file')
    nonlinear_register.connect(nonlinear_reg, 'warped_file',
                               outputspec, 'output')
    nonlinear_register.connect(nonlinear_reg, 'fieldcoeff_file',
                               outputspec, 'nonlinear_xfm')

    nonlinear_register.connect(linear_reg, 'out_matrix_file',
                               inv_flirt_xfm, 'in_file')
    nonlinear_register.connect(inv_flirt_xfm, 'out_file',
                               outputspec, 'linear_invxfm')

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
        inputspec.anat_to_mni_xfm : string (warp file)
            Corresponding anatomical native space to MNI warp file
            
    Workflow Outputs::
    
        outputspec.func_to_anat_xfm : string (mat file)
            Affine transformation from functional to anatomical native space
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
                                                       'anat_to_mni_xfm']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['func_to_anat_xfm',
                                                        'mni_func']),
                         name='outputspec')
    
    linear_reg = pe.Node(interface=fsl.FLIRT(),
                         name='linear_func_to_anat')
    linear_reg.inputs.cost = 'corratio'
    linear_reg.inputs.dof = 6
    linear_reg.inputs.interp = 'nearestneighbour'
    
    mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                       name='mni_warp')
    
    register_func_to_mni.connect(inputspec, 'func',
                                 linear_reg, 'in_file')
    register_func_to_mni.connect(inputspec, 'anat',
                                 linear_reg, 'reference')
    
    register_func_to_mni.connect(inputspec, 'anat_to_mni_xfm',
                                 mni_warp, 'field_file')
    
    register_func_to_mni.connect(linear_reg, 'out_matrix_file',
                                 mni_warp, 'premat')
    
    return register_func_to_mni