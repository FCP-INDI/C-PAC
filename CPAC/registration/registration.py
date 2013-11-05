import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants
from nipype.interfaces.ants import ApplyTransforms
import nipype.interfaces.c3 as c3

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
    
    # FNIRT parameters are specified by FSL config file
    # ${FSLDIR}/etc/flirtsch/TI_2_MNI152_2mm.cnf (or user-specified)
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



def create_register_func_to_anat(name='register_func_to_anat'):

    
    register_func_to_anat = pe.Workflow(name=name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['func',
                                                       'anat',
                                                       'interp']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['func_to_anat_linear_xfm_nobbreg',
                                                        #'func_to_mni_linear_xfm',
                                                        #'mni_to_func_linear_xfm',
                                                        #'anat_wm_edge',
                                                        'anat_func_nobbreg']),
                         name='outputspec')
    
    linear_reg = pe.Node(interface=fsl.FLIRT(),
                         name='linear_func_to_anat')
    linear_reg.inputs.cost = 'corratio'
    linear_reg.inputs.dof = 6

  

    register_func_to_anat.connect(inputspec, 'func',
                                 linear_reg, 'in_file')
    
    register_func_to_anat.connect(inputspec, 'anat',
                                 linear_reg, 'reference')
    
    register_func_to_anat.connect(inputspec, 'interp',
                                 linear_reg, 'interp')

    register_func_to_anat.connect(linear_reg, 'out_matrix_file',
                                 outputspec, 'func_to_anat_linear_xfm_nobbreg')

    register_func_to_anat.connect(linear_reg, 'out_file',
                                 outputspec, 'anat_func_nobbreg')

    
    return register_func_to_anat



def create_bbregister_func_to_anat(name='bbregister_func_to_anat'):
  
    """
    Registers a functional scan in native space to structural.  This is meant to be used 
    after create_nonlinear_register() has been run and relies on some of it's outputs.

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    register_func_to_anat : nipype.pipeline.engine.Workflow

    Notes
    -----
    
    !!!! SOME OF THE BELOW HAS CHANGED, OBSOLETE CODE HAS BEEN REMOVED

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
    
    register_bbregister_func_to_anat = pe.Workflow(name=name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['func',
                                                       'anat_skull',
                                                       'linear_reg_matrix',
                                                       'anat_wm_segmentation',
                                                       'bbr_schedule']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['func_to_anat_linear_xfm',
                                                        #'func_to_mni_linear_xfm',
                                                        #'mni_to_func_linear_xfm',
                                                        #'anat_wm_edge',
                                                        'anat_func']),
                         name='outputspec')
    


    wm_bb_mask = pe.Node(interface=fsl.ImageMaths(),
                         name='wm_bb_mask')
    wm_bb_mask.inputs.op_string = '-thr 0.5 -bin'


    register_bbregister_func_to_anat.connect(inputspec, 'anat_wm_segmentation',
                                 wm_bb_mask, 'in_file')

    def wm_bb_edge_args(mas_file):
        return '-edge -bin -mas ' + mas_file


    #wm_bb_edge = pe.Node(interface=fsl.ImageMaths(),
    #                     name='wm_bb_edge')


    #register_func_to_mni.connect(wm_bb_mask, 'out_file',
    #                             wm_bb_edge, 'in_file')

    #register_func_to_mni.connect(wm_bb_mask, ('out_file', wm_bb_edge_args),
    #                             wm_bb_edge, 'op_string')

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

    register_bbregister_func_to_anat.connect(bbreg_func_to_anat, 'out_matrix_file',
                                 outputspec, 'func_to_anat_linear_xfm')
    
    register_bbregister_func_to_anat.connect(bbreg_func_to_anat, 'out_file',
                                 outputspec, 'anat_func')
   
    
    #register_func_to_mni.connect(wm_bb_edge, 'out_file',
    #                             outputspec, 'anat_wm_edge')
    
    return register_bbregister_func_to_anat
    
    

def create_ants_nonlinear_xfm(name='ants_nonlinear_xfm'):

    """
    Calculates the nonlinear ANTS registration transform.

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    ants_nonlinear_xfm : nipype.pipeline.engine.Workflow

    Notes
    -----
    
    Workflow Inputs::
    
        inputspec.anatomical_brain : string (nifti file)
            File of brain to be normalized (registered)
        inputspec.reference_brain : string (nifti file)
            Target brain file to normalize to

            
    Workflow Outputs::
    
        outputspec.warp_field : string (nifti file)
            Output warp field of registration
        outputspec.affine_transformation : text file
            Affine matrix of nonlinear transformation of brain file
        outputspec.inverse_warp : string (nifti file)
            Inverse of the warp field of the registration
        outputspec.output_brain : string (nifti file)
            Template-registered version of input brain
            
    Registration Procedure:
    
    1. Performs a nonlinear anatomical-to-template registration.
       
    Workflow Graph:
    
    .. image:: 
        :width: 500
    
    Detailed Workflow Graph:
    
    .. image:: 
        :width: 500      
    """

    import nipype.interfaces.ants as ants
    from nipype.interfaces.ants.legacy import GenWarpFields

    ants_nonlinear_xfm = pe.Workflow(name=name)

    inputspec = pe.Node(
        util.IdentityInterface(
            fields=['anatomical_brain', 'reference_brain']),
        name='inputspec')


    # use ANTS to warp the masked anatomical image to a template image
    warp_brain = pe.Node(GenWarpFields(),
                         name='warp_brain')

    outputspec = pe.Node(
        util.IdentityInterface(
            fields=['warp_field', 'affine_transformation',
                'inverse_warp', 'output_brain']),
        name='outputspec')


    ants_nonlinear_xfm.connect(inputspec, 'anatomical_brain', warp_brain, 'input_image')
    ants_nonlinear_xfm.connect(inputspec, 'reference_brain', warp_brain, 'reference_image')

    ants_nonlinear_xfm.connect(warp_brain, 'warp_field', outputspec, 'warp_field')
    ants_nonlinear_xfm.connect(warp_brain, 'affine_transformation', outputspec, 'affine_transformation')
    ants_nonlinear_xfm.connect(warp_brain, 'inverse_warp_field', outputspec, 'inverse_warp')
    ants_nonlinear_xfm.connect(warp_brain, 'output_file', outputspec, 'output_brain')


    return ants_nonlinear_xfm



def create_apply_ants_xfm(dimension, mapnode, name='apply_ants_xfm'):

    """
    Takes in the results of the FSL-based functional-to-anatomical registration,
    and the results of the ANTS-based anatomical-to-template registration, and
    applies these transformations to register functional to template.

    The FSL-based functional-to-anatomical registration output transformations
    are first converted from FSL format to ITK format - this step can and
    should be separated into their own workflows in the future.

    NOTE: The dimension of the input image (3 or 4) must be specified in the
          function call, as well as whether or not the apply warp must be
          treated as a mapnode (0 - no, 1 - yes).

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    apply_ants_xfm : nipype.pipeline.engine.Workflow

    Notes
    -----
    
    Workflow Inputs::

        inputspec.in_file : string (nifti file)
            File of functional brain data to be registered
        inputspec.warp_reference : string (nifti file)
            File of template to be used
        inputspec.use_nearest : boolean (True or False)
            Whether or not to use nearest neighbor interpolation
        inputspec.func_anat_affine : .mat file (affine matrix)
            Output matrix of FSL-based functional to anatomical registration
        inputspec.conversion_reference : string (nifti file)
            File of skull-stripped anatomical brain to be used in affine conversion
        inputspec.conversion_source : string (nifti file)
            Should match the input of the apply warp (in_file) unless you are
            applying the warp to a 4-d file, in which case this file should
            be a mean_functional file
        inputspec.nonlinear_field : string (nifti file)
            Output field file of the anatomical to template ANTS registration
        inputspec.ants_affine : text file
            Output matrix of the anatomical to template ANTS registration
            
    Workflow Outputs::
    
        outputspec.out_file : string (nifti file)
            Normalizion of input functional file

            
    Registration Procedure:
    
    1. Convert the FSL-based functional-to-anatomical output affine matrix into ANTS (ITK) format.
    2. Collect this converted affine, and the ants_affine.txt and nonlinear field file from
       the anatomical-to-template ANTS registration into one transformation series string.
    3. Apply the warp to the input file using WarpImageMultiTransform (for 3d files) or
       WarpTimeSeriesImageMultiTransform (for 4d files with timeseries).
       
    Workflow Graph:
    
    .. image::
        :width: 500
    
    Detailed Workflow Graph:
    
    .. image:: 
        :width: 500    
       
    """

    apply_ants_xfm = pe.Workflow(name=name)

    inputspec = pe.Node(
        util.IdentityInterface(
            fields=['warp_reference', 'in_file', 'use_nearest', 'func_anat_affine',
                'conversion_reference', 'conversion_source', 'nonlinear_field',
                'ants_affine']),
        name='inputspec')


    outputspec = pe.Node(
        util.IdentityInterface(
            fields=['out_file']),
        name='outputspec')


    if dimension == 4:

        if mapnode == 0:
 
            # converts FSL-format .mat affine xfm into ANTS-format .txt
            # .mat affine comes from Func->Anat registration
            fsl_reg_2_itk = pe.Node(interface=c3.C3dAffineTool(), name='fsl_reg_2_itk')
            fsl_reg_2_itk.inputs.itk_transform = True
            fsl_reg_2_itk.inputs.fsl2ras = True

            #collects series of transformations to be applied to the moving images
            collect_transforms = pe.Node(util.Merge(3), name='collect_transforms')

            # performs series of transformations on moving images
            warp_images = pe.Node(interface=ants.WarpTimeSeriesImageMultiTransform(), name='ants_apply_4d_warp')
            warp_images.inputs.dimension = 4

            apply_ants_xfm.connect(inputspec, 'warp_reference', warp_images, 'reference_image')
            apply_ants_xfm.connect(inputspec, 'use_nearest', warp_images, 'use_nearest')

        elif mapnode == 1:

            # in this case with 4-d images (timeseries), fsl_reg_2_itk and collect_transforms are not
            # map nodes due to the fact that the affine conversion (fsl_reg_2_itk) cannot take in a
            # 4-d image as its conversion_source (ordinarily the input image and the conversion source
            # are the same image, however with timeseries, mean_functional should be used as the
            # conversion source and the 4-d image used as the input image to the apply warp.
            # (this is why in_file and conversion_source are separated into two inputs)

            # converts FSL-format .mat affine xfm into ANTS-format .txt
            # .mat affine comes from Func->Anat registration
            fsl_reg_2_itk = pe.Node(interface=c3.C3dAffineTool(), name='fsl_reg_2_itk')
            fsl_reg_2_itk.inputs.itk_transform = True
            fsl_reg_2_itk.inputs.fsl2ras = True

            #collects series of transformations to be applied to the moving images
            collect_transforms = pe.Node(util.Merge(3), name='collect_transforms')

            # performs series of transformations on moving images
            warp_images = pe.MapNode(interface=ants.WarpTimeSeriesImageMultiTransform(), name='ants_apply_4d_warp',
                                     iterfield=['in_file'])
            warp_images.inputs.dimension = 4

            apply_ants_xfm.connect(inputspec, 'warp_reference', warp_images, 'reference_image')
            apply_ants_xfm.connect(inputspec, 'use_nearest', warp_images, 'use_nearest')

    elif dimension == 3:

        if mapnode == 0:

            # converts FSL-format .mat affine xfm into ANTS-format .txt
            # .mat affine comes from Func->Anat registration
            fsl_reg_2_itk = pe.Node(interface=c3.C3dAffineTool(), name='fsl_reg_2_itk')
            fsl_reg_2_itk.inputs.itk_transform = True
            fsl_reg_2_itk.inputs.fsl2ras = True

            #collects series of transformations to be applied to the moving images
            collect_transforms = pe.Node(util.Merge(3), name='collect_transforms')

            # performs series of transformations on moving images
            warp_images = pe.Node(interface=ants.WarpImageMultiTransform(), name='apply_ants_3d_warp')
            warp_images.inputs.dimension = 3

            apply_ants_xfm.connect(inputspec, 'warp_reference', warp_images, 'reference_image')
            apply_ants_xfm.connect(inputspec, 'use_nearest', warp_images, 'use_nearest')

        elif mapnode == 1:

            # converts FSL-format .mat affine xfm into ANTS-format .txt
            # .mat affine comes from Func->Anat registration
            fsl_reg_2_itk = pe.MapNode(interface=c3.C3dAffineTool(), name='fsl_reg_2_itk',
                                       iterfield=['source_file'])
            fsl_reg_2_itk.inputs.itk_transform = True
            fsl_reg_2_itk.inputs.fsl2ras = True

            #collects series of transformations to be applied to the moving images
            collect_transforms = pe.MapNode(util.Merge(3), name='collect_transforms',
                                            iterfield=['in3'])

            # performs series of transformations on moving images
            warp_images = pe.MapNode(interface=ants.WarpImageMultiTransform(), name='apply_ants_3d_warp',
                                     iterfield=['input_image', 'transformation_series'])
            warp_images.inputs.dimension = 3

            apply_ants_xfm.connect(inputspec, 'warp_reference', warp_images, 'reference_image')
            apply_ants_xfm.connect(inputspec, 'use_nearest', warp_images, 'use_nearest')




    # convert the .mat from linear Func->Anat to ANTS format
    apply_ants_xfm.connect(inputspec, 'func_anat_affine', fsl_reg_2_itk, 'transform_file')

    apply_ants_xfm.connect(inputspec, 'conversion_reference', fsl_reg_2_itk, 'reference_file')

    apply_ants_xfm.connect(inputspec, 'conversion_source', fsl_reg_2_itk, 'source_file')


    # Premat from Func->Anat linear reg and bbreg (if bbreg is enabled)
    apply_ants_xfm.connect(fsl_reg_2_itk, 'itk_transform', collect_transforms, 'in3')
    
    # Field file from anatomical nonlinear registration
    apply_ants_xfm.connect(inputspec, 'nonlinear_field', collect_transforms, 'in1')

    # affine transformation from anatomical registration
    apply_ants_xfm.connect(inputspec, 'ants_affine', collect_transforms, 'in2')


    apply_ants_xfm.connect(inputspec, 'in_file', warp_images, 'input_image')

    apply_ants_xfm.connect(collect_transforms, 'out', warp_images, 'transformation_series')

    apply_ants_xfm.connect(warp_images, 'output_image', outputspec, 'out_file')


    return apply_ants_xfm



def create_fsl_to_itk_conversion(name='fsl_to_itk_conversion'):

    """
    Converts an FSL-format output matrix to an ITK-format (ANTS) matrix.

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    fsl_to_itk_conversion : nipype.pipeline.engine.Workflow

    Notes
    -----
    
    Workflow Inputs::
    
        inputspec.transform_file : string (nifti file)
            Output matrix of FSL-based functional to anatomical registration
        inputspec.reference_file : string (nifti file)
            File of skull-stripped anatomical brain to be used in affine conversion
        inputspec.source_file : string (nifti file)
            Should match the input of the apply warp (in_file) unless you are
            applying the warp to a 4-d file, in which case this file should
            be a mean_functional file

            
    Workflow Outputs::
    
        outputspec.itk_transform : string (nifti file)
            Converted affine transform in ITK format usable with ANTS

            
    Registration Procedure:
    
    1. Converts an FSL-format affine matrix into an ITK-format
       affine matrix (for use with ANTS registration tools)
       
    Workflow Graph:
    
    .. image:: 
        :width: 500
    
    Detailed Workflow Graph:
    
    .. image:: 
        :width: 500      
    """

    fsl_to_itk_conversion = pe.Workflow(name=name)


    inputspec = pe.Node(
        util.IdentityInterface(
            fields=['transform_file', 'reference_file', 'source_file']),
        name='inputspec')


    # converts FSL-format .mat affine xfm into ANTS-format .txt
    # .mat affine comes from Func->Anat registration
    fsl_reg_2_itk = pe.Node(c3.C3dAffineTool(), name='fsl_reg_2_itk')
    fsl_reg_2_itk.inputs.itk_transform = True
    fsl_reg_2_itk.inputs.fsl2ras = True

    outputspec = pe.Node(
        util.IdentityInterface(
            fields=['itk_transform']),
        name='outputspec')


    fsl_to_itk_conversion.connect(inputspec, 'transform_file', fsl_reg_2_itk, 'transform_file')

    fsl_to_itk_conversion.connect(inputspec, 'reference_file', fsl_reg_2_itk, 'reference_file')

    fsl_to_itk_conversion.connect(inputspec, 'source_file', fsl_reg_2_itk, 'source_file')


    fsl_to_itk_conversion.connect(warp_brain, 'itk_transform', outputspec, 'itk_transform')



    return fsl_to_itk_conversion
