import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
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
                                                       'ref_mask',
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

    nonlinear_register.connect(inputspec, 'ref_mask',
                               nonlinear_reg, 'refmask_file')
    
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
    
    """
    Registers a functional scan in native space to anatomical space using a
    linear transform and does not include bbregister.

    Parameters
    ----------
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
    
    


def create_wf_calculate_ants_warp(name='create_wf_calculate_ants_warp', mult_input=0):

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
    
        inputspec.anatomical_brain : string (nifti file)
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

    Workflow Graph:
    
    .. image:: 
        :width: 500

    Detailed Workflow Graph:
    
    .. image:: 
        :width: 500      
    '''

    import nipype.interfaces.ants as ants
    from nipype.interfaces.utility import Function
    from CPAC.registration.utils import seperate_warps_list, \
                                        combine_inputs_into_list, \
                                        hardcoded_reg


    calc_ants_warp_wf = pe.Workflow(name=name)


    inputspec = pe.Node(util.IdentityInterface(fields=['anatomical_brain',
            'reference_brain', 'dimension', 'use_histogram_matching',
            'winsorize_lower_quantile', 'winsorize_upper_quantile', 'metric',
            'metric_weight', 'radius_or_number_of_bins', 'sampling_strategy',
            'sampling_percentage', 'number_of_iterations', 
            'convergence_threshold', 'convergence_window_size', 'transforms',
            'transform_parameters', 'shrink_factors', 'smoothing_sigmas',
            'write_composite_transform', 'anatomical_skull',
            'reference_skull']), name='inputspec')#,'wait']),name='inputspec')


    # use ANTS to warp the masked anatomical image to a template image
    '''
    calculate_ants_warp = pe.Node(interface=ants.Registration(),
            name='calculate_ants_warp')

    calculate_ants_warp.inputs.output_warped_image = True
    calculate_ants_warp.inputs.initial_moving_transform_com = 0
    '''
    calculate_ants_warp = pe.Node(interface=util.Function(input_names=['anatomical_brain', 'reference_brain', 'anatomical_skull', 'reference_skull', 'wait'], output_names=['warp_list', 'warped_image'], function=hardcoded_reg), name='calc_ants_warp')

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

    select_forward_warp.inputs.selection = "3Warp"


    select_inverse_warp = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=seperate_warps_list), name='select_inverse_warp')

    select_inverse_warp.inputs.selection = "Inverse"


    outputspec = pe.Node(util.IdentityInterface(fields=['ants_initial_xfm',
            'ants_rigid_xfm', 'ants_affine_xfm', 'warp_field',
            'inverse_warp_field', 'composite_transform', 'wait',
            'normalized_output_brain']), name='outputspec')


    # connections from inputspec

    if mult_input == 1:

        '''
        combine_inputs = pe.Node(util.Function(input_names=['input1', 'input2', 'input3'],
                output_names=['inputs_list'], function=combine_inputs_into_list),
                name='ants_reg_combine_inputs')

        combine_refs = pe.Node(util.Function(input_names=['input1', 'input2', 'input3'],
                output_names=['inputs_list'], function=combine_inputs_into_list),
                name='ants_reg_combine_refs')
        '''

        calc_ants_warp_wf.connect(inputspec, 'anatomical_brain',
                calculate_ants_warp, 'anatomical_brain')

        calc_ants_warp_wf.connect(inputspec, 'anatomical_skull',
                calculate_ants_warp, 'anatomical_skull')

        calc_ants_warp_wf.connect(inputspec, 'reference_brain',
                calculate_ants_warp, 'reference_brain')

        calc_ants_warp_wf.connect(inputspec, 'reference_skull',
                calculate_ants_warp, 'reference_skull')

        #calc_ants_warp_wf.connect(inputspec, 'wait',
        #        calculate_ants_warp, 'wait')

        '''
        calc_ants_warp_wf.connect(inputspec, 'anatomical_brain',
                combine_inputs, 'input1')

        calc_ants_warp_wf.connect(inputspec, 'anatomical_brain',
                combine_inputs, 'input2')

        calc_ants_warp_wf.connect(inputspec, 'anatomical_skull',
                combine_inputs, 'input3')

        calc_ants_warp_wf.connect(combine_inputs, 'inputs_list',
                calculate_ants_warp, 'moving_image')

        calc_ants_warp_wf.connect(inputspec, 'reference_brain',
                combine_refs, 'input1')

        calc_ants_warp_wf.connect(inputspec, 'reference_brain',
                combine_refs, 'input2')

        calc_ants_warp_wf.connect(inputspec, 'reference_skull',
                combine_refs, 'input3')

        calc_ants_warp_wf.connect(combine_refs, 'inputs_list',
                calculate_ants_warp, 'fixed_image') 
        '''

    else:

        '''
        calc_ants_warp_wf.connect(inputspec, 'anatomical_brain',
                calculate_ants_warp, 'moving_image')

        calc_ants_warp_wf.connect(inputspec, 'reference_brain',
                calculate_ants_warp, 'fixed_image')
        '''

        calc_ants_warp_wf.connect(inputspec, 'anatomical_brain',
                calculate_ants_warp, 'anatomical_brain')

        calc_ants_warp_wf.connect(inputspec, 'anatomical_brain',
                calculate_ants_warp, 'anatomical_skull')

        calc_ants_warp_wf.connect(inputspec, 'reference_brain',
                calculate_ants_warp, 'reference_brain')

        calc_ants_warp_wf.connect(inputspec, 'reference_brain',
                calculate_ants_warp, 'reference_skull')

        #calc_ants_warp_wf.connect(inputspec, 'wait',
        #        calculate_ants_warp, 'wait')


    calc_ants_warp_wf.connect(inputspec, 'dimension', calculate_ants_warp,
            'dimension')

    calc_ants_warp_wf.connect(inputspec, 'use_histogram_matching',
            calculate_ants_warp, 'use_histogram_matching')

    calc_ants_warp_wf.connect(inputspec, 'winsorize_lower_quantile',
            calculate_ants_warp, 'winsorize_lower_quantile')

    calc_ants_warp_wf.connect(inputspec, 'winsorize_upper_quantile',
            calculate_ants_warp, 'winsorize_upper_quantile')

    calc_ants_warp_wf.connect(inputspec, 'metric', calculate_ants_warp,
            'metric')

    calc_ants_warp_wf.connect(inputspec, 'metric_weight', calculate_ants_warp,
            'metric_weight')

    calc_ants_warp_wf.connect(inputspec, 'radius_or_number_of_bins',
            calculate_ants_warp, 'radius_or_number_of_bins')

    calc_ants_warp_wf.connect(inputspec, 'sampling_strategy',
            calculate_ants_warp, 'sampling_strategy')

    calc_ants_warp_wf.connect(inputspec, 'sampling_percentage',
            calculate_ants_warp, 'sampling_percentage')

    calc_ants_warp_wf.connect(inputspec, 'number_of_iterations',
            calculate_ants_warp, 'number_of_iterations')

    calc_ants_warp_wf.connect(inputspec, 'convergence_threshold',
            calculate_ants_warp, 'convergence_threshold')

    calc_ants_warp_wf.connect(inputspec, 'convergence_window_size',
            calculate_ants_warp, 'convergence_window_size')

    calc_ants_warp_wf.connect(inputspec, 'transforms', calculate_ants_warp,
            'transforms')

    calc_ants_warp_wf.connect(inputspec, 'transform_parameters',
            calculate_ants_warp, 'transform_parameters')

    calc_ants_warp_wf.connect(inputspec, 'shrink_factors',
            calculate_ants_warp, 'shrink_factors')

    calc_ants_warp_wf.connect(inputspec, 'smoothing_sigmas',
            calculate_ants_warp, 'smoothing_sigmas')

    calc_ants_warp_wf.connect(inputspec, 'write_composite_transform',
            calculate_ants_warp, 'write_composite_transform')

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

#    calc_ants_warp_wf.connect(inputspec, 'wait',
#            outputspec, 'wait')


    return calc_ants_warp_wf



def create_wf_apply_ants_warp(map_node, name='create_wf_apply_ants_warp'):

    """
    Applies previously calculated ANTS registration transforms to input
    images. This workflow employs the antsApplyTransforms tool:

    http://stnava.github.io/ANTs/

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    apply_ants_warp_wf : nipype.pipeline.engine.Workflow

    Notes
    -----

    Workflow Inputs::

        inputspec.input_image : string (nifti file)
            Image file of brain to be registered to reference
        inputspec.reference_image : string (nifti file)
            Image file of brain or template being used as a reference
        inputspec.transforms : list of filepaths (nifti, .mat, .txt)
            List of transforms and warps to be applied to the input image
        inputspec.dimension : integer
            Dimension value of image being registered (2, 3, or 4)
        inputspec.interpolation : string
            Type of interpolation to be used. See antsApplyTransforms
            documentation or Nipype interface documentation for options

            
    Workflow Outputs::
    
        outputspec.output_image : string (nifti file)
            Normalized output file

                 
    Workflow Graph:
    
    .. image::
        :width: 500
    
    Detailed Workflow Graph:
    
    .. image:: 
        :width: 500
       
    """

    import nipype.interfaces.ants as ants

    apply_ants_warp_wf = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['input_image', 
            'reference_image', 'transforms', 'dimension', 'input_image_type', 
            'interpolation']), name='inputspec')

    if map_node == 0:
        apply_ants_warp = pe.Node(interface=ants.ApplyTransforms(),
                name='apply_ants_warp')

    elif map_node == 1:
        apply_ants_warp = pe.MapNode(interface=ants.ApplyTransforms(),
                name='apply_ants_warp_mapnode', iterfield=['input_image', \
                'transforms'])

    apply_ants_warp.inputs.out_postfix = '_antswarp'

    outputspec = pe.Node(util.IdentityInterface(fields=['output_image']),
            name='outputspec')


    # connections from inputspec

    apply_ants_warp_wf.connect(inputspec, 'input_image', apply_ants_warp, 
            'input_image')

    apply_ants_warp_wf.connect(inputspec, 'reference_image', apply_ants_warp, 
            'reference_image')

    apply_ants_warp_wf.connect(inputspec, 'transforms', apply_ants_warp, 
            'transforms')

    apply_ants_warp_wf.connect(inputspec, 'dimension', apply_ants_warp, 
            'dimension')

    apply_ants_warp_wf.connect(inputspec, 'input_image_type', apply_ants_warp, 
            'input_image_type')

    apply_ants_warp_wf.connect(inputspec, 'interpolation', apply_ants_warp, 
            'interpolation')

    # connections to outputspec

    apply_ants_warp_wf.connect(apply_ants_warp, 'output_image',
            outputspec, 'output_image')


    return apply_ants_warp_wf


def create_wf_c3d_fsl_to_itk(map_node, input_image_type=0, name='create_wf_c3d_fsl_to_itk'):

    """
    Converts an FSL-format output matrix to an ITK-format (ANTS) matrix
    for use with ANTS registration tools.

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
    
        inputspec.affine_file : string (nifti file)
            Output matrix of FSL-based functional to anatomical registration
        inputspec.reference_file : string (nifti file)
            File of skull-stripped anatomical brain to be used in affine
            conversion
        inputspec.source_file : string (nifti file)
            Should match the input of the apply warp (in_file) unless you are
            applying the warp to a 4-d file, in which case this file should
            be a mean_functional file

    Workflow Outputs::
    
        outputspec.itk_transform : string (nifti file)
            Converted affine transform in ITK format usable with ANTS
    
    """

    import nipype.interfaces.c3 as c3
    from nipype.interfaces.utility import Function
    from CPAC.registration.utils import change_itk_transform_type
    from nipype.interfaces.afni import preprocess

    fsl_to_itk_conversion = pe.Workflow(name=name)


    inputspec = pe.Node(util.IdentityInterface(fields=['affine_file',
            'reference_file', 'source_file']), name='inputspec')


    # converts FSL-format .mat affine xfm into ANTS-format .txt
    # .mat affine comes from Func->Anat registration

    if map_node == 0:
        fsl_reg_2_itk = pe.Node(c3.C3dAffineTool(), name='fsl_reg_2_itk')

    elif map_node == 1:
        fsl_reg_2_itk = pe.MapNode(c3.C3dAffineTool(),
                name='fsl_reg_2_itk_mapnode', iterfield=['source_file'])
        
    fsl_reg_2_itk.inputs.itk_transform = True
    fsl_reg_2_itk.inputs.fsl2ras = True

    if map_node == 0:
        change_transform = pe.Node(util.Function(\
                input_names=['input_affine_file'],
                output_names=['updated_affine_file'], 
                function=change_itk_transform_type),
                name='change_transform_type')

    elif map_node == 1:
        change_transform = pe.MapNode(util.Function(\
                input_names=['input_affine_file'],
                output_names=['updated_affine_file'], 
                function=change_itk_transform_type),
                name='change_transform_type', iterfield=['input_affine_file'])


    outputspec = pe.Node(util.IdentityInterface(fields=['itk_transform']),
            name='outputspec')


    fsl_to_itk_conversion.connect(inputspec, 'affine_file', fsl_reg_2_itk,
            'transform_file')

    fsl_to_itk_conversion.connect(inputspec, 'reference_file', fsl_reg_2_itk,
            'reference_file')

    # source_file input of the conversion must be a 3D file, so if the source
    # file is 4D (input_image_type=3), average it into a 3D file first
    if input_image_type == 0:

        fsl_to_itk_conversion.connect(inputspec, 'source_file', fsl_reg_2_itk,
                'source_file')

    elif input_image_type == 3:

        tstat_source = pe.Node(interface=preprocess.TStat(),
                name='fsl_to_itk_tcat_source')
        tstat_source.inputs.outputtype = 'NIFTI_GZ'
        tstat_source.inputs.options = '-mean'

        fsl_to_itk_conversion.connect(inputspec, 'source_file', tstat_source,
                'in_file')

        fsl_to_itk_conversion.connect(tstat_source, 'out_file', fsl_reg_2_itk,
                'source_file')


    fsl_to_itk_conversion.connect(fsl_reg_2_itk, 'itk_transform',
            change_transform, 'input_affine_file')

    fsl_to_itk_conversion.connect(change_transform, 'updated_affine_file', 
            outputspec, 'itk_transform')


    return fsl_to_itk_conversion



def create_wf_collect_transforms(map_node, name='create_wf_collect_transforms'):

    """
    DOCSTRINGS

    Parameters
    ----------
    name : string, optional
        Name of the workflow.

    Returns
    -------
    collect_transforms_wf : nipype.pipeline.engine.Workflow

    Notes
    -----
    
    Workflow Inputs::
    
        inputspec.transform_file : string (nifti file)
            Output matrix of FSL-based functional to anatomical registration
        inputspec.reference_file : string (nifti file)
            File of skull-stripped anatomical brain to be used in affine
            conversion
        inputspec.source_file : string (nifti file)
            Should match the input of the apply warp (in_file) unless you are
            applying the warp to a 4-d file, in which case this file should
            be a mean_functional file

    Workflow Outputs::
    
        outputspec.itk_transform : string (nifti file)
            Converted affine transform in ITK format usable with ANTS
    
    """

    collect_transforms_wf = pe.Workflow(name=name)


    inputspec = pe.Node(util.IdentityInterface(fields=['warp_file',
            'linear_initial', 'linear_affine', 'linear_rigid', \
            'fsl_to_itk_affine']), name='inputspec')


    # converts FSL-format .mat affine xfm into ANTS-format .txt
    # .mat affine comes from Func->Anat registration

    if map_node == 0:
        collect_transforms = pe.Node(util.Merge(5), name='collect_transforms')

    elif map_node == 1:
        collect_transforms = pe.MapNode(util.Merge(5),
                name='collect_transforms_mapnode', iterfield=['in5'])

    outputspec = pe.Node(util.IdentityInterface(
            fields=['transformation_series']), name='outputspec')


    # Field file from anatomical nonlinear registration
    collect_transforms_wf.connect(inputspec, 'warp_file', collect_transforms,
            'in1')

    # affine transformation from anatomical registration
    collect_transforms_wf.connect(inputspec, 'linear_affine',
            collect_transforms, 'in2')

    # rigid transformation from anatomical registration
    collect_transforms_wf.connect(inputspec, 'linear_rigid',
            collect_transforms, 'in3')

    # initial transformation from anatomical registration
    collect_transforms_wf.connect(inputspec, 'linear_initial',
            collect_transforms, 'in4')

    # Premat from Func->Anat linear reg and bbreg (if bbreg is enabled)
    collect_transforms_wf.connect(inputspec, 'fsl_to_itk_affine',
            collect_transforms, 'in5')

    collect_transforms_wf.connect(collect_transforms, 'out', outputspec,
            'transformation_series')


    return collect_transforms_wf



