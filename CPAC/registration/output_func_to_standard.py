import nipype.interfaces.fsl as fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from registration import create_wf_apply_ants_warp, \
    create_wf_c3d_fsl_to_itk, \
    create_wf_collect_transforms

# Todo: CC distcor is not implement for fsl apply xform func to mni, why ??
def fsl_apply_transform_func_to_mni(workflow, output_name, func_key, ref_key, num_strat, strat, interpolation_method, distcor=False, map_node=False):

    strat_nodes = strat.get_nodes_names()

    if map_node == True:
        # func_mni_warp
        func_mni_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                name='func_mni_fsl_warp_{0}_{1:d}'.format(output_name, num_strat),
                iterfield=['in_file'],
                mem_gb=1.5)
    else:
        # func_mni_warp
        func_mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                name='func_mni_fsl_warp_{0}_{1:d}'.format(output_name, num_strat))

        
    func_mni_warp.inputs.interp = interpolation_method

    if func_key == 'leaf':
        func_node, func_file = strat.get_leaf_properties()
    else:
        func_node, func_file = strat[func_key]

    workflow.connect(func_node, func_file,
                     func_mni_warp, 'in_file')

    ref_node, ref_out_file = strat['template_brain_for_func_preproc']

    ref_node, ref_out_file = strat[ref_key]
    workflow.connect(ref_node, ref_out_file,
                     func_mni_warp, 'ref_file')

    if 'anat_mni_fnirt_register' in strat_nodes:

        node, out_file = strat['functional_to_anat_linear_xfm']
        workflow.connect(node, out_file,
                         func_mni_warp, 'premat')

        node, out_file = strat['anatomical_to_mni_nonlinear_xfm']
        workflow.connect(node, out_file,
                         func_mni_warp, 'field_file')

    elif 'anat_mni_flirt_register' in strat_nodes:

        if 'functional_to_mni_linear_xfm' not in strat:

            combine_transforms = pe.Node(interface=fsl.ConvertXFM(),
                name='combine_fsl_xforms_func_mni_{0}_{1:d}'.format(output_name, num_strat))

            combine_transforms.inputs.concat_xfm = True

            node, out_file = strat['anatomical_to_mni_linear_xfm']
            workflow.connect(node, out_file,
                             combine_transforms, 'in_file2')

            node, out_file = strat['functional_to_anat_linear_xfm']
            workflow.connect(node, out_file,
                             combine_transforms, 'in_file')

            strat.update_resource_pool({ 'functional_to_mni_linear_xfm': (combine_transforms, 'out_file')}) 
            strat.append_name(combine_transforms.name)

        combine_transforms, outfile = strat['functional_to_mni_linear_xfm']

        workflow.connect(combine_transforms, outfile,
                         func_mni_warp, 'premat')

        strat.update_resource_pool({ output_name: (func_mni_warp, 'out_file')}) 
        strat.append_name(func_mni_warp.name)

    else:
        raise ValueError('Could not find anat_mni_flirt_register or anat_mni_fnirt register in nodes')
    
    return workflow


# use preproc
def ants_apply_warps_func_mni(
        workflow, strat, num_strat, num_ants_cores,
        input_key, ref_key, func_name,
        interp='LanczosWindowedSinc',
        template_brain_name='template_brain_for_func_preproc',
        input_image_type=0, 
        distcor=False,
        map_node=False,
        inverse=False
    ):
    """Apply the functional-to-structural and structural-to-template warps to
    the 4D functional time-series to warp it to template space.

    Parameters
    ----------
    workflow: Nipype workflow object
        the workflow containing the resources involved
    strat: C-PAC Strategy object
        a strategy with one or more resource pools
    num_strat: int
        the number of strategy objects
    num_ants_cores: int
        the number of CPU cores dedicated to ANTS anatomical-to-standard
        registration
    input_key: string
        resource pool key correspoding to the node containing the 3D or 4D
        functional file to be written into MNI space, use 'leaf' for a 
        leaf node
    ref_key: string
        resource pool key correspoding to the node containing the reference 
        volume for the C3D FSL-to-ITK affine conversion (often the mean of 
        the functional time-series, which is a single volume)
    template_brain_name: str
        resource pool key correspoding to the file path to the template brain
        used for functional-to-template registration
    func_name: str
        what the name of the warped functional should be when written to the
        resource pool
    interp: str
        which interpolation to use when applying the warps, commonly used
        options are 'Linear', 'Bspline', 'LanczosWindowedSinc' (default) 
        for derivatives and image data 'NearestNeighbor' for masks
    input_image_type: int
        argument taken by the ANTs apply warp tool; in this case, should be
        0 for scalars (default) and 3 for 4D functional time-series
    distcor: boolean
        indicates whether a distortion correction transformation should be 
        added to the transforms, this of course requires that a distortion
        correction map exist in the resource pool
    """

    inverse_string = ''
    if inverse is True:
        inverse_string = '_inverse'

    # make sure that resource pool has some required resources before proceeding
    if 'fsl_mat_as_itk' not in strat:

        # converts FSL-format .mat affine xfm into ANTS-format
        # .txt; .mat affine comes from Func->Anat registration
        fsl_to_itk_func_mni = create_wf_c3d_fsl_to_itk(
            name='fsl_to_itk_{0}'.format(num_strat)
        )

        # convert the .mat from linear Func->Anat to
        # ANTS format
        node, out_file = strat['functional_to_anat_linear_xfm']
        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                         'inputspec.affine_file')

        node, out_file = strat["anatomical_brain"]
        workflow.connect(node, out_file, fsl_to_itk_func_mni,
                         'inputspec.reference_file')

        ref_node, ref_out = strat[ref_key]
        workflow.connect(ref_node, ref_out,
                            fsl_to_itk_func_mni,
                            'inputspec.source_file')

        strat.update_resource_pool({
            'fsl_mat_as_itk': (fsl_to_itk_func_mni, 'outputspec.itk_transform')
        })

        strat.append_name(fsl_to_itk_func_mni.name)

    # stack of transforms to be combined to acheive the desired transformation
    collect_transforms_key = 'collect_transforms_func_mni{0}{1}'.format('_distcor' if distcor else '', inverse_string)

    if collect_transforms_key not in strat:

        # collects series of warps to be applied
        collect_transforms_func_mni = \
            create_wf_collect_transforms(
                inverse=inverse,
                name='collect_transforms_func_mni{0}_{1}'.format(inverse_string, num_strat)
            )

        # transforms to be concatenated, the first element of each tuple is the resource
        # pool key related to the resource that should be connected in, and the second
        # element is the input to which it should be connected
        transforms_to_combine = [('anatomical_to_mni_nonlinear_xfm', 'inputspec.warp_file'),
                                 ('ants_initial_xfm', 'inputspec.linear_initial'),
                                 ('ants_affine_xfm', 'inputspec.linear_affine'),
                                 ('ants_rigid_xfm', 'inputspec.linear_rigid'),
                                 ('fsl_mat_as_itk', 'inputspec.fsl_to_itk_affine')]

        if distcor is True:
            if inverse is True:
                transforms_to_combine.append(('blip_warp_inverse', 'inputspec.distortion_unwarp'))
            else:
                transforms_to_combine.append(('blip_warp', 'inputspec.distortion_unwarp'))

        for transform_key, inputspec_port in transforms_to_combine:

             # Field file from anatomical nonlinear registration
             node, out_file = strat[transform_key]
             workflow.connect(node, out_file,
                 collect_transforms_func_mni,
                 inputspec_port)

        strat.update_resource_pool({
            collect_transforms_key: (collect_transforms_func_mni, 'outputspec.transformation_series')
        })

        strat.append_name(collect_transforms_func_mni.name)

    #### now we add in the apply ants warps node
    apply_ants_warp_func_mni = \
        create_wf_apply_ants_warp(map_node, name='apply_ants_warp_{0}{1}_{2}'.format(func_name, inverse_string, num_strat),
                                  ants_threads=int(num_ants_cores), inverse=inverse)

    # input_image_type:
    # (0 or 1 or 2 or 3)
    # Option specifying the input image type of scalar
    # (default), vector, tensor, or time series.
    apply_ants_warp_func_mni.inputs.inputspec.input_image_type = input_image_type
    apply_ants_warp_func_mni.inputs.inputspec.dimension = 3
    apply_ants_warp_func_mni.inputs.inputspec.interpolation = interp

    node, out_file = strat[template_brain_name]
    workflow.connect(node, out_file,
            apply_ants_warp_func_mni, 'inputspec.reference_image')

    collect_node, collect_out = strat[collect_transforms_key]
    workflow.connect(collect_node, collect_out,
                     apply_ants_warp_func_mni, 'inputspec.transforms')

    if input_key == "leaf":
        input_node, input_out = strat.get_leaf_properties()
    else:
        input_node, input_out = strat[input_key]

    workflow.connect(input_node, input_out,
                     apply_ants_warp_func_mni, 'inputspec.input_image')

    strat.update_resource_pool({
        func_name: (apply_ants_warp_func_mni, 'outputspec.output_image')
    })

    strat.append_name(apply_ants_warp_func_mni.name)

    return workflow

def output_func_to_standard(workflow, func_key, ref_key, output_name, strat, num_strat, pipeline_config_obj,
                            input_image_type='derivative', inverse=True):

    image_types = ['func_derivative', 'func_derivative_multi', 'func_4d', 'func_mask']

    if input_image_type not in image_types:
        raise ValueError('Input image type {0} should be one of {1}'.format(input_image_type, ', '.join(image_types)))

    nodes = strat.get_nodes_names()

    map_node = True if input_image_type == 'derivative_multi' else False

    distcor = True if 'epi_distcorr' in nodes or 'blip_correct' in nodes else False

    if 'anat_mni_fnirt_register' in nodes or 'anat_mni_flirt_register' in nodes:

        if input_image_type == 'map':
            interp = 'nn'
        else:
            interp = pipeline_config_obj.funcRegFSLinterpolation

        fsl_apply_transform_func_to_mni(workflow, output_name, func_key, ref_key, num_strat,
                strat, interp, distcor=distcor, map_node=map_node)

    elif 'ANTS' in pipeline_config_obj.regOption:

        if input_image_type == 'map':
            interp = 'nn'
        else:
            interp = pipeline_config_obj.funcRegANTSinterpolation

        image_type = 1 if input_image_type == 'func4d' else 0

        ants_apply_warps_func_mni(workflow, strat, num_strat,
                            pipeline_config_obj.num_ants_threads,
                            func_key, 'mean_functional', output_name,
                            interp=interp,
                            template_brain_name=ref_key,
                            input_image_type=image_type,
                            distcor=distcor,
                            map_node=map_node,
                            inverse=inverse
                        )

    return strat


