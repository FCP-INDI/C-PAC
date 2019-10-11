

import nipype.interfaces.fsl as fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.afni import preprocess
from CPAC.registration import create_fsl_fnirt_nonlinear_reg, \
    create_register_func_to_anat, \
    create_bbregister_func_to_anat, \
    create_wf_calculate_ants_warp, \
    create_wf_apply_ants_warp, \
    create_wf_c3d_fsl_to_itk, \
    create_wf_collect_transforms

from CPAC.utils import Configuration, function, find_files
from CPAC.utils.utils import (
    set_gauss,
    get_scan_params,
    extract_output_mean,
    get_zscore,
    get_fisher_zscore)


# use preproc
def ants_apply_warps_func_mni(
        workflow, strat, num_strat, num_ants_cores,
        input_key, ref_key, func_name,
        interp='LanczosWindowedSinc',
        template_brain_name='template_brain_for_func_preproc',
        input_image_type=0, 
        distcor=False
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
    collect_transforms_key = 'collect_transforms_func_mni{0}'.format('_distcor' if distcor else '')

    if collect_transforms_key not in strat:

        # collects series of warps to be applied
        collect_transforms_func_mni = \
            create_wf_collect_transforms(
                name='collect_transforms_func_mni_{0}'.format(num_strat)
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
        create_wf_apply_ants_warp(name='apply_ants_warp_{0}_{1}'.format(func_name, num_strat),
                                  ants_threads=int(num_ants_cores))

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
