
import nipype.interfaces.fsl as fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.afni import preprocess
from CPAC.registration import create_nonlinear_register, \
    create_register_func_to_anat, \
    create_bbregister_func_to_anat, \
    create_wf_calculate_ants_warp, \
    create_wf_apply_ants_warp, \
    create_wf_c3d_fsl_to_itk, \
    create_wf_collect_transforms

from CPAC.utils import Configuration, function, find_files
from CPAC.utils.utils import extract_one_d, set_gauss, \
    process_outputs, get_scan_params, \
    get_tr, extract_txt, create_log, \
    extract_output_mean, create_output_mean_csv, get_zscore, \
    get_fisher_zscore, dbg_file_lineno, add_afni_prefix

# Apply warps, Z-scoring, Smoothing, Averages

def output_to_standard(workflow, output_name, strat, num_strat, pipeline_config_obj,
                        map_node=False, input_image_type=0):

    nodes = strat.get_nodes_names()

    if 'apply_ants_warp_functional_to_standard' in nodes:

        # ANTS WARP APPLICATION

        # convert the func-to-anat linear warp from FSL FLIRT to
        # ITK (ANTS) format
        fsl_to_itk_convert = create_wf_c3d_fsl_to_itk(input_image_type,
                                                      map_node,
                                                      name='{0}_fsl_to_itk_{1}'.format(output_name, num_strat))

        # collect the list of warps into a single stack to feed into the
        # ANTS warp apply tool
        collect_transforms = create_wf_collect_transforms(map_node,
                                                          name='{0}_collect_transforms_{1}'.format(output_name, num_strat))

        # ANTS apply warp
        apply_ants_warp = create_wf_apply_ants_warp(map_node,
                                                    name='{0}_to_standard_{1}'.format(
                                                        output_name, num_strat),
                                                    ants_threads=int(pipeline_config_obj.num_ants_threads))

        apply_ants_warp.inputs.inputspec.dimension = 3
        apply_ants_warp.inputs.inputspec.interpolation = 'Linear'
        apply_ants_warp.inputs.inputspec.reference_image = \
            pipeline_config_obj.template_brain_only_for_func

        apply_ants_warp.inputs.inputspec.input_image_type = \
            input_image_type

        # affine from FLIRT func->anat linear registration
        node, out_file = strat['functional_to_anat_linear_xfm']
        workflow.connect(node, out_file, fsl_to_itk_convert,
                            'inputspec.affine_file')

        # reference used in FLIRT func->anat linear registration
        node, out_file = strat['anatomical_brain']
        workflow.connect(node, out_file, fsl_to_itk_convert,
                            'inputspec.reference_file')

        # output file to be converted
        node, out_file = \
            strat[output_name]
        workflow.connect(node, out_file, fsl_to_itk_convert,
                            'inputspec.source_file')

        # nonlinear warp from anatomical->template ANTS registration
        node, out_file = strat['anatomical_to_mni_nonlinear_xfm']
        workflow.connect(node, out_file, collect_transforms,
                            'inputspec.warp_file')

        # linear initial from anatomical->template ANTS registration
        node, out_file = strat['ants_initial_xfm']
        workflow.connect(node, out_file, collect_transforms,
                            'inputspec.linear_initial')

        # linear affine from anatomical->template ANTS registration
        node, out_file = strat['ants_affine_xfm']
        workflow.connect(node, out_file, collect_transforms,
                            'inputspec.linear_affine')

        # rigid affine from anatomical->template ANTS registration
        node, out_file = strat['ants_rigid_xfm']
        workflow.connect(node, out_file, collect_transforms,
                            'inputspec.linear_rigid')

        # converted FLIRT func->anat affine, now in ITK (ANTS) format
        workflow.connect(fsl_to_itk_convert,
                            'outputspec.itk_transform',
                            collect_transforms,
                            'inputspec.fsl_to_itk_affine')

        # output file to be converted
        node, out_file = strat[output_name]
        workflow.connect(node, out_file, apply_ants_warp,
                            'inputspec.input_image')

        # collection of warps to be applied to the output file
        workflow.connect(collect_transforms,
                            'outputspec.transformation_series',
                            apply_ants_warp,
                            'inputspec.transforms')

        strat.update_resource_pool({
            '{0}_to_standard'.format(output_name): (apply_ants_warp, 'outputspec.output_image')
        })

        strat.append_name(apply_ants_warp.name)

        num_strat += 1

    else:
        # FSL WARP APPLICATION
        if map_node:
            apply_fsl_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                                        name='{0}_to_standard_{1}'.format(output_name, num_strat),
                                        iterfield=['in_file'])
        else:
            apply_fsl_warp = pe.Node(interface=fsl.ApplyWarp(),
                                        name='{0}_to_standard_{1}'.format(output_name,
                                                                        num_strat))

        apply_fsl_warp.inputs.ref_file = \
            pipeline_config_obj.template_skull_for_func

        # output file to be warped
        node, out_file = strat[output_name]
        workflow.connect(node, out_file, apply_fsl_warp, 'in_file')

        # linear affine from func->anat linear FLIRT registration
        node, out_file = strat['functional_to_anat_linear_xfm']
        workflow.connect(node, out_file, apply_fsl_warp, 'premat')

        # nonlinear warp from anatomical->template FNIRT registration
        node, out_file = strat['anatomical_to_mni_nonlinear_xfm']
        workflow.connect(node, out_file, apply_fsl_warp, 'field_file')

        strat.update_resource_pool({'{0}_to_standard'.format(output_name): (apply_fsl_warp, 'out_file')})
        strat.append_name(apply_fsl_warp.name)

    return strat


def z_score_standardize(workflow, output_name, mask_name,
                        strat, num_strat, map_node=False):

    # call the z-scoring sub-workflow builder
    z_score_std = get_zscore(output_name, map_node,
                             'z_score_std_%s_%d' % (output_name, num_strat))

    node, out_file = strat[output_name]
    workflow.connect(node, out_file,
                     z_score_std, 'inputspec.input_file')

    # get the mask
    if type(mask_name) == str:
        node, out_file = strat[mask_name]
        workflow.connect(node, out_file,
                         z_score_std, 'inputspec.mask_file')
    else:
        # mask_name is a direct file path and not the name of a
        # resource pool key
        workflow.connect(mask_name, 'local_path',
                         z_score_std, 'inputspec.mask_file')

    strat.append_name(z_score_std.name)
    strat.update_resource_pool({'{0}_zstd'.format(output_name): (z_score_std, 'outputspec.z_score_img')})

    return strat


def fisher_z_score_standardize(workflow, output_name, timeseries_oned_file,
                               strat, num_strat, map_node=False):

    # call the fisher r-to-z sub-workflow builder
    fisher_z_score_std = get_fisher_zscore(output_name, map_node,
                                            'fisher_z_score_std_%s_%d' \
                                            % (output_name, num_strat))

    node, out_file = strat[output_name]

    workflow.connect(node, out_file, fisher_z_score_std,
                        'inputspec.correlation_file')

    node, out_file = strat[timeseries_oned_file]
    workflow.connect(node, out_file, fisher_z_score_std,
                        'inputspec.timeseries_one_d')

    strat.append_name(fisher_z_score_std.name)
    strat.update_resource_pool({'{0}_fisher_zstd'.format(output_name): (fisher_z_score_std, 'outputspec.fisher_z_score_img')})

    return strat


def output_smooth(workflow, output_name, mask_name, fwhm,
                  strat, num_strat, map_node=False):

    if map_node:
        output_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                   name='{0}_smooth_{1}'.format(output_name,
                                                                num_strat),
                                   iterfield=['in_file'])
    else:
        output_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                                name='{0}_smooth_{1}'.format(output_name,
                                                             num_strat))

    # TODO review connetion to config, is the node really necessary?
    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input_{0}_{1}'.format(output_name, num_strat))
    inputnode_fwhm.iterables = ("fwhm", fwhm)

    # get the resource to be smoothed
    node, out_file = strat[output_name]

    workflow.connect(node, out_file, output_smooth, 'in_file')

    # get the parameters for fwhm
    workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                     output_smooth, 'op_string')

    # get the mask
    if type(mask_name) == str:
        node, out_file = strat[mask_name]
        workflow.connect(node, out_file,
                         output_smooth, 'operand_files')
    else:
        # mask_name is a direct file path and not the name of a
        # resource pool key
        workflow.connect(mask_name, 'local_path',
                         output_smooth, 'operand_files')

    strat.append_name(output_smooth.name)
    strat.update_resource_pool({'{0}_smooth'.format(output_name): (output_smooth, 'out_file')})

    return strat


def calc_avg(workflow, output_name, strat, num_strat, map_node=False):
    """Calculate the average of an output using AFNI 3dmaskave."""

    if map_node:
        calc_average = pe.MapNode(interface=preprocess.Maskave(),
                                  name='{0}_mean_{1}'.format(output_name,
                                                             num_strat),
                                  iterfield=['in_file'])

        mean_to_csv = pe.MapNode(function.Function(input_names=['in_file',
                                                                'output_name'],
                                                   output_names=[
                                                       'output_mean'],
                                                   function=extract_output_mean,
                                                   as_module=True),
                                 name='{0}_mean_to_txt_{1}'.format(output_name,
                                                                   num_strat),
                                 iterfield=['in_file'])
    else:
        calc_average = pe.Node(interface=preprocess.Maskave(),
                               name='{0}_mean_{1}'.format(output_name,
                                                          num_strat))

        mean_to_csv = pe.Node(function.Function(input_names=['in_file',
                                                             'output_name'],
                                                output_names=['output_mean'],
                                                function=extract_output_mean,
                                                as_module=True),
                              name='{0}_mean_to_txt_{1}'.format(output_name,
                                                                num_strat))

    mean_to_csv.inputs.output_name = output_name

    node, out_file = strat[output_name]
    workflow.connect(node, out_file, calc_average, 'in_file')
    workflow.connect(calc_average, 'out_file', mean_to_csv, 'in_file')

    strat.append_name(calc_average.name)
    strat.update_resource_pool({
        'output_means.@{0}_average'.format(output_name): (mean_to_csv, 'output_mean')
    })

    return strat


def ants_apply_warps_func_mni(
        workflow, strat, num_strat, num_ants_cores,
        input_node, input_outfile,
        ref_node, ref_outfile, standard,
        func_name, interp,
        input_image_type
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
    input_node: Nipype pointer
        pointer to the node containing the 4D functional time-series (often
        the leaf node)
    input_outfile: Nipype pointer
        pointer to the output of the node, i.e. the 4D functional time-series
        itself
    ref_node: Nipype pointer
        pointer to the node containing the reference volume for the C3D
        FSL-to-ITK affine conversion (often the mean of the functional
        time-series, which is a single volume)
    ref_outfile: Nipype pointer
        pointer to the output of ref_node, i.e. the reference volume itself
    standard: str
        file path to the template brain used for functional-to-template
        registration
    func_name: str
        what the name of the warped functional should be when written to the
        resource pool
    interp: str
        which interpolation to use when applying the warps
    input_image_type: int
        argument taken by the ANTs apply warp tool; in this case, should be
        3 for 4D functional time-series
    """

    # converts FSL-format .mat affine xfm into ANTS-format
    # .txt; .mat affine comes from Func->Anat registration
    fsl_to_itk_func_mni = create_wf_c3d_fsl_to_itk(
        name='fsl_to_itk_%s_%d' % (func_name, num_strat)
    )

    # collects series of warps to be applied
    collect_transforms_func_mni = \
        create_wf_collect_transforms(
            name='collect_transforms_%s_%d' % (func_name, num_strat)
        )

    # apply ants warps
    apply_ants_warp_func_mni = \
        create_wf_apply_ants_warp(name='apply_ants_warp_%s_%d' % (func_name, num_strat),
                                  ants_threads=int(num_ants_cores))

    apply_ants_warp_func_mni.inputs.inputspec.reference_image = standard
    apply_ants_warp_func_mni.inputs.inputspec.dimension = 3
    apply_ants_warp_func_mni.inputs.inputspec.interpolation = interp

    # input_image_type:
    # (0 or 1 or 2 or 3)
    # Option specifying the input image type of scalar
    # (default), vector, tensor, or time series.
    apply_ants_warp_func_mni.inputs.inputspec. \
        input_image_type = input_image_type

    # convert the .mat from linear Func->Anat to
    # ANTS format
    node, out_file = strat['functional_to_anat_linear_xfm']
    workflow.connect(node, out_file, fsl_to_itk_func_mni,
                     'inputspec.affine_file')

    node, out_file = strat["anatomical_brain"]
    workflow.connect(node, out_file, fsl_to_itk_func_mni,
                     'inputspec.reference_file')

    workflow.connect(ref_node, ref_outfile,
                     fsl_to_itk_func_mni,
                     'inputspec.source_file')

    # Field file from anatomical nonlinear registration
    node, out_file = strat['anatomical_to_mni_nonlinear_xfm']
    workflow.connect(node, out_file,
                     collect_transforms_func_mni,
                     'inputspec.warp_file')

    # initial transformation from anatomical registration
    node, out_file = strat['ants_initial_xfm']
    workflow.connect(node, out_file,
                     collect_transforms_func_mni,
                     'inputspec.linear_initial')

    # affine transformation from anatomical registration
    node, out_file = strat['ants_affine_xfm']
    workflow.connect(node, out_file,
                     collect_transforms_func_mni,
                     'inputspec.linear_affine')

    # rigid transformation from anatomical registration
    node, out_file = strat['ants_rigid_xfm']
    workflow.connect(node, out_file,
                     collect_transforms_func_mni,
                     'inputspec.linear_rigid')

    # Premat from Func->Anat linear reg and bbreg
    # (if bbreg is enabled)
    workflow.connect(fsl_to_itk_func_mni,
                     'outputspec.itk_transform',
                     collect_transforms_func_mni,
                     'inputspec.fsl_to_itk_affine')

    # this <node, out_file> pulls in directly because
    # it pulls in the leaf in some instances
    workflow.connect(input_node,
                     input_outfile,
                     apply_ants_warp_func_mni,
                     'inputspec.input_image')

    workflow.connect(collect_transforms_func_mni,
                     'outputspec.transformation_series',
                     apply_ants_warp_func_mni,
                     'inputspec.transforms')

    strat.update_resource_pool({
        func_name: (apply_ants_warp_func_mni, 'outputspec.output_image')
    })

    strat.append_name(apply_ants_warp_func_mni.name)

    return apply_ants_warp_func_mni


def ants_apply_inverse_warps_template_to_func(
        workflow, strat, num_strat, num_ants_cores, input_node, input_outfile,
        ref_node, ref_outfile, func_name, interp, input_image_type
):
    """Apply the functional-to-structural and structural-to-template warps
    inversely to functional time-series in template space to warp it back to
    native functional space.

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
    input_node: Nipype pointer
        pointer to the node containing the 4D functional time-series (often
        the leaf node)
    input_outfile: Nipype pointer
        pointer to the output of the node, i.e. the 4D functional time-series
        itself
    ref_node: Nipype pointer
        pointer to the node containing the reference volume for the C3D
        FSL-to-ITK affine conversion (often the mean of the functional
        time-series, which is a single volume)
    ref_outfile: Nipype pointer
        pointer to the output of ref_node, i.e. the reference volume itself
    func_name: str
        what the name of the warped functional should be when written to the
        resource pool
    interp: str
        which interpolation to use when applying the warps
    input_image_type: int
        argument taken by the ANTs apply warp tool; in this case, should be
        3 for 4D functional time-series
    """

    # converts FSL-format .mat affine xfm into ANTS-format
    # .txt; .mat affine comes from Func->Anat registration
    fsl_to_itk_mni_func = create_wf_c3d_fsl_to_itk(
        name='fsl_to_itk_%s_%d' % (func_name, num_strat)
    )

    # collects series of warps to be applied
    collect_transforms_mni_func = \
        create_wf_collect_transforms(
            inverse=True,
            name='collect_transforms_%s_%d' % (func_name, num_strat)
        )

    # apply ants warps
    apply_ants_warp_mni_func = \
        create_wf_apply_ants_warp(
            inverse=True,
            name='apply_ants_warp_%s_%d' % (func_name, num_strat),
            ants_threads=int(num_ants_cores))

    apply_ants_warp_mni_func.inputs.inputspec.dimension = 3
    apply_ants_warp_mni_func.inputs.inputspec.interpolation = interp

    # input_image_type:
    # (0 or 1 or 2 or 3)
    # Option specifying the input image type of scalar
    # (default), vector, tensor, or time series.
    apply_ants_warp_mni_func.inputs.inputspec. \
        input_image_type = input_image_type

    # convert the .mat from linear Func->Anat to
    # ANTS format
    node, out_file = strat['functional_to_anat_linear_xfm']
    workflow.connect(node, out_file, fsl_to_itk_mni_func,
                     'inputspec.affine_file')

    node, out_file = strat["anatomical_brain"]
    workflow.connect(node, out_file, fsl_to_itk_mni_func,
                     'inputspec.reference_file')

    workflow.connect(ref_node, ref_outfile,
                     fsl_to_itk_mni_func,
                     'inputspec.source_file')

    workflow.connect(ref_node, ref_outfile,
                     apply_ants_warp_mni_func, 'inputspec.reference_image')

    # Field file from anatomical nonlinear registration
    node, out_file = strat['mni_to_anatomical_nonlinear_xfm']
    workflow.connect(node, out_file,
                     collect_transforms_mni_func,
                     'inputspec.warp_file')

    # initial transformation from anatomical registration
    node, out_file = strat['ants_initial_xfm']
    workflow.connect(node, out_file,
                     collect_transforms_mni_func,
                     'inputspec.linear_initial')

    # affine transformation from anatomical registration
    node, out_file = strat['ants_affine_xfm']
    workflow.connect(node, out_file,
                     collect_transforms_mni_func,
                     'inputspec.linear_affine')

    # rigid transformation from anatomical registration
    node, out_file = strat['ants_rigid_xfm']
    workflow.connect(node, out_file,
                     collect_transforms_mni_func,
                     'inputspec.linear_rigid')

    # Premat from Func->Anat linear reg and bbreg
    # (if bbreg is enabled)
    workflow.connect(fsl_to_itk_mni_func,
                     'outputspec.itk_transform',
                     collect_transforms_mni_func,
                     'inputspec.fsl_to_itk_affine')

    # this <node, out_file> pulls in directly because
    # it pulls in the leaf in some instances
    workflow.connect(input_node,
                     input_outfile,
                     apply_ants_warp_mni_func,
                     'inputspec.input_image')

    workflow.connect(collect_transforms_mni_func,
                     'outputspec.transformation_series',
                     apply_ants_warp_mni_func,
                     'inputspec.transforms')

    strat.update_resource_pool({
        func_name: (apply_ants_warp_mni_func, 'outputspec.output_image')
    })

    strat.append_name(apply_ants_warp_mni_func.name)

    return apply_ants_warp_mni_func
