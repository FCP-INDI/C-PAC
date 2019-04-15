import nipype.interfaces.fsl as fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.afni import preprocess
from CPAC.utils import Configuration, function, find_files
from CPAC.nuisance import create_nuisance_workflow, bandpass_voxels, NuisanceRegressor
import copy


def connect_nuisance(workflow,strat, c, regressors_selector, has_segmentation, use_ants,num_strat):

    new_strat = strat.fork()

    # to guarantee immutability
    regressors_selector = NuisanceRegressor(
        copy.deepcopy(regressors_selector),
        copy.deepcopy(c.Regressors)
    )

    # remove tissue regressors when there is no segmentation
    # on the strategy
    if not has_segmentation:
        for reg in ['aCompCor',
                    'WhiteMatter',
                    'GreyMatter',
                    'CerebrospinalFluid']:

            if reg in regressors_selector:
                del regressors_selector[reg]

    nuisance_regression_workflow = create_nuisance_workflow(
        regressors_selector,
        use_ants=use_ants,
        name='nuisance_{0}_{1}'.format(regressors_selector, num_strat)
    )

    node, out_file = new_strat['anatomical_brain']
    workflow.connect(
        node, out_file,
        nuisance_regression_workflow, 'inputspec.anatomical_file_path'
    )

    if has_segmentation:
        workflow.connect(
            c.lateral_ventricles_mask, 'local_path',
            nuisance_regression_workflow, 'inputspec.lat_ventricles_mask_file_path'
        )

        node, out_file = new_strat['anatomical_gm_mask']
        workflow.connect(
            node, out_file,
            nuisance_regression_workflow, 'inputspec.gm_mask_file_path'
        )

        node, out_file = new_strat['anatomical_wm_mask']
        workflow.connect(
            node, out_file,
            nuisance_regression_workflow, 'inputspec.wm_mask_file_path'
        )

        node, out_file = new_strat['anatomical_csf_mask']
        workflow.connect(
            node, out_file,
            nuisance_regression_workflow, 'inputspec.csf_mask_file_path'
        )

    node, out_file = new_strat['movement_parameters']
    workflow.connect(
        node, out_file,
        nuisance_regression_workflow,
        'inputspec.motion_parameters_file_path'
    )

    node, out_file = new_strat['functional_to_anat_linear_xfm']
    workflow.connect(
        node, out_file,
        nuisance_regression_workflow,
        'inputspec.func_to_anat_linear_xfm_file_path'
    )

    node, out_file = new_strat.get_leaf_properties()
    workflow.connect(
        node, out_file,
        nuisance_regression_workflow,
        'inputspec.functional_file_path'
    )

    node, out_file = new_strat['frame_wise_displacement_jenkinson']
    workflow.connect(
        node, out_file,
        nuisance_regression_workflow,
        'inputspec.fd_j_file_path'
    )

    node, out_file = new_strat['frame_wise_displacement_power']
    workflow.connect(
        node, out_file,
        nuisance_regression_workflow,
        'inputspec.fd_p_file_path'
    )

    node, out_file = new_strat['dvars']
    workflow.connect(
        node, out_file,
        nuisance_regression_workflow,
        'inputspec.dvars_file_path'
    )

    node, out_file = new_strat['functional_brain_mask']
    workflow.connect(
        node, out_file,
        nuisance_regression_workflow,
        'inputspec.functional_brain_mask_file_path'
    )

    nuisance_regression_workflow.get_node('inputspec').iterables = ([
        ('selector', [regressors_selector]),
    ])

    if use_ants:

        # pass the ants_affine_xfm to the input for the
        # INVERSE transform, but ants_affine_xfm gets inverted
        # within the workflow

        node, out_file = new_strat['ants_initial_xfm']
        workflow.connect(
            node, out_file,
            nuisance_regression_workflow,
            'inputspec.anat_to_mni_initial_xfm_file_path'
        )

        node, out_file = new_strat['ants_rigid_xfm']
        workflow.connect(
            node, out_file,
            nuisance_regression_workflow,
            'inputspec.anat_to_mni_rigid_xfm_file_path'
        )

        node, out_file = new_strat['ants_affine_xfm']
        workflow.connect(
            node, out_file,
            nuisance_regression_workflow,
            'inputspec.anat_to_mni_affine_xfm_file_path'
        )
    else:
        node, out_file = new_strat['mni_to_anatomical_linear_xfm']
        workflow.connect(
            node, out_file,
            nuisance_regression_workflow,
            'inputspec.mni_to_anat_linear_xfm_file_path'
        )

    new_strat.append_name(nuisance_regression_workflow.name)

    new_strat.set_leaf_properties(
        nuisance_regression_workflow,
        'outputspec.residual_file_path'
    )

    new_strat.update_resource_pool({
        'nuisance_regression_selector': regressors_selector,

        'functional_nuisance_residuals': (
            nuisance_regression_workflow,
            'outputspec.residual_file_path')
        ,
        'functional_nuisance_regressors': (
            nuisance_regression_workflow,
            'outputspec.regressors_file_path'
        ),
    })

    return new_strat
