import os
from nipype.interfaces import utility as util
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.configuration import Configuration
from CPAC.utils.datasource import resolve_resolution
from CPAC.utils.interfaces.function import Function
from CPAC.utils.strategy import Strategy

def file_node(path, file_node_num=0):
    input_node = pe.Node(
        util.IdentityInterface(fields=['file']), name='file_node_{0}'.format(file_node_num)
    )
    input_node.inputs.file = path
    return input_node, 'file'

def configuration_strategy_mock( method = 'FSL' ):
    fsldir = os.environ.get("FSLDIR")
    # mock the config dictionary
    c = Configuration({
        "num_ants_threads": 4,
        "workingDirectory": "/scratch/pipeline_tests",
        "crashLogDirectory": "/scratch",
        "outputDirectory": "/output/output/pipeline_analysis_nuisance/sub-M10978008_ses-NFB3",
        "resolution_for_func_preproc": "3mm",
        "resolution_for_func_derivative": "3mm",
        "template_for_resample": f"{fsldir}/data/standard/"
                                 "MNI152_T1_1mm_brain.nii.gz",
        "template_brain_only_for_func": f"{fsldir}/data/standard/"
                                        r"MNI152_T1_${func_resolution}_"
                                        "brain.nii.gz",
        "template_skull_for_func":  f"{fsldir}/data/standard/"
                                    r"MNI152_T1_${func_resolution}.nii.gz",
        "identityMatrix":  f"{fsldir}/etc/flirtsch/ident.mat",
        "funcRegFSLinterpolation": "sinc",
        "funcRegANTSinterpolation": "LanczosWindowedSinc"
    })

    if method == 'ANTS':
        c.update('regOption', 'ANTS')
    else:
        c.update('regOption', 'FSL')

    # mock the strategy
    strat = Strategy()

    resource_dict = {
            "mean_functional": os.path.join(c.outputDirectory,
                "mean_functional/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat.nii.gz"),
            "motion_correct": os.path.join(c.outputDirectory,
                "motion_correct/_scan_test/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg.nii.gz"),
            "anatomical_brain": os.path.join(c.outputDirectory,
                "anatomical_brain/sub-M10978008_ses-NFB3_acq-ao_brain_resample.nii.gz"),
            "ants_initial_xfm": os.path.join(c.outputDirectory,
                "ants_initial_xfm/transform0DerivedInitialMovingTranslation.mat"),
            "ants_affine_xfm": os.path.join(c.outputDirectory,
                "ants_affine_xfm/transform2Affine.mat"),
            "ants_rigid_xfm": os.path.join(c.outputDirectory,
                "ants_rigid_xfm/transform1Rigid.mat"),
            "anatomical_to_mni_linear_xfm": os.path.join(c.outputDirectory,
                "anatomical_to_mni_linear_xfm/sub-M10978008_ses-NFB3_T1w_resample_calc_flirt.mat"),
            "functional_to_anat_linear_xfm": os.path.join(c.outputDirectory,
                "functional_to_anat_linear_xfm/_scan_test/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_flirt.mat"),
            'ants_symm_warp_field': os.path.join(c.outputDirectory,
                "anatomical_to_symmetric_mni_nonlinear_xfm/transform3Warp.nii.gz"),
            'ants_symm_affine_xfm': os.path.join(c.outputDirectory,
                "ants_symmetric_affine_xfm/transform2Affine.mat"),
            'ants_symm_rigid_xfm': os.path.join(c.outputDirectory,
                "ants_symmetric_rigid_xfm/transform1Rigid.mat"),
            'ants_symm_initial_xfm': os.path.join(c.outputDirectory,
                "ants_symmetric_initial_xfm/transform0DerivedInitialMovingTranslation.mat"),
            "dr_tempreg_maps_files": [os.path.join('/scratch', 'resting_preproc_sub-M10978008_ses-NFB3_cpac105', 'temporal_dual_regression_0/_scan_test/_selector_CSF-2mmE-M_aC-WM-2mmE-DPC5_G-M_M-SDB_P-2/_spatial_map_PNAS_Smith09_rsn10_spatial_map_file_..cpac_templates..PNAS_Smith09_rsn10.nii.gz/split_raw_volumes/temp_reg_map_000{0}.nii.gz'.format(n)) for n in range(10)]
    }

    if method == 'ANTS':
        resource_dict["anatomical_to_mni_nonlinear_xfm"] = os.path.join(c.outputDirectory,
            "anatomical_to_mni_nonlinear_xfm/transform3Warp.nii.gz")
    else:
        resource_dict["anatomical_to_mni_nonlinear_xfm"] = os.path.join(c.outputDirectory,
            "anatomical_to_mni_nonlinear_xfm/sub-M10978008_ses-NFB3_T1w_resample_fieldwarp.nii.gz")
   
    file_node_num = 0
    for resource, filepath in resource_dict.items():
        strat.update_resource_pool({
            resource: file_node(filepath, file_node_num)
        })
        strat.append_name(resource+'_0')
        file_node_num += 1

    templates_for_resampling = [
        (c.resolution_for_func_preproc, c.template_brain_only_for_func,
            'template_brain_for_func_preproc', 'resolution_for_func_preproc'),
        (c.resolution_for_func_preproc, c.template_brain_only_for_func,
            'template_skull_for_func_preproc', 'resolution_for_func_preproc')
    ]

    for resolution, template, template_name, tag in templates_for_resampling:
        resampled_template = pe.Node(Function(input_names = ['resolution', 'template', 'template_name', 'tag'],
                                              output_names = ['resampled_template'],
                                              function = resolve_resolution,
                                              as_module = True),
                                        name = 'resampled_' + template_name)

        resampled_template.inputs.resolution = resolution
        resampled_template.inputs.template = template
        resampled_template.inputs.template_name = template_name
        resampled_template.inputs.tag = tag

        strat.update_resource_pool({template_name: (resampled_template, 'resampled_template')})
        strat.append_name('resampled_template_0')

    return c, strat
