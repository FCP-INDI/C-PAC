import os

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.registration import fsl_apply_transform_func_to_mni
from CPAC.utils import Configuration, Strategy, Outputs

from CPAC.pipeline.cpac_group_runner import gather_outputs
from CPAC.utils.interfaces.function import Function
from CPAC.utils.datasource import resolve_resolution

def file_node(path, file_node_num=0):
    input_node = pe.Node(
        util.IdentityInterface(fields=['file']), name='file_node_{0}'.format(file_node_num)
    )
    input_node.inputs.file = path
    return input_node, 'file'

#templates_for_resampling = [
    #(c.resolution_for_func_preproc, c.template_brain_only_for_func, 'template_brain_for_func_preproc', 'resolution_for_func_preproc'),
    #(c.resolution_for_func_preproc, c.template_skull_for_func, 'template_skull_for_func_preproc', 'resolution_for_func_preproc'),
    #(c.resolution_for_func_derivative, c.template_brain_only_for_func, 'template_brain_for_func_derivative', 'resolution_for_func_preproc'),
    #(c.resolution_for_func_derivative, c.template_skull_for_func, 'template_skull_for_func_derivative', 'resolution_for_func_preproc')
#]

def test_fsl_apply_transform_func_to_mni_nonlinear():


    # resources needed by workflow
    #
    #  template_brain_for_func_preproc
    #  anat_mni_fnirt_register
    #  anatomical_to_mni_nonlinear_xfm
    #  functional_to_anat_linear_xfm
    #
    #  anat_mni_flirt_register
    #  functional_to_mni_linear_xfm [?]
    #  anatomical_to_mni_linear_xfm
    #  functional_to_anat_linear_xfm
    # 

    outputs = Outputs()

    # mock the config dictionary
    c = Configuration({
        "workingDirectory": "/scratch/pipeline_tests",
        "crashLogDirectory": "/scratch",
        "outputDirectory": "/output/output/pipeline_analysis_nuisance/sub-M10978008_ses-NFB3",
        "resolution_for_func_preproc": "3mm",
        "resolution_for_func_derivative": "3mm",
        "template_for_resample": "/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz",
        "template_brain_only_for_func": "/usr/share/fsl/5.0/data/standard/MNI152_T1_${resolution_for_func_preproc}_brain.nii.gz",
        "template_skull_for_func":  "/usr/share/fsl/5.0/data/standard/MNI152_T1_${resolution_for_func_preproc}.nii.gz",
        "identityMatrix":  "/usr/share/fsl/5.0/etc/flirtsch/ident.mat",
        "funcRegFSLinterpolation": "sinc"

    })

    # mock the strategy
    strat = Strategy()

    strat.append_name('anat_mni_fnirt_register_0')

    resource_dict = {
            "mean_functional": os.path.join(c.outputDirectory,
                "mean_functional/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat.nii.gz"),
            "motion_correct": os.path.join(c.outputDirectory,
                "motion_correct/_scan_test/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg.nii.gz"),
            "anatomical_to_mni_nonlinear_xfm": os.path.join(c.outputDirectory,
                "anatomical_to_mni_nonlinear_xfm/sub-M10978008_ses-NFB3_T1w_resample_fieldwarp.nii.gz"),
            "functional_to_anat_linear_xfm": os.path.join(c.outputDirectory,
                "functional_to_anat_linear_xfm/_scan_test/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_flirt.mat")
    }
   
    file_node_num = 0
    for resource, filepath in resource_dict.items():
        print('resource: {0}, filename: {1}'.format(resource, filepath))
        strat.update_resource_pool({
            resource: file_node(filepath, file_node_num)
        })
        file_node_num += 1

    templates_for_resampling = [
        (c.resolution_for_func_preproc, c.template_brain_only_for_func,
            'template_brain_for_func_preproc', 'resolution_for_func_preproc')
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

    # build the workflow
    workflow = pe.Workflow(name='test_fsl_apply_transform_func_to_mni_nonlinear')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = fsl_apply_transform_func_to_mni(workflow, 'mean_functional_to_standard',
            'mean_functional', 'template_brain_for_func_preproc', 0, strat, c.funcRegFSLinterpolation)

    workflow.run()

def test_fsl_apply_transform_func_to_mni_linear():


    # resources needed by workflow
    #
    #  anat_mni_flirt_register
    #  functional_to_mni_linear_xfm [?]
    #  anatomical_to_mni_linear_xfm
    #  functional_to_anat_linear_xfm
    # 

    outputs = Outputs()

    # mock the config dictionary
    c = Configuration({
        "workingDirectory": "/scratch/pipeline_tests",
        "crashLogDirectory": "/scratch",
        "outputDirectory": "/output/output/pipeline_analysis_nuisance/sub-M10978008_ses-NFB3",
        "resolution_for_func_preproc": "3mm",
        "resolution_for_func_derivative": "3mm",
        "template_for_resample": "/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz",
        "template_brain_only_for_func": "/usr/share/fsl/5.0/data/standard/MNI152_T1_${resolution_for_func_preproc}_brain.nii.gz",
        "template_skull_for_func":  "/usr/share/fsl/5.0/data/standard/MNI152_T1_${resolution_for_func_preproc}.nii.gz",
        "identityMatrix":  "/usr/share/fsl/5.0/etc/flirtsch/ident.mat",
        "funcRegFSLinterpolation": "sinc"

    })

    # mock the strategy
    strat = Strategy()

    strat.append_name('anat_mni_flirt_register_0')

    resource_dict = {
            "mean_functional": os.path.join(c.outputDirectory,
                "mean_functional/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat.nii.gz"),
            "motion_correct": os.path.join(c.outputDirectory,
                "motion_correct/_scan_test/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg.nii.gz"),
            "anatomical_to_mni_linear_xfm": os.path.join(c.outputDirectory,
                "anatomical_to_mni_linear_xfm/sub-M10978008_ses-NFB3_T1w_resample_calc_flirt.mat"),
            "functional_to_anat_linear_xfm": os.path.join(c.outputDirectory,
                "functional_to_anat_linear_xfm/_scan_test/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_flirt.mat")
    }
   
    file_node_num = 0
    for resource, filepath in resource_dict.items():
        print('resource: {0}, filename: {1}'.format(resource, filepath))
        strat.update_resource_pool({
            resource: file_node(filepath, file_node_num)
        })
        file_node_num += 1

    templates_for_resampling = [
        (c.resolution_for_func_preproc, c.template_brain_only_for_func,
            'template_brain_for_func_preproc', 'resolution_for_func_preproc')
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

    # build the workflow
    workflow = pe.Workflow(name='test_fsl_apply_transform_func_to_mni_linear')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = fsl_apply_transform_func_to_mni(workflow, 'mean_functional_to_standard',
            'mean_functional', 'template_brain_for_func_preproc', 0, strat, c.funcRegFSLinterpolation)

    workflow.run()
