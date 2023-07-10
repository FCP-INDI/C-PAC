
import os

import nibabel as nb
import numpy as np

from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
from CPAC.utils.interfaces.fsl import Merge as fslMerge

from CPAC.utils.interfaces.function import Function

def length(it):
    return len(it)

def detect_qpp(datasets,
               joint_datasets, joint_mask,
               window_length, permutations,
               lower_correlation_threshold, higher_correlation_threshold,
               correlation_threshold_iteration,
               iterations, convergence_iterations):
    
    from CPAC.qpp.qpp import detect_qpp

    joint_mask_img = nb.load(joint_mask)
    joint_mask = joint_mask_img.get_fdata().astype(bool)
    joint_datasets_img = nb.load(joint_datasets)
    joint_datasets = joint_datasets_img.get_fdata()[joint_mask]

    correlation_threshold = lambda i: \
        higher_correlation_threshold \
        if i > correlation_threshold_iteration else \
        lower_correlation_threshold

    best_template_segment, _, _ = detect_qpp(
        joint_datasets,
        datasets,
        window_length,
        permutations,
        correlation_threshold,
        iterations,
        convergence_iterations
    )

    qpp = np.zeros(joint_datasets_img.shape[0:3] + (window_length,))
    qpp[joint_mask] = best_template_segment

    qpp_img = nb.Nifti1Image(qpp, joint_mask_img.affine)
    qpp_img.to_filename('./qpp.nii.gz')

    return os.path.abspath('./qpp.nii.gz')


def create_qpp(name='qpp', working_dir=None, crash_dir=None):
    
    if not working_dir:
        working_dir = os.path.join(os.getcwd(), 'QPP_work_dir')
    if not crash_dir:
        crash_dir = os.path.join(os.getcwd(), 'QPP_crash_dir')

    workflow = pe.Workflow(name=name)
    workflow.base_dir = working_dir
    workflow.config['execution'] = {'hash_method': 'timestamp',
                                    'crashdump_dir': os.path.abspath(crash_dir)}

    inputspec = pe.Node(util.IdentityInterface(fields=[
        'datasets',
        'window_length',
        'permutations',
        'lower_correlation_threshold',
        'higher_correlation_threshold',
        'correlation_threshold_iteration',
        'iterations',
        'convergence_iterations',
    ]), name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['qpp']),
                         name='outputspec')

    merge = pe.Node(fslMerge(), name='joint_datasets')
    merge.inputs.dimension = 't'
    merge.inputs.output_type = 'NIFTI_GZ'

    mask = pe.Node(interface=fsl.ImageMaths(), name='joint_mask')
    mask.inputs.op_string = '-abs -Tmin -bin'

    detect = pe.Node(Function(input_names=['datasets',
                                           'joint_datasets',
                                           'joint_mask',
                                           'window_length',
                                           'permutations',
                                           'lower_correlation_threshold',
                                           'higher_correlation_threshold',
                                           'correlation_threshold_iteration',
                                           'iterations',
                                           'convergence_iterations'],
                                output_names=['qpp'],
                                function=detect_qpp,
                                as_module=True),
                     name='detect_qpp')
    
    workflow.connect([
        (inputspec, merge, [('datasets', 'in_files')]),
        (merge, mask, [('merged_file', 'in_file')]),
        (merge, detect, [('merged_file', 'joint_datasets')]),
        (mask, detect, [('out_file', 'joint_mask')]),
        (inputspec, detect, [
            (('datasets', length), 'datasets'),
            ('window_length' ,'window_length'),
            ('permutations' ,'permutations'),
            ('lower_correlation_threshold' ,'lower_correlation_threshold'),
            ('higher_correlation_threshold' ,'higher_correlation_threshold'),
            ('correlation_threshold_iteration' ,'correlation_threshold_iteration'),
            ('iterations' ,'iterations'),
            ('convergence_iterations' ,'convergence_iterations'),
        ]),
        (detect, outputspec, [('qpp', 'qpp')]),
    ])

    return workflow
