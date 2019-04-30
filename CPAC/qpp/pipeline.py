
import os
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl

from CPAC.utils.function import Function


def create_qpp(name='qpp', working_dir=None, crash_dir=None):
    
    if not working_dir:
        working_dir = os.path.join(os.getcwd(), 'QPP_work_dir')
    if not crash_dir:
        crash_dir = os.path.join(os.getcwd(), 'QPP_crash_dir')

    workflow = pe.Workflow(name=name)
    workflow.base_dir = working_dir
    workflow.config['execution'] = {'hash_method': 'timestamp',
                                    'crashdump_dir': os.path.abspath(crash_dir)}

    inputspec = pe.Node(util.IdentityInterface(fields=['datasets',
                                                       'window_length',
                                                       'permutations',
                                                       'lower_correlation_threshold',
                                                       'higher_correlation_threshold',
                                                       'max_iterations',
                                                       'convergence_iterations']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['qpp']),
                         name='outputspec')

    merge = pe.Node(fsl.Merge()), name='joint_datasets')
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
                                           'max_iterations',
                                           'convergence_iterations']),
                                output_names=['qpp'],
                                function=detect_qpp,
                                as_module=True),
                       name='detect_qpp')
    
    wf.connect([
        (inputspec, merge, [('datasets', 'in_files')]),
        (merge, mask, [('merged_file', 'joint_datasets')]),
        (merge, detect, [('out_file', 'joint_datasets')]),
        (mask, detect, [('out_file', 'joint_mask')]),
        (inputspec, detect, [
            ('datasets', 'datasets'),
            ('window_length' ,'window_length')
            ('permutations' ,'permutations')
            ('lower_correlation_threshold' ,'lower_correlation_threshold')
            ('higher_correlation_threshold' ,'higher_correlation_threshold')
            ('max_iterations' ,'max_iterations')
            ('convergence_iterations' ,'convergence_iterations')
        ]),
        (detect, outputspec, [('qpp', 'qpp')]),
    ])

    return workflow
