import nipype.interfaces.fsl as fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.afni import preprocess
from CPAC.registration import create_fsl_fnirt_nonlinear_reg, \
    create_register_func_to_anat, \
    create_bbregister_func_to_anat, \
    create_wf_calculate_ants_warp

from CPAC.utils import Configuration, function, find_files
from CPAC.utils.utils import (
    get_scan_params,
    extract_output_mean,
    get_zscore,
    get_fisher_zscore
)


def z_score_standardize(wf_name, input_image_type='func_derivative',
                        opt=None):

    wf = pe.Workflow(name=wf_name)

    map_node = False
    if input_image_type == 'func_derivative_multi':
        map_node = True

    inputnode = pe.Node(util.IdentityInterface(fields=['in_file',
                                                       'mask']),
                        name='inputspec')

    z_score_std = get_zscore(map_node, 'z_score_std')

    wf.connect(inputnode, 'in_file', z_score_std, 'inputspec.input_file')
    wf.connect(inputnode, 'mask', z_score_std, 'inputspec.mask_file')

    outputnode = pe.Node(util.IdentityInterface(fields=['out_file']),
                         name='outputspec')

    wf.connect(z_score_std, 'outputspec.z_score_img', outputnode, 'out_file')

    return wf


def fisher_z_score_standardize(wf_name, label,
                               input_image_type='func_derivative', opt=None):

    wf = pe.Workflow(name=wf_name)

    map_node = False
    if input_image_type == 'func_derivative_multi':
        map_node = True

    inputnode = pe.Node(util.IdentityInterface(fields=['correlation_file',
                                                       'timeseries_oned']),
                        name='inputspec')

    fisher_z_score_std = get_fisher_zscore(label, map_node,
                                           'fisher_z_score_std')
    wf.connect(inputnode, 'correlation_file',
               fisher_z_score_std, 'inputspec.correlation_file')

    wf.connect(inputnode, 'timeseries_oned',
               fisher_z_score_std, 'inputspec.timeseries_one_d')

    outputnode = pe.Node(util.IdentityInterface(fields=['out_file']),
                         name='outputspec')

    wf.connect(fisher_z_score_std, 'outputspec.fisher_z_score_img',
               outputnode, 'out_file')

    return wf


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
