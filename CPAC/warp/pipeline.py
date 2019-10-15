
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
    set_gauss,
    get_scan_params,
    extract_output_mean,
    get_zscore,
    get_fisher_zscore
)

# Apply  Z-scoring, Smoothing, Averages
# use derivative
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
