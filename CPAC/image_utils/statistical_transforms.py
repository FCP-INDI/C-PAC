# Copyright (C) 2018-2022  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
from nipype.interfaces import utility as util
from nipype.interfaces.afni import preprocess
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils import function
from CPAC.utils.utils import (
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
