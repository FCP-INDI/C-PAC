from nipype import logging
logger = logging.getLogger('workflow')

import nipype.pipeline.engine as pe

import nipype.interfaces.afni as afni

from CPAC.utils.interfaces.function import Function
from CPAC.utils.utils import (
    get_scan_params
)

from CPAC.utils.datasource import (
    create_func_datasource,
    create_fmap_datasource,
    get_fmap_phasediff_metadata,
    calc_deltaTE_and_asym_ratio
)


def connect_func_ingress(workflow, strat_list, c, sub_dict, subject_id,
                         input_creds_path, unique_id=None):

    for num_strat, strat in enumerate(strat_list):

        if 'func' in sub_dict:
            func_paths_dict = sub_dict['func']
        else:
            func_paths_dict = sub_dict['rest']

        if unique_id is None:
            workflow_name=f'func_gather_{num_strat}'
        else:
            workflow_name=f'func_gather_{unique_id}_{num_strat}'

        func_wf = create_func_datasource(func_paths_dict,
                                        workflow_name)

        func_wf.inputs.inputnode.set(
            subject=subject_id,
            creds_path=input_creds_path,
            dl_dir=c.pipeline_setup['working_directory']['path']
        )
        func_wf.get_node('inputnode').iterables = \
            ("scan", list(func_paths_dict.keys()))

        strat.update_resource_pool({
            'subject': (func_wf, 'outputspec.subject'),
            'scan': (func_wf, 'outputspec.scan')
        })

        # Grab field maps
        diff = False
        blip = False
        fmap_rp_list = []
        fmap_TE_list = []
        if "fmap" in sub_dict:
            for key in sub_dict["fmap"]:
                gather_fmap = create_fmap_datasource(sub_dict["fmap"],
                                                     "fmap_gather_"
                                                     "{0}".format(key))
                gather_fmap.inputs.inputnode.set(
                    subject=subject_id,
                    creds_path=input_creds_path,
                    dl_dir=c.pipeline_setup['working_directory']['path']
                )
                gather_fmap.inputs.inputnode.scan = key
                strat.update_resource_pool({
                    key: (gather_fmap, 'outputspec.rest'),
                    "{0}_scan_params".format(key): (gather_fmap,
                                                    'outputspec.scan_params')
                })

                fmap_rp_list.append(key)

                if key == "diff_phase" or key == "diff_mag_one" or \
                                key == "diff_mag_two":
                    diff = True

                    get_fmap_metadata_imports = ['import json']
                    get_fmap_metadata = pe.Node(Function(
                        input_names=['data_config_scan_params'],
                        output_names=['echo_time',
                                      'dwell_time',
                                      'pe_direction'],
                        function=get_fmap_phasediff_metadata,
                        imports=get_fmap_metadata_imports),
                        name='{0}_get_metadata_{1}'.format(key,
                                                           num_strat))

                    node, out_file = strat["{}_scan_params".format(key)]
                    workflow.connect(node, out_file, get_fmap_metadata,
                                     'data_config_scan_params')

                    strat.update_resource_pool({
                        "{}_TE".format(key): (get_fmap_metadata,
                                              'echo_time'),
                        "{}_dwell".format(key): (get_fmap_metadata,
                                                 'dwell_time'),
                        "{}_pedir".format(key): (get_fmap_metadata,
                                                 'pe_direction')
                    })
                    fmap_TE_list.append("{}_TE".format(key))

                if key == "epi_AP" or key == "epi_PA":
                    blip = True

            if diff:
                calc_delta_ratio = pe.Node(Function(
                    input_names=['dwell_time',
                                 'echo_time_one',
                                 'echo_time_two',
                                 'echo_time_three'],
                    output_names=['deltaTE',
                                  'dwell_asym_ratio'],
                    function=calc_deltaTE_and_asym_ratio),
                    name='diff_distcor_calc_delta_{}'.format(num_strat))

                node, out_file = strat['diff_phase_dwell']
                workflow.connect(node, out_file, calc_delta_ratio,
                                 'dwell_time')

                node, out_file = strat[fmap_TE_list[0]]
                workflow.connect(node, out_file, calc_delta_ratio,
                                 'echo_time_one')

                node, out_file = strat[fmap_TE_list[1]]
                workflow.connect(node, out_file, calc_delta_ratio,
                                 'echo_time_two')

                if len(fmap_TE_list) > 2:
                    node, out_file = strat[fmap_TE_list[2]]
                    workflow.connect(node, out_file, calc_delta_ratio,
                                     'echo_time_three')

                strat.update_resource_pool({
                    'deltaTE': (calc_delta_ratio, 'deltaTE'),
                    'dwell_asym_ratio': (calc_delta_ratio,
                                         'dwell_asym_ratio')
                })

        # Add in nodes to get parameters from configuration file
        # a node which checks if scan_parameters are present for each scan
        if unique_id is None:
            workflow_name=f'scan_params_{num_strat}'
        else:
            workflow_name=f'scan_params_{unique_id}_{num_strat}'

        scan_params = \
            pe.Node(Function(
                input_names=['data_config_scan_params',
                             'subject_id',
                             'scan',
                             'pipeconfig_tr',
                             'pipeconfig_tpattern',
                             'pipeconfig_start_indx',
                             'pipeconfig_stop_indx'],
                output_names=['tr',
                              'tpattern',
                              'ref_slice',
                              'start_indx',
                              'stop_indx',
                              'pe_direction'],
                function=get_scan_params,
                as_module=True
            ), name=workflow_name)

        if "Selected Functional Volume" in c.functional_registration['1-coregistration']['func_input_prep']['input']:
            get_func_volume = pe.Node(interface=afni.Calc(),
                                      name='get_func_volume_{0}'.format(
                                          num_strat))

            get_func_volume.inputs.set(
                expr='a',
                single_idx=c.functional_registration['1-coregistration']['func_input_prep']['Selected Functional Volume']['func_reg_input_volume'],
                outputtype='NIFTI_GZ'
            )
            workflow.connect(func_wf, 'outputspec.rest',
                             get_func_volume, 'in_file_a')

        # wire in the scan parameter workflow
        workflow.connect(func_wf, 'outputspec.scan_params',
                         scan_params, 'data_config_scan_params')

        workflow.connect(func_wf, 'outputspec.subject',
                         scan_params, 'subject_id')

        workflow.connect(func_wf, 'outputspec.scan',
                         scan_params, 'scan')

        # connect in constants
        scan_params.inputs.set(
            pipeconfig_start_indx=c.functional_preproc['truncation']['start_tr'],
            pipeconfig_stop_indx=c.functional_preproc['truncation']['stop_tr']
        )

        strat.update_resource_pool({
            'raw_functional': (func_wf, 'outputspec.rest'),
            'scan_id': (func_wf, 'outputspec.scan'),
            'tr': (scan_params, 'tr'),
            'tpattern': (scan_params, 'tpattern'),
            'start_idx': (scan_params, 'start_indx'),
            'stop_idx': (scan_params, 'stop_indx'),
            'pe_direction': (scan_params, 'pe_direction'),
        })

        strat.set_leaf_properties(func_wf, 'outputspec.rest')

        if "Selected Functional Volume" in c.functional_registration['1-coregistration']['func_input_prep']['input']:
            strat.update_resource_pool({
                'selected_func_volume': (get_func_volume, 'out_file')
            })

    return (workflow, diff, blip, fmap_rp_list)
