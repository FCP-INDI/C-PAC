# Copyright (C) 2020-2022  C-PAC Developers

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
from nipype import logging
from CPAC.utils.datasource import (
    create_func_datasource,
    ingress_func_metadata)
logger = logging.getLogger('nipype.workflow')


def connect_func_ingress(workflow, strat_list, c, sub_dict, subject_id,
                         input_creds_path, unique_id=None):
    for num_strat, strat in enumerate(strat_list):

        if 'func' in sub_dict:
            func_paths_dict = sub_dict['func']
        else:
            func_paths_dict = sub_dict['rest']

        if unique_id is None:
            workflow_name = f'func_gather_{num_strat}'
        else:
            workflow_name = f'func_gather_{unique_id}_{num_strat}'

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

        (workflow, strat.rpool, diff, blip, fmap_rp_list
         ) = ingress_func_metadata(workflow, c, strat.rpool, sub_dict,
                                   subject_id, input_creds_path, unique_id)

    return (workflow, diff, blip, fmap_rp_list)
