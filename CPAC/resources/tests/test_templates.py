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
"""Tests for packaged templates"""
import os
import pytest
from CPAC.pipeline import ALL_PIPELINE_CONFIGS
from CPAC.pipeline.engine import ingress_pipeconfig_paths, ResourcePool
from CPAC.utils.configuration import Preconfiguration
from CPAC.utils.datasource import get_highest_local_res


@pytest.mark.parametrize('pipeline', ALL_PIPELINE_CONFIGS)
def test_packaged_path_exists(pipeline):
    """
    Check that all local templates are included in image at at
    least one resolution
    """
    rpool = ingress_pipeconfig_paths(Preconfiguration(pipeline),
                                     ResourcePool(), 'pytest')
    for resource in rpool.rpool.values():
        node = list(resource.values())[0].get('data')[0]
        if hasattr(node.inputs, 'template'
                   ) and not node.inputs.template.startswith('s3:'):
            if not pipeline == 'rodent' and node.inputs.template.startswith(
                '/template/study_based'):
                assert (os.path.exists(node.inputs.template) or
                        get_highest_local_res(node.inputs.template,
                                              node.inputs.resolution).exists())
