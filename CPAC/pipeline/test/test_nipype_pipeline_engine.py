# Copyright (C) 2021-2024  C-PAC Developers

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
"""Tests for C-PAC customizations to nipype.pipeline.engine."""
import os
import pytest
from nibabel.testing import data_path
from nipype.interfaces.utility import IdentityInterface
from nipype.utils.profiler import get_system_total_memory_gb
from traits.trait_base import Undefined

from CPAC.nuisance.utils.compcor import cosine_filter
from CPAC.pipeline.nipype_pipeline_engine import (
    DEFAULT_MEM_GB, get_data_size, Node, MapNode, Workflow)
from CPAC.pipeline.nipype_pipeline_engine.plugins import (
    LegacyMultiProcPlugin,
    MultiProcPlugin,
)
from CPAC.utils.interfaces.function import Function


def get_sample_data(filepath):
    # pylint: disable=missing-function-docstring
    return filepath


def square_func(x):  # pylint: disable=invalid-name
    # pylint: disable=missing-function-docstring
    return x ** 2


square = Function(["x"], ["f_x"], square_func)


def test_MapNode():  # pylint: disable=invalid-name
    square_node = MapNode(square, name="square", iterfield=["x"])
    square_node.inputs.x = [0, 1, 2, 3]
    res = square_node.run()
    assert res.outputs.f_x == [0, 1, 4, 9]


def test_Node():  # pylint: disable=invalid-name
    square_node = Node(square, name="square")
    square_node.inputs.x = 2
    res = square_node.run()
    assert res.outputs.f_x == 4


@pytest.mark.parametrize("plugin", [LegacyMultiProcPlugin, MultiProcPlugin])
@pytest.mark.parametrize("throttle", [True, False])
def test_throttle(plugin: type, throttle: bool) -> None:
    """Test Node throttling."""
    small_estimate = 0.0003
    cosfilter_node = Node(
        Function(
            input_names=["input_image_path", "timestep"],
            output_names=["cosfiltered_img"],
            function=cosine_filter,
            as_module=True,
        ),
        name="aCompCor_cosine_filter",
        mem_gb=small_estimate,
        throttle=throttle,
    )
    wf = Workflow(name="aCompCor_cosine_filter_wf")
    wf.add_nodes([cosfilter_node])
    plugin(plugin_args={"raise_insufficient": False})._prerun_check(wf._graph)
    if throttle:
        assert cosfilter_node.mem_gb > small_estimate
        assert cosfilter_node.mem_gb <= get_system_total_memory_gb()
    else:
        assert cosfilter_node.mem_gb == small_estimate


def test_Workflow(tmpdir):  # pylint: disable=invalid-name
    example_filepath = os.path.join(data_path, 'example4d.nii.gz')
    pass_in_filepath = IdentityInterface(fields=['x'])
    pass_in_filepath.inputs.x = example_filepath
    node_1 = Node(pass_in_filepath, 'pass_in_filepath')

    # This node just returns the filepath given. We're testing the
    # just-in-time memory allocation
    node_2 = Node(Function(['filepath'], ['filepath'], get_sample_data),
                  name='get_sample_data',
                  inputs=['filepath'],
                  outputs=['filepath'])

    assert node_2.mem_gb == DEFAULT_MEM_GB

    node_2 = Node(Function(['filepath'], ['filepath'], get_sample_data),
                  name='get_sample_data',
                  inputs=['filepath'],
                  outputs=['filepath'],
                  mem_x=(0.1, 'filepath'))
    with pytest.raises(FileNotFoundError):
        assert node_2.mem_gb

    wf = Workflow('example_workflow', base_dir=tmpdir)
    wf.connect(node_1, 'x', node_2, 'filepath')

    assert wf.get_node('get_sample_data').inputs.filepath is Undefined

    out = wf.run()
    out_node = list(out.nodes)[0]
    assert out_node.mem_gb == DEFAULT_MEM_GB + get_data_size(
        example_filepath, 'xyzt') * 0.1
