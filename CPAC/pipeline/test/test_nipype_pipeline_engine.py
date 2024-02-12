import os
import pytest
from nibabel.testing import data_path
from nipype import Function
from nipype.interfaces.utility import IdentityInterface
from traits.trait_base import Undefined
from CPAC.pipeline.nipype_pipeline_engine import (
    DEFAULT_MEM_GB, get_data_size, Node, MapNode, Workflow)


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
