from nipype import Function
from CPAC.pipeline.nipype_pipeline_engine import Node, MapNode


def square_func(x):
    return x ** 2


square = Function(["x"], ["f_x"], square_func)


def test_MapNode():
    square_node = MapNode(square, name="square", iterfield=["x"])
    square_node.inputs.x = [0, 1, 2, 3]
    res = square_node.run()
    assert res.outputs.f_x == [0, 1, 4, 9]


def test_Node():
    square_node = Node(square, name="square")
    square_node.inputs.x = 2
    res = square_node.run()
    assert res.outputs.f_x == 4
