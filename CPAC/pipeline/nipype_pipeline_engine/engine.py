'''Module to import Nipype Pipeline engine and override some Classes.
See https://fcp-indi.github.com/docs/developer/nodes
for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.pipeline.engine.html
for Nipype's documentation.'''  # noqa E501
from nipype.pipeline import engine as pe
from functools import partialmethod

# set global default mem_gb
DEFAULT_MEM_GB = 2.0


class Node(pe.Node):
    __doc__ = pe.Node.__doc__

    __init__ = partialmethod(pe.Node.__init__, mem_gb=DEFAULT_MEM_GB)


class MapNode(pe.MapNode):
    __doc__ = f'kwarg default: mem_gb={DEFAULT_MEM_GB}\n\n{pe.MapNode.__doc__}'

    __init__ = partialmethod(pe.MapNode.__init__, mem_gb=DEFAULT_MEM_GB)
