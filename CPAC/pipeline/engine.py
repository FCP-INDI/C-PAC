'''Module to import Nipype Pipeline engine and override some Classes.

See https://fcp-indi.github.com/docs/developer/nodes for
C-PAC-specific documentation.

See https://nipype.readthedocs.io/en/latest/api/generated/nipype.pipeline.engine.html
for Nipype's documentation.'''  # noqa E501
from CPAC.utils import set_new_default_parameter
from nipype.pipeline.engine import *  # noqa: F403

DEFAULT_MEM_GB = 2.0

Node = set_new_default_parameter(Node, 'mem_gb', DEFAULT_MEM_GB)  # noqa F405
MapNode = set_new_default_parameter(MapNode, 'mem_gb', DEFAULT_MEM_GB)  # noqa F405
