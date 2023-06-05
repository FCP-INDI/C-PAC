# Copyright (C) 2022  C-PAC Developers

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
'''Module to import Nipype Pipeline engine and override some Classes.
See https://fcp-indi.github.io/docs/developer/nodes
for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.pipeline.engine.html
for Nipype's documentation.'''  # noqa: E501  # pylint: disable=line-too-long
from nipype.pipeline import engine as pe
# import everything in nipype.pipeline.engine.__all__
from nipype.pipeline.engine import *  # noqa: F401,F403
# import our DEFAULT_MEM_GB and override Node, MapNode
from .engine import DEFAULT_MEM_GB, export_graph, get_data_size, Node, \
                    MapNode, UNDEFINED_SIZE, Workflow

__all__ = [
    interface for interface in dir(pe) if not interface.startswith('_')
] + ['DEFAULT_MEM_GB', 'export_graph', 'get_data_size', 'Node', 'MapNode',
     'UNDEFINED_SIZE', 'Workflow']

del pe
