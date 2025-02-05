# Copyright (C) 2025  C-PAC Developers

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
# pylint: disable=too-many-lines,ungrouped-imports,wrong-import-order
"""Longitudial registration workflows and utilities."""

from nipype.pipeline.engine import Node, Workflow
from CPAC.pipeline.engine import ResourcePool
from CPAC.pipeline.nodeblock import NODEBLOCK_RETURN


def t1w_to_longitudinal_to_template(
    wf: Workflow, strat_pool: ResourcePool
) -> NODEBLOCK_RETURN:
    """Combine and apply T1w transforms from native to longitudinal to template."""
    outputs: dict[str, tuple[Node | Workflow, str]] = {}
    return wf, outputs
