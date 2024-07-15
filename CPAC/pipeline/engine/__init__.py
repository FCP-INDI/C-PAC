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
"""C-PAC engine."""

from .engine import (
    NodeBlock,
    PIPELINE_BLOCKS,
    run_node_blocks,
    wrap_block,
)
from .resource import ResourcePool, StratPool

__all__ = [
    "NodeBlock",
    "PIPELINE_BLOCKS",
    "ResourcePool",
    "StratPool",
    "run_node_blocks",
    "wrap_block",
]
