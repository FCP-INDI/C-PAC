# Copyright (C) 2012-2024  C-PAC Developers

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
"""General utilities for C-PAC."""

from . import build_data_config, create_fsl_flame_preset, versioning
from .configuration import check_pname, Configuration, set_subject
from .datatypes import ListFromItem
from .interfaces import function
from .sklearn import check_random_state
from .utils import (
    correlation,
    find_files,
    get_zscore,
    repickle,
    safe_shape,
)

__all__ = [
    "build_data_config",
    "check_pname",
    "check_random_state",
    "Configuration",
    "correlation",
    "create_fsl_flame_preset",
    "find_files",
    "function",
    "get_zscore",
    "ListFromItem",
    "repickle",
    "safe_shape",
    "set_subject",
    "versioning",
]
