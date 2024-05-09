# Copyright (C) 2023-2024  C-PAC Developers

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
"""Custom serialization for flowdump."""

from typing import Callable

from traits.has_traits import HasTraits
from nipype.interfaces.base.support import Bunch

from . import Configuration


def _truncate_large_cpac_internals(obj: object) -> object:
    """Recursively replace json_data values and Config objects."""
    if isinstance(obj, Configuration):
        return "[C-PAC config]"
    if isinstance(obj, list):
        return [_truncate_large_cpac_internals(i) for i in obj]
    if isinstance(obj, tuple):
        return (_truncate_large_cpac_internals(i) for i in obj)
    if isinstance(obj, HasTraits):
        return _truncate_large_cpac_internals(obj.trait_get())
    if isinstance(obj, (Bunch, dict)):
        return {
            str(k): _truncate_large_cpac_internals(v)
            if k != "json_data"
            else "[Truncated]"
            for k, v in obj.items()
        }
    return obj


def cpac_flowdump_serializer(
    flowdump_serializer: Callable[[object], object], obj: object
) -> object:
    """Remove `json_data` fields.

    Custom flowdump serializer that removes `json_data` fields
    and CPAC configs from workflow json files as these are repeated
    for every node (and increase file size dramatically).
    """
    return flowdump_serializer(_truncate_large_cpac_internals(obj))
