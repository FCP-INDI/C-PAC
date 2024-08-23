# Copyright (C) 2024  C-PAC Developers

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
"""Type aliases for C-PAC."""

from typing import ForwardRef

from pandas import DataFrame

LIST_OF_LIST_OF_STR = str | list[ForwardRef("LIST_OF_LIST_OF_STR")]
# _PIPE_IDX = list[ForwardRef("PIPE_IDX")] | str | tuple[ForwardRef("PIPE_IDX"), ...]
# PIPE_IDX = TypeVar("PIPE_IDX", bound=_PIPE_IDX)
PIPE_IDX = list[str | tuple] | str | tuple
SUB_GROUP = tuple[tuple[str, str], DataFrame]
