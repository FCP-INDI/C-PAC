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
# """Specific typing stubs for bids2table."""
from typing import Literal, Optional, overload

from elbow.typing import StrOrPath
from bids2table.table import BIDSTable

@overload
def bids2table(
    root: StrOrPath,
    *,
    with_meta: bool = True,
    persistent: bool = False,
    index_path: Optional[StrOrPath] = None,
    exclude: Optional[list[str]] = None,
    incremental: bool = False,
    overwrite: bool = False,
    workers: Optional[int] = None,
    worker_id: Optional[int] = None,
    return_table: Literal[True] = True,
) -> BIDSTable: ...
@overload
def bids2table(
    root: StrOrPath,
    *,
    with_meta: bool = True,
    persistent: bool = False,
    index_path: Optional[StrOrPath] = None,
    exclude: Optional[list[str]] = None,
    incremental: bool = False,
    overwrite: bool = False,
    workers: Optional[int] = None,
    worker_id: Optional[int] = None,
    return_table: Literal[False],
) -> None: ...
