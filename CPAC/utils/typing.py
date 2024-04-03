# Copyright (C) 2023  C-PAC Developers

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
"""
Helpers and aliases for handling typing in main and variant Python versions.

Once all variants (see {DOCS_URL_PREFIX}/user/versions#variants)
run Python ≥ 3.10, these global variables can be replaced with the
current preferred syntax.
"""
from pathlib import Path
import sys
from typing import Union

from CPAC.utils.docs import DOCS_URL_PREFIX

# Set the version-specific documentation URL in the module docstring:
__doc__ = __doc__.replace(r"{DOCS_URL_PREFIX}", DOCS_URL_PREFIX)  # noqa: A001

if sys.version_info >= (3, 8):
    from typing import Literal as LITERAL
else:
    from typing_extensions import Literal as LITERAL
if sys.version_info >= (3, 9):
    from builtins import dict as DICT, list as LIST
    from collections.abc import Iterable as ITERABLE
else:
    from typing import Dict as DICT, Iterable as ITERABLE, List as LIST
if sys.version_info >= (3, 10):
    from builtins import tuple as TUPLE

    LIST_OR_STR = LIST[str] | str  # pylint: disable=invalid-name
else:
    from typing import Tuple as TUPLE

    LIST_OR_STR = Union[LIST[str], str]  # pylint: disable=invalid-name

PATHSTR = Union[Path, str]
ConfigKeyType = Union[str, LIST[str]]
__all__ = [
    "ConfigKeyType",
    "DICT",
    "ITERABLE",
    "LIST",
    "LIST_OR_STR",
    "LITERAL",
    "PATHSTR",
    "TUPLE",
]
