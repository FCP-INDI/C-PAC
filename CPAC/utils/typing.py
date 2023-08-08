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
Helpers and aliases for handling typing in main and variant Python versions

Once all variants (see {DOCS_URL_PREFIX}/user/versions#variants)
run Python ≥ 3.10, these global variables can be replaced with the
current preferred syntax.
"""
import sys
from CPAC.utils.docs import DOCS_URL_PREFIX
__doc__ = __doc__.replace(r'{DOCS_URL_PREFIX}', DOCS_URL_PREFIX)

if sys.version_info >= (3, 10):
    LIST_OR_STR = list | str
    TUPLE = tuple
else:
    from typing import Tuple, Union
    LIST_OR_STR = Union[list, str]
    TUPLE = Tuple

__all__ = ['LIST_OR_STR', 'TUPLE']
