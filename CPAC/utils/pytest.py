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
"""Utilities for Pytest integration"""
try:
    import pytest
    HAVE_PYTEST = True
except ImportError:
    HAVE_PYTEST = False


def skipif(condition, reason):
    """Skip test if we have Pytest, ignore test entirely if not"""
    def decorator(func):
        """Skip test if we have Pytest"""
        if HAVE_PYTEST:
            return pytest.mark.skipif(condition, reason)(func)
        return func  # return undecorated function
    return decorator  # return conditionally decorated function
