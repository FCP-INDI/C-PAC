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
"""Tests for requisite surface prerequisites"""
import os
import pytest
from CPAC.utils.tests.test_utils import _installation_check


@pytest.mark.parametrize("executable", ["bc", "csh"])
@pytest.mark.skipif("FREESURFER_HOME" not in os.environ or
                    not os.path.exists(os.environ['FREESURFER_HOME']),
                    reason="We don't need these dependencies if we don't"
                           "have FreeSurfer.")
def test_executable(executable):
    """Make sure executable is installed"""
    _installation_check(executable, "--version")
