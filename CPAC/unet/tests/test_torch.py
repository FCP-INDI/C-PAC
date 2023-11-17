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
"""Test torch installation"""
import os
from pathlib import Path


def test_import_torch():
    """
    Test that ``torch`` successfully imports after being installed dynamically.

    This test is necessarily slow because it involves dynamically
    installing ``torch``.
    """
    # pylint: disable=import-error,unused-import,wrong-import-order
    from CPAC import unet
    import torch


def test_writable_homedir():
    """
    Test that ``/tmp/home/c-pac_user`` should be writable
    """
    home_dir = Path('/tmp/home/c-pac_user')
    owner = home_dir.owner()
    assert owner == 'c-pac_user', \
           f"{home_dir} is owned by {owner} instead of 'c-pac_user'"
    assert os.access(home_dir, os.W_OK), f'{home_dir} is not writable'
