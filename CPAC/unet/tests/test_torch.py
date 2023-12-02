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
from unittest.mock import MagicMock, patch
import pytest


@pytest.mark.parametrize('readonly', [False, True])
@pytest.mark.parametrize('workdir', [False, True])
def test_import_torch(monkeypatch, readonly, tmp_path, workdir):
    """
    Test that ``torch`` successfully imports after being installed dynamically.

    This test is necessarily slow because it involves dynamically
    installing ``torch``.
    """
    if readonly:
        # set PYTHONUSERBASE to a readonly directory
        os.chmod(tmp_path, 0o444)
        monkeypatch.setenv('PYTHONUSERBASE', tmp_path)
    if workdir:
        os.environ['CPAC_WORKDIR'] = str(tmp_path)
    else:
        if 'CPAC_WORKDIR' in os.environ:
            del os.environ['CPAC_WORKDIR']

    # pylint: disable=import-error,unused-import,wrong-import-order
    from CPAC import unet
    import torch


@pytest.mark.parametrize('error', [ImportError, ModuleNotFoundError, None])
def test_validate_unet(error):
    """Test that pipeline validation throws error if torch is not
    installable"""
    if error:
        import_module = MagicMock(side_effect=error())
        with patch('importlib.import_module', import_module):
            with pytest.raises(OSError) as os_error:
                from CPAC.utils.configuration import Preconfiguration
                monkey = Preconfiguration('monkey')
            assert "U-Net" in str(os_error)
    else:
        from CPAC.utils.configuration import Preconfiguration
        monkey = Preconfiguration('monkey')
        assert 'unet' in [using.lower() for using in monkey[
            'anatomical_preproc', 'brain_extraction', 'using']]
