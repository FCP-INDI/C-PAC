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
"""Dynamically install torch iff we're going to use it"""
# pylint: disable=import-error,redefined-outer-name,ungrouped-imports,unused-import
from importlib import invalidate_caches
import os
import site
from subprocess import CalledProcessError
import sys
from typing import Optional
from CPAC.info import UNET_REQUIREMENTS
from CPAC.utils.monitoring.custom_logging import log_subprocess


def _custom_pip_install(env_var: Optional[str] = None) -> None:
    """
    ``pip install --user torch``, in a custom ``--user`` directory if
    one is provided, then ``import torch``

    Parameters
    ----------
    env_var : str, optional
        environment variable to define custom ``--user`` directory
    """
    if env_var is not None:
        if env_var not in os.environ:
            raise FileNotFoundError(f'${env_var}')
        site.USER_BASE = os.environ['PYTHONUSERBASE'] = os.path.join(
            os.environ[env_var], '.local')
        py_version = '.'.join(str(getattr(sys.version_info, attr)) for
                              attr in ['major', 'minor'])
        pythonpath = f'{site.USER_BASE}/lib/python{py_version}/site-packages'
        sys.path.append(pythonpath)
        os.environ['PYTHONPATH'] = ':'.join([os.environ['PYTHONPATH'],
                                             pythonpath]).replace('::', ':')
    invalidate_caches()
    log_subprocess(['pip', 'install', '--user', *UNET_REQUIREMENTS])
    invalidate_caches()
    import torch


try:
    import torch  # just import if it's available
except (ImportError, ModuleNotFoundError):
    torch = NotImplemented
    try:
        _custom_pip_install()  # pip install in default user directory
    except (CalledProcessError, FileNotFoundError, ImportError,
            ModuleNotFoundError, OSError):
        try:
            _custom_pip_install('CPAC_WORKDIR')  # pip install in $CPAC_WORKDIR
        except (CalledProcessError, FileNotFoundError, ImportError,
                ModuleNotFoundError, OSError):
            _custom_pip_install('PWD')  # pip install in $PWD
if torch is not NotImplemented:
    __all__ = ['torch']
