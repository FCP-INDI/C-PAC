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
# pylint: disable=ungrouped-imports
try:
    import torch
except (ImportError, ModuleNotFoundError):
    from importlib import invalidate_caches
    from CPAC.info import UNET_REQUIREMENTS
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    invalidate_caches()
    log_subprocess(['pip', 'install', '--user', *UNET_REQUIREMENTS])
    invalidate_caches()
    import torch

__all__ = ['torch']
