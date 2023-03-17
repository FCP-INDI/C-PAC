# Copyright (C) 2012-2023  C-PAC Developers

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
"""Functions for generating motion statistics"""
from .generate_motion_statistics import (motion_power_statistics,
                                        calculate_FD_P,
                                        calculate_FD_J,
                                        gen_motion_parameters,
                                        gen_power_parameters,
                                        calculate_DVARS,
                                        ImageTo1D)
from .utils import affine_file_from_params_file, affine_from_params

__all__ = [
    'affine_file_from_params_file',
    'affine_from_params',
    'calculate_DVARS',
    'calculate_FD_P',
    'calculate_FD_J',
    'gen_motion_parameters',
    'gen_power_parameters',
    'ImageTo1D',
    'motion_power_statistics'
]
