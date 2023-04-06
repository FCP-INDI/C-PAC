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
"""Functional preprocessing"""
from .func_motion import calc_motion_stats, func_motion_correct, \
                         func_motion_correct_only, func_motion_estimates, \
                         get_motion_ref, motion_estimate_filter
from .func_preproc import slice_timing_wf, \
                          get_idx

__all__ = ['calc_motion_stats', 'func_motion_correct',
           'func_motion_correct_only', 'func_motion_estimates', 'get_idx',
           'get_motion_ref', 'motion_estimate_filter', 'slice_timing_wf']
