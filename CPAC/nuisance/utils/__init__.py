# Copyright (C) 2019-2025  C-PAC Developers

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
"""Utilities for nuisance regression."""

from . import compcor
from .compcor import calc_compcor_components
from .utils import (
    find_offending_time_points,
    generate_summarize_tissue_mask,
    load_censor_tsv,
    NuisanceRegressor,
    temporal_variance_mask,
)

__all__ = [
    "calc_compcor_components",
    "compcor",
    "find_offending_time_points",
    "generate_summarize_tissue_mask",
    "load_censor_tsv",
    "NuisanceRegressor",
    "temporal_variance_mask",
]
