# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Handles local paths
#     * Docstrings updated accordingly
#     * Style modifications

# ORIGINAL WORK'S ATTRIBUTION NOTICE:
#     Copyright (c) 2009-2016, Nipype developers

#     Licensed under the Apache License, Version 2.0 (the "License");
#     you may not use this file except in compliance with the License.
#     You may obtain a copy of the License at

#         http://www.apache.org/licenses/LICENSE-2.0

#     Unless required by applicable law or agreed to in writing, software
#     distributed under the License is distributed on an "AS IS" BASIS,
#     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#     See the License for the specific language governing permissions and
#     limitations under the License.

#     Prior to release 0.12, Nipype was licensed under a BSD license.

# Modifications copyright (C) 2021 - 2024  C-PAC Developers
# This file is part of C-PAC.
"""FSL Nipype interfaces with customized functionality.

Modified from https://github.com/nipy/nipype/blob/f64bf33/nipype/interfaces/fsl/utils.py#L373-L413
"""

import os

from nipype.interfaces.base import isdefined
from nipype.interfaces.fsl.utils import Merge as fslMerge


class Merge(fslMerge):  # noqa: D101
    def _format_arg(self, name, spec, value):
        """Can use relative paths for input files."""
        if name == "tr":
            if self.inputs.dimension != "t":
                msg = "When TR is specified, dimension must be t"
                raise ValueError(msg)
            return spec.argstr % value
        if name == "dimension":
            if isdefined(self.inputs.tr):
                return "-tr"
            return spec.argstr % value
        # Begin custom code
        # -----------------
        if name == "in_files":
            return " ".join([os.path.relpath(in_file) for in_file in value])
        # ---------------
        # End custom code
        return super(Merge, self)._format_arg(name, spec, value)


Merge.__doc__ = fslMerge.__doc__
__all__ = ["Merge"]
