# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Removes interfaces and functions we're not using in C-PAC
#     * Adjusts imports for C-PAC
#     * Takes final value from list instead of trying to cast list to float if `min_val` is a list
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

# Modifications Copyright (C) 2019-2024  C-PAC Developers

# This file is part of C-PAC.
"""Customization of BrickStat.

Modified from https://github.com/nipy/nipype/blob/f64bf33/nipype/interfaces/afni/utils.py#L280-L330
"""

import os

from nipype.interfaces.afni.utils import BrickStat as NipypeBrickStat
from nipype.utils.filemanip import load_json, save_json


class BrickStat(NipypeBrickStat):  # noqa: D101
    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        """Populate outputs."""
        outputs = self._outputs()

        outfile = os.path.join(os.getcwd(), "stat_result.json")

        if runtime is None:
            try:
                min_val = load_json(outfile)["stat"]
            except IOError:
                return self.run().outputs
        else:
            min_val = []
            for line in runtime.stdout.split("\n"):
                if line:
                    values = line.split()
                    if len(values) > 1:
                        min_val.append([float(val) for val in values])
                    else:
                        min_val.extend([float(val) for val in values])

            if len(min_val) == 1:
                min_val = min_val[0]
            save_json(outfile, {"stat": min_val})
        outputs.min_val = min_val[-1] if isinstance(min_val, list) else min_val

        return outputs


BrickStat.__doc__ = NipypeBrickStat.__doc__
__all__ = ["BrickStat"]
