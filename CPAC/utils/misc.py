# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Removes functions not used by C-PAC
#     * Removes calls to templateflow
#     * Removes templateflow doctests
#     * Docstrings updated accordingly

# ORIGINAL WORK'S ATTRIBUTION NOTICE:
#    Copyright 2020 The NiPreps Developers
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at

#        http://www.apache.org/licenses/LICENSE-2.0

#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

# Modifications copyright (C) 2019 - 2024  C-PAC Developers
# This file is part of C-PAC.
"""Miscellaneous utilities for Nipype translation of ANTs workflows.

This functionality is adapted from poldracklab/niworkflows:
https://github.com/poldracklab/niworkflows/blob/994dd2dc/niworkflows/utils/misc.py
https://fmriprep.readthedocs.io/
https://poldracklab.stanford.edu/
We are temporarily maintaining our own copy for more granular control.
"""

__all__ = ["get_template_specs"]


def get_template_specs(in_template, template_spec=None, default_resolution=1):
    """Parse template specifications."""
    # Massage spec (start creating if None)
    template_spec = template_spec or {}
    template_spec["desc"] = template_spec.get("desc", None)
    template_spec["atlas"] = template_spec.get("atlas", None)
    template_spec["resolution"] = template_spec.pop(
        "res", template_spec.get("resolution", default_resolution)
    )

    common_spec = {"resolution": template_spec["resolution"]}
    if "cohort" in template_spec:
        common_spec["cohort"] = template_spec["cohort"]


if __name__ == "__main__":
    pass
