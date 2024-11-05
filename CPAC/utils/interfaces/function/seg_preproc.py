# Copyright (C) 2022-2023  C-PAC Developers

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
"""Function interfaces for seg_preproc."""

from CPAC.utils.interfaces.function.function import Function


def pick_tissue_from_labels_file_interface(input_names=None):
    """Create a Function interface for ~CPAC.seg_preproc.utils.pick_tissue_from_labels_file.

    Parameters
    ----------
    input_names : list, optional

    Returns
    -------
    nipype.interfaces.base.core.Interface
    """
    # pylint: disable=import-outside-toplevel
    from CPAC.seg_preproc.utils import pick_tissue_from_labels_file

    if input_names is None:
        input_names = ["multiatlas_Labels", "csf_label", "gm_label", "wm_label"]
    return Function(
        input_names=input_names,
        output_names=["csf_mask", "gm_mask", "wm_mask"],
        function=pick_tissue_from_labels_file,
    )
