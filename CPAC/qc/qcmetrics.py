# Modifications: Copyright (C) 2022  C-PAC Developers

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

# Original code: BSD 3-Clause License

# Copyright (c) 2020, Lifespan Informatics and Neuroimaging Center

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""QC metrics from XCP-D v0.0.9

Ref: https://github.com/PennLINC/xcp_d/tree/0.0.9
"""
# LGPL-3.0-or-later: Module docstring and lint exclusions
# pylint: disable=invalid-name, redefined-outer-name
# BSD-3-Clause: imports and unspecified sections
import nibabel as nb
import numpy as np


# BSD-3-Clause
def coverage(input1, input2):
    """Estimate the coverage between two masks."""
    input1 = nb.load(input1).get_fdata()
    input2 = nb.load(input2).get_fdata()
    input1 = np.atleast_1d(input1.astype(np.bool))
    input2 = np.atleast_1d(input2.astype(np.bool))
    intsec = np.count_nonzero(input1 & input2)
    if np.sum(input1) > np.sum(input2):
        smallv = np.sum(input2)
    else:
        smallv = np.sum(input1)
    cov = float(intsec)/float(smallv)
    return cov


# BSD-3-Clause
def crosscorr(input1, input2):
    r"""cross correlation: compute cross correction bewteen input masks"""
    input1 = nb.load(input1).get_fdata()
    input2 = nb.load(input2).get_fdata()
    input1 = np.atleast_1d(input1.astype(np.bool)).flatten()
    input2 = np.atleast_1d(input2.astype(np.bool)).flatten()
    cc = np.corrcoef(input1, input2)[0][1]
    return cc


# BSD-3-Clause
def dc(input1, input2):
    r"""
    Dice coefficient
    Computes the Dice coefficient (also known as Sorensen index)
    between the binary objects in two images.
    The metric is defined as

    .. math::
        DC=\frac{2|A\cap B|}{|A|+|B|}

    , where :math:`A` is the first and :math:`B` the second set of
    samples (here: binary objects).

    Parameters
    ----------
    input1 : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.

    input2 : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.

    Returns
    -------
    dc : float
        The Dice coefficient between the object(s) in ```input1``` and the
        object(s) in ```input2```. It ranges from 0 (no overlap) to 1
        (perfect overlap).

    Notes
    -----
    This is a real metric.
    """
    input1 = nb.load(input1).get_fdata()
    input2 = nb.load(input2).get_fdata()
    input1 = np.atleast_1d(input1.astype(np.bool))
    input2 = np.atleast_1d(input2.astype(np.bool))

    intersection = np.count_nonzero(input1 & input2)

    size_i1 = np.count_nonzero(input1)
    size_i2 = np.count_nonzero(input2)

    try:
        dc = 2. * intersection / float(size_i1 + size_i2)
    except ZeroDivisionError:
        dc = 0.0

    return dc


# BSD-3-Clause
def jc(input1, input2):
    r"""
    Jaccard coefficient
    Computes the Jaccard coefficient between the binary objects in two images.
    Parameters
    ----------
    input1: array_like
            Input data containing objects. Can be any type but will be
            converted into binary: background where 0, object everywhere else.
    input2: array_like
            Input data containing objects. Can be any type but will be
            converted into binary: background where 0, object everywhere else.
    Returns
    -------
    jc: float
        The Jaccard coefficient between the object(s) in `input1` and the
        object(s) in `input2`. It ranges from 0 (no overlap) to 1
        (perfect overlap).
    Notes
    -----
    This is a real metric.
    """
    input1 = nb.load(input1).get_fdata()
    input2 = nb.load(input2).get_fdata()
    input1 = np.atleast_1d(input1.astype(np.bool))
    input2 = np.atleast_1d(input2.astype(np.bool))

    intersection = np.count_nonzero(input1 & input2)
    union = np.count_nonzero(input1 | input2)

    jc = float(intersection) / float(union)

    return jc


# LGPL-3.0-or-later
def _prefix_regqc_keys(qc_dict: dict, prefix: str) -> str:
    """Prepend string to each key in a qc dict

    Parameters
    ----------
    qc_dict : dict
        output of ``qc_masks``

    prefix : str
        string to prepend

    Returns
    -------
    dict
    """
    return {f'{prefix}{_key}': _value for _key, _value in qc_dict.items()}


# BSD-3-Clause: logic
# LGPL-3.0-or-later: docstring and refactored function
def qc_masks(registered_mask: str, native_mask: str) -> dict:
    """Return QC measures for coregistration

    Parameters
    ----------
    registered_mask : str
        path to registered mask

    native_mask : str
        path to native-space mask

    Returns
    -------
    dict
    """
    return {'Dice': [dc(registered_mask, native_mask)],
            'Jaccard': [jc(registered_mask, native_mask)],
            'CrossCorr': [crosscorr(registered_mask, native_mask)],
            'Coverage': [coverage(registered_mask, native_mask)]}


# BSD-3-Clause: name and signature
# LGPL-3.0-or-later: docstring and refactored function
def regisQ(bold2t1w_mask: str, t1w_mask: str, bold2template_mask: str,
           template_mask: str) -> dict:
    """Collect coregistration QC measures

    Parameters
    ----------
    bold2t1w_mask, t1w_mask, bold2template_mask, template_mask : str

    Returns
    -------
    dict
    """
    return {**_prefix_regqc_keys(qc_masks(bold2t1w_mask, t1w_mask), 'coreg'),
            **_prefix_regqc_keys(qc_masks(bold2template_mask, template_mask),
                                 'norm')}
