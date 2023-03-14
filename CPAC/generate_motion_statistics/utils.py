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
"""Utilities for motion parameters"""


def affine_from_params(params):
    """Convert a 6-DOF motion parameters array into a 4x4 affine matrix

    Parameters
    ----------
    params : np.ndarray
        2-dimensional array (t x 6) with the first dimension as timepoints and
        the second dimension containing these 6 elements
        (as output by `3dvolreg -1Dfile`):
        roll, pitch, yaw in degrees counterclockwise, and
        displacement in millimeters in the respective
        superior, left and posterior directions

    Returns
    -------
    affine : np.ndarray
        t x 4 x 4 affine matrix with
        the upper-left 3 x 3 matrix encoding rotation,
        the top three values in the fourth column encoding translation
        and the bottom row being [0, 0, 0, 1] for each timepoint
"""
    from nipype.algorithms.rapidart import _get_affine_matrix
    from scipy.spatial.transform import Rotation
    out = []
    for i in range(params.shape[0]):
        affine = _get_affine_matrix(params=params[i], source='AFNI')
        affine[:3, :3] = Rotation.from_euler("ZXY", -params[i, :3],
                                             degrees=True).as_matrix()
        out.append(affine)
    return out
