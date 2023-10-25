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
import numpy as np


def affine_file_from_params_file(params_file: str, affine_file: str = None
                                 ) -> str:
    """Convert a 6-DOF motion parameters array into a 4x4 affine matrix.

    Parameters
    ----------
    params_file : str
        path to a motion parameter file (6 DOF, one timepoint per line)

    affine_file : str
        path to a 4x4 affine matrix file (the first 12 values,
        one timepoint per line), optional
        If included, comments will be passed on to filtered affine file

    Returns
    -------
    affine_file : str
        path to a 4x4 affine matrix file (the first 12 values,
        one timepoint per line)
    """
    import os
    import numpy as np
    from CPAC.generate_motion_statistics.utils import affine_from_params

    ## load parameters file into array
    affine = affine_from_params(np.genfromtxt(params_file))
    header = ''
    if affine_file:
        # grab comment(s), if any
        with open(affine_file, 'r', encoding='utf-8') as _f:
            header = '\n'.join([line for line in _f.readlines()
                                if line.lstrip().startswith('#')])
    basename = (os.path.basename(affine_file) if affine_file
                else f'affine_{os.path.basename(params_file)}')
    affine_file = f'{os.getcwd()}/filtered_{basename}'

    # drop bottom [0, 0, 0, 1] row from each matrix
    # insert original comments, if any, as header
    np.savetxt(affine_file, affine[:, :3].reshape(affine.shape[0], 12),
               header=header, comments='')
    return affine_file


def affine_from_params(params: np.ndarray) -> np.ndarray:
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
    return np.array(out)


def load_mats(mat_dir: str) -> np.ndarray:
    """
    Given a directory of affince matrices as output by MCFLIRT,
    return an array of these matrices

    Parameters
    ----------
    mat_dir : str
        path to MCFLIRT matrix output

    Returns
    -------
    np.ndarray
        t x 4 x 4 affine matrix
    """
    from pathlib import Path
    import numpy as np
    mats = []
    mat_paths = sorted(list(Path(mat_dir).glob("MAT_*")))
    for path in mat_paths:
        mat = np.loadtxt(path)
        mats.append(mat)
    mats = np.stack(mats)
    return mats
