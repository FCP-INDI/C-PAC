"""Convert formats"""
import os

import numpy as np


def one_d_to_mat(one_d_filename):
    """Convert a .1D file to a .mat directory

    Parameters
    ----------
    one_d_filename : str
        The filename of the .1D file to convert

    Returns
    -------
    mat_dirname : str
    """
    mat_dirname = one_d_filename.replace('.1D', '.mat')
    with open(one_d_filename, 'r') as one_d_file:
        rows = [np.reshape(row, (4, 4)).astype('float') for row in [[
            term.strip() for term in row.split(' ') if term.strip()
        ] + [0, 0, 0, 1] for row in [
            line.strip() for line in one_d_file.readlines() if
            not line.startswith('#')]]]
    try:
        os.mkdir(mat_dirname)
    except FileExistsError:
        pass
    for i, row in enumerate(rows):
        np.savetxt(os.path.join(mat_dirname, f'MAT_{i:04}'),
                   row, fmt='%.5f', delimiter=' ')
    return mat_dirname
