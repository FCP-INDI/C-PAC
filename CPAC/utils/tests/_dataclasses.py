"""Dataclasses for pickled test data"""
from dataclasses import dataclass
import nibabel as nib
import numpy as np

__all__ = ['_FdjTestData', '_MockImageHeaderOnly']


@dataclass
class _FdjTestData:
    """Class for storing test data for FD-J functions"""
    affine: np.ndarray
    img: nib.Nifti1Image
    rels_rms: np.ndarray
