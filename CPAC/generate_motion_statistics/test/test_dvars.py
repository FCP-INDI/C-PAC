import os
import nibabel as nb
import numpy as np
from CPAC.generate_motion_statistics import ImageTo1D, calculate_DVARS

np.random.seed(10)


def test_dvars():

    os.chdir('/tmp')

    img_data = np.random.uniform(-3000, 3000, (30, 35, 30, 120))
    img = nb.Nifti1Image(img_data, np.eye(4))
    nb.save(img, 'dvars_data.nii.gz')

    mask_data = np.ones((30, 35, 30))
    mask = nb.Nifti1Image(mask_data, np.eye(4))
    nb.save(mask, 'dvars_mask.nii.gz')

    node = ImageTo1D(method='dvars')
    node.inputs.in_file = 'dvars_data.nii.gz'
    node.inputs.mask = 'dvars_mask.nii.gz'
    node.run()

    calculate_DVARS('dvars_data.nii.gz', 'dvars_mask.nii.gz')

    afni_result = np.loadtxt('dvars_data_3DtoT1.1D')[1:]
    python_result = np.loadtxt('DVARS.txt')

    assert all(np.isclose(afni_result, python_result, 1e-4))
