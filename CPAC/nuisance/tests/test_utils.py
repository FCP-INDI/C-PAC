import os
import tempfile
import pkg_resources as p
import numpy as np
from CPAC.nuisance.utils import find_offending_time_points
from CPAC.nuisance.utils import calc_compcor_components

mocked_outputs = \
    p.resource_filename(
        "CPAC",
        os.path.join(
            'nuisance',
            'tests',
            'motion_statistics'
        )
    )


def test_find_offending_time_points():

    dl_dir = tempfile.mkdtemp()
    os.chdir(dl_dir)

    censored = find_offending_time_points(
        os.path.join(mocked_outputs, 'FD_J.1D'),
        os.path.join(mocked_outputs, 'FD_P.1D'),
        os.path.join(mocked_outputs, 'DVARS.1D'),
        2.0,
        2.0,
        '1.5SD'
    )

    censored = np.loadtxt(censored).astype(bool)

    assert set(np.where(np.logical_not(censored))[0].tolist()) == set([1, 3, 7])

def test_calc_compcor_components():

    data_filename = "/cc_dev/cpac_working/old_compcor/nuisance_0_0/_scan_test/_selector_CSF-2mmE-M_aC-WM-2mm-DPC5_G-M_M-SDB_P-2_BP-B0.01-T0.1/Functional_2mm_flirt/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_maths_flirt.nii.gz"
    mask_filename = "/cc_dev/cpac_working/old_compcor/nuisance_0_0/_scan_test/_selector_CSF-2mmE-M_aC-WM-2mm-DPC5_G-M_M-SDB_P-2_BP-B0.01-T0.1/aCompCor_union_masks/segment_seg_2_maths_flirt_mask.nii.gz"


    compcor_filename = calc_compcor_components(data_filename, 5, mask_filename)

    print('compcor components written to {0}'.format(compcor_filename))
    assert 0 == 1
