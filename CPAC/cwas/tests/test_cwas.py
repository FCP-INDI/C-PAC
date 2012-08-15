import numpy as np
from ..cwas import nifti_cwas

def test_nifti_cwas():

def test_cwas():
    from CPAC.cwas import create_cwas
    import numpy as np
    import os, glob

    g_string = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/func_in_mnioutputs/fmri_mnioutputs/_session_id_NYU_TRT_session1_subject_id_sub*/_csf_threshold_0.4/_gm_threshold_0.2/_wm_threshold_0.66/_run_scrubbing_False/_nc_5/_selector_6.7/apply_warp/mapflow/_apply_warp0/residual_warp.nii.gz'
    roi_file = '/home/data/Projects/nuisance_reliability_paper/seed_files/basil_ganglia/LEFT_BG_3_numbered+tlrc..nii.gz'
    subjects_list = glob.glob(g_string)
    
    c = create_cwas()
    c.base_dir = os.getcwd()
    c.inputs.inputspec.roi = roi_file
    c.inputs.inputspec.subjects = subjects_list[:4]
    c.inputs.inputspec.regressor = np.arange(len(subjects_list[:4]))[:,np.newaxis]
    c.inputs.inputspec.parallel_nodes = 3
    c.run()