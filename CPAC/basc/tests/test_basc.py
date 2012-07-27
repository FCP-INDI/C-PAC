import numpy as np

from ..utils import timeseries_bootstrap, \
                    standard_bootstrap, \
                    cluster_timeseries, \
                    adjacency_matrix, \
                    individual_stability_matrix

def test_sample_bootstrap():
    """
    Tests the sample_bootstrap method of BASC workflow
    """
    np.random.seed(27)
    x = np.arange(50).reshape((5,10)).T
    actual = sample_bootstrap(x,3)
    desired = np.array([[ 4, 14, 24, 34, 44],
                       [ 5, 15, 25, 35, 45],
                       [ 6, 16, 26, 36, 46],
                       [ 8, 18, 28, 38, 48],
                       [ 9, 19, 29, 39, 49],
                       [ 0, 10, 20, 30, 40],
                       [ 7, 17, 27, 37, 47],
                       [ 8, 18, 28, 38, 48],
                       [ 9, 19, 29, 39, 49],
                       [ 8, 18, 28, 38, 48]])
    np.testing.assert_equal(actual, desired)
    

def generate_blobs():
    np.random.seed(27)
    x1 = np.random.randn(200,2) + np.array([1.4, 1.8])
    x2 = np.random.randn(100,2) + np.array([4.7, 4.0])
    x3 = np.random.randn(400,2) + np.array([100.7, 100.0])
    blobs = np.vstack((x1,x2,x3))
    return blobs

def test_individual_stability_matrix():
    """
    Tests individual_stability_matrix method on three gaussian blobs.
    """
    
    blobs = generate_blobs()
    ism = individual_stability_matrix(blobs.T, 10, 3)
    
    assert False
    
def test_group_stability_matrix():
    """
    Tests group_stability_matrix method.  This creates a dataset of blobs varying only by additive zero-mean gaussian
    noise and calculates the group stability matrix.
    """
    blobs = generate_blobs()
    
    ism_dataset = np.zeros((5, blobs.shape[0], blobs.shape[0]))
    for i in range(ism_dataset.shape[0]):
        ism_dataset[i] = individual_stability_matrix(blobs.T + 0.2*np.random.randn(blobs.shape[1], blobs.shape[0]), 10, 3)
        
    gsm_stratified = group_stability_matrix(ism_dataset, 10, 100, [0,1,1,1,0])
    gsm = group_stability_matrix(ism_dataset)
    
    assert False
    
def test_basc():
    import glob, os
    g_string = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/func_in_mnioutputs/fmri_mnioutputs/_session_id_NYU_TRT_session1_subject_id_sub*/_csf_threshold_0.4/_gm_threshold_0.2/_wm_threshold_0.66/_run_scrubbing_False/_nc_5/_selector_6.7/apply_warp/mapflow/_apply_warp0/residual_warp.nii.gz'
    roi_file = '/home/data/Projects/nuisance_reliability_paper/seed_files/basil_ganglia/LEFT_BG_3_numbered+tlrc..nii.gz'
    
    subjects_list = glob.glob(g_string)
    b = basc.create_basc()
    b.base_dir = os.getcwd()
    b.inputs.inputspec.roi = roi_file
    b.inputs.inputspec.subjects = subjects_list
    b.inputs.inputspec.k_clusters = 6
    b.inputs.inputspec.dataset_bootstraps = 10
    b.inputs.inputspec.timeseries_bootstraps = 1000
    