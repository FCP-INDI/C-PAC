import numpy as np

from ..utils import timeseries_bootstrap, \
                    standard_bootstrap, \
                    cluster_timeseries, \
                    adjacency_matrix, \
                    individual_stability_matrix

def test_timeseries_bootstrap():
    """
    Tests the timeseries_bootstrap method of BASC workflow
    """
    np.random.seed(27)
    # Create a 10x5 matrix which counts up by column-wise
    x = np.arange(50).reshape((5,10)).T
    actual = timeseries_bootstrap(x,3)
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

def test_adjacency_matrix():
    """
    Tests the adjacency_matrix of BASC workflow
    """
    x = np.asarray([1, 2, 2, 3, 1])[:,np.newaxis]
    actual = adjacency_matrix(x).astype(int)
    desired = np.array([[1, 0, 0, 0, 1],
                       [0, 1, 1, 0, 0],
                       [0, 1, 1, 0, 0],
                       [0, 0, 0, 1, 0],
                       [1, 0, 0, 0, 1]])
    np.testing.assert_equal(actual, desired)
    
def generate_blobs():
    np.random.seed(27)
    offset = np.random.randn(30)
    
    x1 = np.random.randn(200,30) + 2*offset
    x2 = np.random.randn(100,30) + 44*np.random.randn(30)
    x3 = np.random.randn(400,30)
    blobs = np.vstack((x1,x2,x3))
    return blobs

def generate_blobs_3d():
    np.random.seed(27)
    x1 = np.random.randn(200,3) + np.array([1.4, 1.8, 22.2])
    x2 = np.random.randn(100,3) + np.array([4.7, 4.0, 9.6])
    x3 = np.random.randn(400,3) + np.array([100.7, 100.0, 100.8])
    blobs = np.vstack((x1,x2,x3))
    return blobs
    
def test_cluster_timeseries():
    """
    Tests the cluster_timeseries method on three blobs in three dimensions (to make correlation possible)
    """
    blobs = generate_blobs_3d()
    y_predict = cluster_timeseries(blobs, 3, similarity_metric = 'correlation')

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
    ism_list = []
    for i in range(ism_dataset.shape[0]):
        ism_dataset[i] = individual_stability_matrix(blobs.T + 0.2*np.random.randn(blobs.shape[1], blobs.shape[0]), 10, 3, affinity_threshold = 0.0)
        f = 'ism_dataset_%i.npy' % i
        ism_list.append(f)
        np.save(f, ism_dataset[i])
    
    G, cluster_G, cluster_voxel_scores = group_stability_matrix(ism_list, 10, 3, [0,1,1,1,0])

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
    