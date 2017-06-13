import os
import time
from os.path import expanduser
import matplotlib
import numpy as np
import nibabel as nb
import pandas as pd
import nilearn.image as image
import scipy as sp


from nilearn import datasets
from nilearn.input_data import NiftiMasker
from nilearn.plotting import plot_roi, show
from nilearn.image.image import mean_img
from nilearn.image import resample_img
from matplotlib import pyplot as plt


from sklearn import cluster, datasets
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
#
#import CPAC
#from CPAC.basc.utils import standard_bootstrap, adjacency_matrix, cluster_timeseries, cluster_matrix_average, individual_stability_matrix
#from CPAC.basc import group_stability_matrix, individual_group_clustered_maps, individual_stability_matrix, nifti_individual_stability, ndarray_to_vol, create_basc
#from CPAC.utils import safe_shape

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


matplotlib.style.use('ggplot')
home = expanduser("~")

#from ..utils import timeseries_bootstrap, \
#                    standard_bootstrap, \
#                    cluster_timeseries, \
#                    adjacency_matrix, \
#                    individual_stability_matrix

def test_timeseries_bootstrap():
    """
    Tests the timeseries_bootstrap method of BASC workflow
    """
    np.random.seed(27)
    #np.set_printoptions(threshold=np.nan)
    
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


def test_standard_bootstrap():
    """
    Tests the standard_bootstrap method of BASC workflow
    """
    np.random.seed(27)
    x = np.arange(50).reshape((5,10)).T
    actual = standard_bootstrap(x)
    desired = np.array([[ 3, 13, 23, 33, 43],
                        [ 8, 18, 28, 38, 48],
                        [ 8, 18, 28, 38, 48],
                        [ 8, 18, 28, 38, 48],
                        [ 0, 10, 20, 30, 40],
                        [ 5, 15, 25, 35, 45],
                        [ 8, 18, 28, 38, 48],
                        [ 1, 11, 21, 31, 41],
                        [ 2, 12, 22, 32, 42],
                        [ 1, 11, 21, 31, 41]])
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


def generate_simple_blobs(x):
    np.random.seed(x)
    offset = np.random.randn(30)

    x1 = np.random.randn(200,30) + 2*offset
    x2 = np.random.randn(100,30) + 44*np.random.randn(30)+ 2*offset

    blobs = np.vstack((x1,x2))
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


def test_cross_cluster_timeseries():
    np.random.seed(30)
    offset = np.random.randn(30)
    x1 = np.random.randn(20,30) + 10*offset
    x2 = np.random.randn(10,30) + 44*np.random.randn(30)
    sampledata1 = np.vstack((x1,x2))
    sampledata2 = sampledata1
    actual = cross_cluster_timeseries(sampledata1, sampledata2, 2, 'correlation')
    desired = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1])
    np.testing.assert_equal(actual,desired)


def test_individual_stability_matrix():
    """
    Tests individual_stability_matrix method on three gaussian blobs.
    """

    blobs = generate_blobs()
    ism = individual_stability_matrix(blobs, 10, 3)

    assert False

def test_cross_cluster_individual_stability_matrix():
    """
    Tests individual_stability_matrix method on three gaussian blobs.
    """

    blobs1 = generate_simple_blobs(27)
    blobs2 = generate_simple_blobs(27)
    blobs2 = blobs2[0:150,:]
    ism = individual_stability_matrix(blobs1, 10, 2, Y2 = blobs2, cross_cluster = True)

    assert False

def test_nifti_individual_stability():

    #subject_file= home + '/Dropbox/1_Projects/1_Research/2_BASC/Data/Test_Data/residual_antswarp.nii.gz'
    subject_file= home + '/C-PAC/CPAC/basc/sampledata/subjects/sub1/Func_Quarter_Res.nii.gz'
    roi_mask_file= home + '/C-PAC/CPAC/basc/sampledata/masks/LC_Quarter_Res.nii.gz'
    #roi_mask_file= home + '/Dropbox/1_Projects/1_Research/2_BASC/Data/BasalGanglia_MNI2mm/Left_Caudate.nii.gz'
    n_bootstraps=100
    n_clusters=2
    cross_cluster=True
    roi2_mask_file= home + '/C-PAC/CPAC/basc/sampledata/masks/RC_Quarter_Res.nii.gz'
    #roi2_mask_file= home + '/Dropbox/1_Projects/1_Research/2_BASC/Data/BasalGanglia_MNI2mm/Right_Caudate.nii.gz'
    cbb_block_size=None
    affinity_threshold=0.5
    nifti_individual_stability(subject_file, roi_mask_file, n_bootstraps, n_clusters, cross_cluster, roi2_mask_file, cbb_block_size, affinity_threshold)


def new_test_group_stability_matrix():
    """
    Tests group_stability_matrix method.  This creates a dataset of blobs varying only by additive zero-mean gaussian
    noise and calculates the group stability matrix.
    """
    blobs = generate_blobs()

    ism_dataset = np.zeros((5, blobs.shape[0], blobs.shape[0]))
    ism_list = []
    for i in range(ism_dataset.shape[0]):
        ism_dataset[i] = individual_stability_matrix(blobs + 0.2*np.random.randn(blobs.shape[0], blobs.shape[1]), 10, 3, affinity_threshold = 0.0)
        f = 'ism_dataset_%i.npy' % i
        ism_list.append(f)
        np.save(f, ism_dataset[i])

    indiv_stability_list=ism_list
    n_bootstraps=10
    n_clusters=3
    G, cluster_G, cluster_voxel_scores = group_stability_matrix(ism_list, 10, 3, [0,1,1,1,0])
    G, cluster_G, cluster_voxel_scores = group_stability_matrix(ism_list, 10, 3)

#def group_stability_matrix(indiv_stability_list, n_bootstraps, k_clusters, stratification=None):
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
        ism_dataset[i] = individual_stability_matrix(blobs + 0.2*np.random.randn(blobs.shape[0], blobs.shape[1]), 10, 3, affinity_threshold = 0.0)
        f = 'ism_dataset_%i.npy' % i
        ism_list.append(f)
        np.save(f, ism_dataset[i])

    G, cluster_G, cluster_voxel_scores = group_stability_matrix(ism_list, 10, 3, [0,1,1,1,0])

    return G, cluster_g, cluster_voxel_scores


def test_basc_workflow_runner():

    from basc_workflow_runner import run_basc_workflow
    import utils
    subject_file_list= [home + '/C-PAC/CPAC/basc/sampledata/subjects/sub1/Func_Quarter_Res.nii.gz',
                        home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
                        home + '/C-PAC/CPAC/basc/sampledata/subjects/sub3/Func_Quarter_Res.nii.gz']

    roi_mask_file= home + '/C-PAC/CPAC/basc/sampledata/masks/LC_Quarter_Res.nii.gz'
    dataset_bootstraps=50
    timeseries_bootstraps=10
    n_clusters=2
    cross_cluster=True
    roi2_mask_file= home + '/C-PAC/CPAC/basc/sampledata/masks/RC_Quarter_Res.nii.gz'
    affinity_threshold= [0.5, 0.5, 0.5]
    out_dir= home + '/BASC_outputs'
    run=True
    
    

    basc_test= run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)


def NKI_Ned_test():
    
    from basc_workflow_runner import run_basc_workflow
    import utils
    subject_file_list_Ned= ['/data/rockland_sample/A00060603/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00060503/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00060429/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00060384/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00060280/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00059935/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00059875/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00059734/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00059733/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz']
    
    subject_file_list=['/Users/aki.nikolaidis/BGDev_SampleData/A00060846/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/Users/aki.nikolaidis/BGDev_SampleData/A00060603/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/Users/aki.nikolaidis/BGDev_SampleData/A00060503/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/Users/aki.nikolaidis/BGDev_SampleData/A00060429/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/Users/aki.nikolaidis/BGDev_SampleData/A00060384/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/Users/aki.nikolaidis/BGDev_SampleData/A00060280/bandpassed_demeaned_filtered_antswarp.nii.gz']

    roi_mask_file='/Users/aki.nikolaidis/C-PAC/CPAC/basc/sampledata/masks/BG.nii.gz'
    roi2_mask_file='/Users/aki.nikolaidis/C-PAC/CPAC/basc/sampledata/masks/yeo_2.nii.gz'
    output_size = 2000


    dataset_bootstraps=10
    timeseries_bootstraps=10
    n_clusters=2
    cross_cluster=True
    affinity_threshold= [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    out_dir= home + '/BASC_outputs'
    run=True
    
    basc_test= run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, output_size, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)


def bruteforce_workflow_test():

    #Updates to bruteforce test
    #Change  functional data to higher resolution
    #Change ROIs to larger ones.
    #work out the transformation of the Yeo to the correct size

#    subject_file_list=  [home + '/C-PAC/CPAC/basc/sampledata/subjects/sub1/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub2/Func_Quarter_Res.nii.gz',
#                         home + '/C-PAC/CPAC/basc/sampledata/subjects/sub3/Func_Quarter_Res.nii.gz']
#    
#    /Users/aki.nikolaidis/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz
    
    import time

    subject_file_list = [home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz',
                         home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz']
    
    sample_file = home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/fixedfunc_2thirds_res.nii.gz'
    filename = home + '/C-PAC/CPAC/basc/sampledata/Striatum_GroupLevel_MotorCluster.nii.gz'
    
    

    roi_mask_file= home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/Striatum_2thirdsRes.nii.gz'
    dataset_bootstraps=50
    timeseries_bootstraps=50
    n_clusters=4
    cross_cluster=True
    roi2_mask_file= home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/Yeo_LowRes/yeo_2_2thirdsRes_bin.nii.gz'
    output_size = 500

    ism_list = []
    ism_dataset = np.zeros((len(subject_file_list), output_size, output_size))

    cbb_block_size=None
    affinity_threshold=0.5
    n_bootstraps=timeseries_bootstraps

    roi_mask_file_nb = nb.load(roi_mask_file).get_data().astype('float32').astype('bool')
    roi2_mask_file_nb = nb.load(roi2_mask_file).get_data().astype('float32').astype('bool')


    start = time.time()
    for i in range(len(subject_file_list)):
        data = nb.load(subject_file_list[int(i)]).get_data().astype('float32')
        roi1data = data[roi_mask_file_nb]
        roi2data = data[roi2_mask_file_nb]


        Y1_compressed = data_compression(roi1data, roi_mask_file_nb, output_size).T
        Y2_compressed = data_compression(roi2data, roi2_mask_file_nb, output_size).T
        
        print '(%i voxels, %i timepoints and %i bootstraps)' % (Y1_compressed.shape[0], Y1_compressed.shape[1], n_bootstraps)
        print '(%i voxels, %i timepoints and %i bootstraps)' % (Y2_compressed.shape[0], Y2_compressed.shape[1], n_bootstraps)

        
        ism_dataset[int(i)] = individual_stability_matrix(Y1_compressed, n_bootstraps, n_clusters, Y2_compressed, cross_cluster, cbb_block_size, affinity_threshold)

        #ism_dataset[i] = individual_stability_matrix(blobs.T + 0.2*np.random.randn(blobs.shape[1], blobs.shape[0]), 10, 3, affinity_threshold = 0.0)
        f = home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Results/Testing/ism_dataset_%i.npy' % i
        ism_list.append(f)
        np.save(f, ism_dataset[i])

    print((time.time() - start))
    G, cluster_G, cluster_voxel_scores = group_stability_matrix(ism_list, dataset_bootstraps, n_clusters=n_clusters)
    #Need to figure out how to get low res cluster version outputted in normal space- or how to output 
   # ndarray_to_vol(cluster_G, roi_mask_file, sample_file, filename)
    #loop over cluster_voxel_scores[i] and save each cluster map to nifti file.

def test_basc():
    import glob, os
    g_string = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/func_in_mnioutputs/fmri_mnioutputs/_session_id_NYU_TRT_session1_subject_id_sub*/_csf_threshold_0.4/_gm_threshold_0.2/_wm_threshold_0.66/_run_scrubbing_False/_nc_5/_selector_6.7/apply_warp/mapflow/_apply_warp0/residual_warp.nii.gz'
    roi_file = '/home/data/Projects/nuisance_reliability_paper/seed_files/basil_ganglia/LEFT_BG_3_numbered+tlrc..nii.gz'

    subjects_list = glob.glob(g_string)
    b = basc.create_basc()
    b.base_dir = os.getcwd()
    b.inputs.inputspec.roi = roi_file
    b.inputs.inputspec.subjects = subjects_list
    b.inputs.inputspec.n_clusters = 6
    b.inputs.inputspec.dataset_bootstraps = 10
    b.inputs.inputspec.timeseries_bootstraps = 1000
