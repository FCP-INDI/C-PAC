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

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


matplotlib.style.use('ggplot')
home = expanduser("~")
proc_mem= [2,4] #first is number of proc , second total number of mem
    
#%% TEST UTILS.PY
#Remaining Tests:
    #compare_stability_matrix
    #expand_ism
    #data_compression

def test_timeseries_bootstrap():
    """
    Tests the timeseries_bootstrap method of BASC workflow
    """
    np.random.seed(27)
    #np.set_printoptions(threshold=np.nan)
    
    # Create a 10x5 matrix which counts up by column-wise
    x = np.arange(50).reshape((5,10)).T
    actual, other= timeseries_bootstrap(x,3)
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
    roi_mask_nparray='empty'
    blobs = generate_blobs_3d()
    n_clusters=2
    y_predict = cluster_timeseries(blobs, roi_mask_nparray, n_clusters, similarity_metric = 'correlation', affinity_threshold=0.0, neighbors=10)


def test_cross_cluster_timeseries():
    np.random.seed(30)
    offset = np.random.randn(30)
    x1 = np.random.randn(20,30) + 10*offset
    x2 = np.random.randn(10,30) + 44*np.random.randn(30)
    data1 = np.vstack((x1,x2))
    data2 = data1
    actual = cross_cluster_timeseries(data1, data2, n_clusters=2, similarity_metric='correlation', affinity_threshold=0.0)
    desired = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1])
    np.testing.assert_equal(actual,desired)
    print('Correlation equals ', 1-sp.spatial.distance.correlation(actual,desired))


def test_individual_stability_matrix():
    """
    Tests individual_stability_matrix method on three gaussian blobs.
    """
    import utils
    import numpy as np
    import scipy as sp
    desired = np.load(home + '/git_repo/PyBASC/tests/ism_test.npy')
    blobs = generate_blobs()
    ism = utils.individual_stability_matrix(blobs, 20, 3, similarity_metric='correlation')
    #how to use test here?
#    np.corrcoef(ism.flatten(),desired.flatten())
#    np.testing.assert_equal(ism,desired)
#    
#    corr=np.array(sp.spatial.distance.cdist(ism, desired, metric = 'correlation'))
#    
    assert False

def test_cross_cluster_individual_stability_matrix():
    """
    Tests individual_stability_matrix method on three gaussian blobs.
    """

    blobs1 = generate_simple_blobs(27)
    blobs2 = generate_simple_blobs(27)
    blobs2 = blobs2[0:150,:]
    ism = individual_stability_matrix(blobs1, 10, 2, Y2 = blobs2, cross_cluster = True)

    return ism

    
def test_expand_ism_options():
    
    import time
    import random
    import pandas as pd
    import numpy as np
    
    n=60
    k=55
    i=0
    vec1=[]
    for x in range(0, n):
       vec1.append(random.randint(0, k-1))
    
    temp=np.random.random((k,k))
    vec1=np.array(vec1)
    sizevec1=len(vec1)
    matshape=(sizevec1,sizevec1)
    target_mat=np.zeros(matshape)
    
    
    source_mat=temp*temp.T
    np.fill_diagonal(source_mat, 1)
    transform_mat=np.zeros((len(source_mat),len(target_mat)))
    
    
    #Slow Solution
    matrixtime = time.time()
    for row in range(0,target_mat.shape[0]):
        #print 'row is ', row
        for column in range(0,target_mat.shape[1]):
    
            #print 'column is', column
            if (row == column):
                target_mat[row,column]=1
            else:
                target_mat[row,column] = source_mat.item(int(vec1[row]), int(vec1[column]))
            
    print((time.time() - matrixtime))
    target_mat_slow=target_mat
    
    #XU MACKENZIE SOLUTION

    target_mat=np.zeros(matshape)
    matrixtime = time.time()

    for i in range(0,len(target_mat)):
      transform_mat[vec1[i],i]=1
    
    temp=np.dot(source_mat,transform_mat)
    target_mat=np.dot(temp.T,transform_mat)
    target_mat_XM=target_mat
    target_mat=np.zeros(matshape)
    XM_time= time.time() - matrixtime
    print((time.time() - matrixtime))
    
    #Older 'fast' solution
    matrixtime = time.time()
    for row in range(0,source_mat.shape[0]):
        #print('row is ', row)
        #for column in range(0, source_mat.shape[1]):
        for column in range(0, row):   
            rowmatch = np.array([vec1==row])
            rowmatch = rowmatch*1
    
            colmatch = np.array([vec1==column])
            colmatch = colmatch*1
    
            match_matrix=rowmatch*colmatch.T
            target_mat=target_mat+(match_matrix*source_mat[row,column])
    
    print((time.time() - matrixtime))
    target_mat_fast=target_mat
    target_mat=np.zeros(matshape)
    
    target_mat_slow==target_mat_fast
    target_mat_fast==target_mat_XM
    target_mat_slow==target_mat_XM
    
def test_data_compress_expand():
    
    import os
    import numpy as np
    import nibabel as nb
    import utils
    import pandas as pd
    import sklearn as sk
    
    #Setup
    subject_file = home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz'
    roi_mask_file= home + '/git_repo/PyBASC/masks/LC_Quarter_Res.nii.gz'
    roi2_mask_file= home + '/git_repo/PyBASC/masks/RC_Quarter_Res.nii.gz'
    n_bootstraps=100
    n_clusters=10
    output_size=20
    cross_cluster=True
    cbb_block_size=None
    affinity_threshold=0.5
    
    print( 'Calculating individual stability matrix of:', subject_file)


    data = nb.load(subject_file).get_data().astype('float32')
    print( 'Data Loaded')

    if (roi2_mask_file != None):
        print( 'Setting up NIS')
        roi_mask_file_nb = nb.load(roi_mask_file)
        roi2_mask_file_nb= nb.load(roi2_mask_file)

        roi_mask_nparray = nb.load(roi_mask_file).get_data().astype('float32').astype('bool')
        roi2_mask_nparray = nb.load(roi2_mask_file).get_data().astype('float32').astype('bool')


        roi1data = data[roi_mask_nparray]
        roi2data = data[roi2_mask_nparray]
        
        #add code that uploads the roi1data and roi2data, divides by the mean and standard deviation of the timeseries
        roi1data=sk.preprocessing.normalize(roi1data, norm='l2')
        roi2data=sk.preprocessing.normalize(roi2data, norm='l2')
        
        print( 'Compressing data')
        data_dict1 = utils.data_compression(roi1data.T, roi_mask_file_nb, roi_mask_nparray, output_size)
        Y1_compressed = data_dict1['data']
        Y1_compressed = Y1_compressed.T
        Y1_labels = pd.DataFrame(data_dict1['labels'])
        Y1_labels=np.array(Y1_labels)
        print( 'Y1 compressed')
        
        print( 'Compressing Y2')

        data_dict2 = utils.data_compression(roi2data.T, roi2_mask_file_nb, roi2_mask_nparray, output_size)
        Y2_compressed = data_dict2['data']
        Y2_compressed=Y2_compressed.T
        Y2_labels = pd.DataFrame(data_dict2['labels'])
        print( 'Y2 compressed')
        
        print('Going into ism')
        ism = utils.individual_stability_matrix(Y1_compressed, n_bootstraps, n_clusters, Y2_compressed, cross_cluster, cbb_block_size, affinity_threshold)
        #ism=ism/n_bootstraps # was already done in ism

        
        print('Expanding ism')
        voxel_num=roi1data.shape[0]
        voxel_ism = utils.expand_ism(ism, Y1_labels)
      
        
        #voxel_ism=voxel_ism*100 # was already done in ism
        voxel_ism=voxel_ism.astype("uint8")
    


def test_nifti_individual_stability():

    subject_file = home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz'

    roi_mask_file= home + '/git_repo/PyBASC/masks/LC_Quarter_Res.nii.gz'
    roi2_mask_file= None#home + '/git_repo/PyBASC/masks/RC_Quarter_Res.nii.gz'
    
   
    n_bootstraps=100
    n_clusters=2
    output_size=20
    cross_cluster=False
    similarity_metric='correlation'
    cbb_block_size=None
    affinity_threshold=0.5
    nifti_individual_stability(subject_file, roi_mask_file, n_bootstraps, n_clusters, output_size, similarity_metric, cross_cluster, roi2_mask_file, cbb_block_size, affinity_threshold)


def test_cluster_matrix_average():
    
    import utils
    import basc
    import matplotlib.pyplot as plt
    
    roi_mask_nparray='empty'
    blobs = generate_blobs()
    n_clusters=3
    similarity_metric='correlation'
    ism = utils.individual_stability_matrix(blobs, 100, n_clusters, similarity_metric)

    y_predict = utils.cluster_timeseries(blobs, roi_mask_nparray, n_clusters, similarity_metric, affinity_threshold=0.0, neighbors = 10)
    cluster_voxel_scores, K_mask = utils.cluster_matrix_average(ism, y_predict)


    
    plt.imshow(K_mask)
    

#%% TEST BASC.PY
#Remaining Tests to write:
    #Join_group_stability
    #cluster_selection
    #individual_group_clustered_maps
    #ndarray_to_vol
    
def new_test_group_stability_matrix():
    """
    Tests group_stability_matrix method.  This creates a dataset of blobs varying only by additive zero-mean gaussian
    noise and calculates the group stability matrix.
    """
    
    import utils
    import basc
    
    bootstrap=20
    blobs = generate_blobs()

    ism_dataset = np.zeros((5, blobs.shape[0], blobs.shape[0]))
    
    indiv_stability_list=[]
    
    for i in range(ism_dataset.shape[0]):
        ism_dataset[i] = utils.individual_stability_matrix(blobs + 0.2*np.random.randn(blobs.shape[0], blobs.shape[1]), bootstrap, 3, similarity_metric='correlation',affinity_threshold = 0.0)
        f = 'ism_dataset_%i.npy' % i
        indiv_stability_list.append(f)
        np.save(f, ism_dataset[i])

    #indiv_stability_list=ism_list
    n_bootstraps=10
    n_clusters=3

    
    G = basc.map_group_stability(indiv_stability_list, n_bootstraps, n_clusters)
    
    return G
    
  

def test_group_stability_matrix():
    """
    Tests group_stability_matrix method.  This creates a dataset of blobs varying only by additive zero-mean gaussian
    noise and calculates the group stability matrix.
    """
    #def map_group_stability(indiv_stability_list, n_clusters, bootstrap_list, stratification=None):

    
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


def test_individual_group_clustered_maps():
#    indiv_stability_list
#    clusters_G 
#    roi_mask_file
#
#    import utils
#    import basc
#    
#    bootstrap=20
#    blobs = generate_blobs()
#
#    ism_dataset = np.zeros((5, blobs.shape[0], blobs.shape[0]))
#    
#    indiv_stability_list=[]
#    
#    for i in range(ism_dataset.shape[0]):
#        ism_dataset[i] = utils.individual_stability_matrix(blobs + 0.2*np.random.randn(blobs.shape[0], blobs.shape[1]), 10, 3, affinity_threshold = 0.0)
#        f = 'ism_dataset_%i.npy' % i
#        indiv_stability_list.append(f)
#        np.save(f, ism_dataset[i])
#
#    G, cluster_G, cluster_voxel_scores = group_stability_matrix(ism_list, 10, 3, [0,1,1,1,0])
    
    
    
    import basc
    import utils
    subject_file_list= [home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub3/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz']

    roi_mask_file= home + '/git_repo/PyBASC/masks/LC_Quarter_Res.nii.gz'
    dataset_bootstraps=5
    timeseries_bootstraps=5
    n_clusters=3
    output_size=10
    similarity_metric = 'correlation'
    bootstrap_list=list(range(0,dataset_bootstraps))
    cross_cluster=True
    blocklength=1
    roi2_mask_file= home + '/git_repo/PyBASC/masks/RC_Quarter_Res.nii.gz'
    cbb_block_size=None
    affinity_threshold= 0.5 #* len(subject_file_list)
    out_dir= home + '/PyBASC_outputs/IGCM_HowsItWork3'
    run=True
    indiv_stability_list=[]
    for i in range(0,len(subject_file_list)):
        temp = basc.nifti_individual_stability(subject_file_list[i], roi_mask_file, timeseries_bootstraps, n_clusters, output_size,similarity_metric, cross_cluster, roi2_mask_file, blocklength, cbb_block_size, affinity_threshold)
        #def nifti_individual_stability(subject_file, roi_mask_file, n_bootstraps, n_clusters, output_size, similarity_metric, cross_cluster=False, roi2_mask_file=None, blocklength=1, cbb_block_size=None, affinity_threshold=0.5):

        #temp=temp/timeseries_bootstraps
        indiv_stability_list.append(temp)
        
    G_file=[]
    for i in range(0,dataset_bootstraps):
        temp2= map_group_stability(indiv_stability_list, n_clusters, bootstrap_list, roi_mask_file)
        G_file.append(temp2)
        
    G, clusters_G, ism_gsm_corr, gsm_file, clusters_G_file, ism_gsm_corr_file= basc.join_group_stability(indiv_stability_list, G_file, dataset_bootstraps, n_clusters, roi_mask_file)
    #k_mask,k_mask_file, icvs, cluster_voxel_scores,
    for i in range(0,len(subject_file_list)):
        icvs_file, cluster_voxel_scores_file, k_mask_file, ind_group_cluster_stability_file, individualized_group_clusters_img_file =basc.individual_group_clustered_maps(indiv_stability_list[i], clusters_G, roi_mask_file)


    return icvs_file, cluster_voxel_scores_file, k_mask_file, ind_group_cluster_stability_file, individualized_group_clusters_img_file #  G, clusters_G, cluster_voxel_scores, ism_gsm_corr, gsm_file, clusters_G_file, , ism_gsm_corr_file

#output_names=['icvs_file',
#                                               'cluster_voxel_scores_file',
#                                               'k_mask_file',
#                                               'ind_group_cluster_stability_file',
#                                               'individualized_group_clusters_img_file'],


def test_save_igcm_nifti(cluster_voxel_scores_file):
    
    import basc
    import utils
    subject_file_list= [home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub3/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz']

    roi_mask_file= home + '/git_repo/PyBASC/masks/LC_Quarter_Res.nii.gz'
    dataset_bootstraps=5
    timeseries_bootstraps=5
    n_clusters=3
    output_size=10
    bootstrap_list=list(range(0,dataset_bootstraps))
    cross_cluster=True
    roi2_mask_file= home + '/git_repo/PyBASC/masks/RC_Quarter_Res.nii.gz'
    cbb_block_size=None
    affinity_threshold= 0.5 #* len(subject_file_list)
    out_dir= home + '/PyBASC_outputs/SaveIGCMDebug2'
    run=True
    indiv_stability_list=[]
    for i in range(0,len(subject_file_list)):
        temp = basc.nifti_individual_stability(subject_file_list[i], roi_mask_file, timeseries_bootstraps, n_clusters, output_size, cross_cluster, roi2_mask_file, cbb_block_size, affinity_threshold)
        #temp=temp/timeseries_bootstraps
        indiv_stability_list.append(temp)
        
    G_file=[]
    for i in range(0,dataset_bootstraps):
        temp2= map_group_stability(indiv_stability_list, n_clusters, bootstrap_list)
        G_file.append(temp2)
        
    G, clusters_G, ism_gsm_corr, gsm_file, clusters_G_file, ism_gsm_corr_file= basc.join_group_stability(indiv_stability_list, G_file, dataset_bootstraps, n_clusters)
    #k_mask,k_mask_file, icvs, cluster_voxel_scores,
    for i in range(0,len(subject_file_list)):
        icvs_file, cluster_voxel_scores_file, k_mask_file, ind_group_cluster_stability_file =basc.individual_group_clustered_maps(indiv_stability_list[i], clusters_G, roi_mask_file)
 
    basc.save_igcm_nifti(cluster_voxel_scores_file,clusters_G_file,roi_mask_file)

#%% TEST BASC WORKFLOW



def test_basc_workflow_runner():

    from basc_workflow_runner import run_basc_workflow
    #import utils
    subject_file_list= [home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub3/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz',
                        home + '/git_repo/PyBASC/sample_data/sub2/Func_Quarter_Res.nii.gz']

    roi_mask_file= home + '/git_repo/PyBASC/masks/LC_Quarter_Res.nii.gz'
    dataset_bootstraps=20
    timeseries_bootstraps=20
    n_clusters=4
    output_size=10  
    blocklength=1
    bootstrap_list=list(range(0,dataset_bootstraps))
    cross_cluster=True
    similarity_metric='correlation'
    roi2_mask_file= home + '/git_repo/PyBASC/masks/RC_Quarter_Res.nii.gz'
    affinity_threshold= [0.0] * len(subject_file_list)
    out_dir= home + '/PyBASC_outputs/Testing_spatialconstraint'
    run=True
    
    

    basc_test= run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, output_size, bootstrap_list, proc_mem, similarity_metric, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, blocklength=blocklength, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)
    #PyBASC_test=run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, output_size, bootstrap_list, proc_mem, similarity_metric, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, blocklength=blocklength, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)


#%%
def heavy_test_basc_workflow_runner():
#%%
    from basc_workflow_runner import run_basc_workflow
    import utils
    subject_file_list=    ['/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060280/3mm_bandpassed_demeaned_filtered_antswarp.nii.gz',
                           '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060384/3mm_bandpassed_demeaned_filtered_antswarp.nii.gz']#,
#                           '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060429/3mm_bandpassed_demeaned_filtered_antswarp.nii.gz',
#                           '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060503/3mm_bandpassed_demeaned_filtered_antswarp.nii.gz',
#                           '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060603/3mm_bandpassed_demeaned_filtered_antswarp.nii.gz',
#                           '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060864/3mm_bandpassed_demeaned_filtered_antswarp.nii.gz']

    proc_mem= [3,6] #first is number of proc , second total number of mem

    roi_mask_file=home + '/git_repo/PyBASC/masks/Yeo7_3mmMasks/BilateralStriatumThalamus_3mm.nii.gz'
    dataset_bootstraps=2
    timeseries_bootstraps=100
    n_clusters=8
    output_size=500
    bootstrap_list=list(range(0,dataset_bootstraps))
    cross_cluster=True
    blocklength=2
    similarity_metric='correlation'
    roi2_mask_file=home + '/git_repo/PyBASC/masks/Yeo7_3mmMasks/YeoTest2.nii.gz'
    affinity_threshold= [0.0] * len(subject_file_list)
    out_dir= home + '/PyBASC_outputs/NewWOrkerTest'
    run=True
    

    basc_test= run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, output_size, bootstrap_list, proc_mem, similarity_metric, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, blocklength=blocklength, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)


#%%

def test_compare_stability_matrices():
    
    import utils
    import basc
    
    
    bootstrap=20
    blobs = generate_blobs()
    n_bootstraps=10
    n_clusters=5
    subjects=20
    
    
    ism_dataset = np.zeros((subjects, blobs.shape[0], blobs.shape[0]))
    ism_list = []
    for i in range(ism_dataset.shape[0]):
        ism_dataset[i] = utils.individual_stability_matrix(blobs + 0.2*np.random.randn(blobs.shape[0], blobs.shape[1]), n_bootstraps, n_clusters, affinity_threshold = 0.0)
        f = 'ism_dataset_%i.npy' % i
        ism_list.append(f)
        np.save(f, ism_dataset[i])

    indiv_stability_list=ism_list
   
    
    G = basc.map_group_stability(ism_list, n_bootstraps, n_clusters)
    
    gsm=np.load(G)
    gsm=gsm.astype("float64")

    corr= []
    corr= np.zeros((subjects,1))
    for j in range(ism_dataset.shape[0]):
        ism=ism_dataset[j].astype("float64")
        corr[j] = utils.compare_stability_matrices(gsm,ism)

    meandist5 = corr.mean()
    vardist5 = corr.var()
    sumdist5 = corr.cumsum()

#%%






def NED_heavy_basc_workflow_test():
#%%

    
    from basc_workflow_runner import run_basc_workflow
    import utils
    import time
    matrixtime = time.time()

#    subject_file_list= ['/data/rockland_sample/A00060603/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
#                         '/data/rockland_sample/A00060503/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
#                         '/data/rockland_sample/A00060429/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
#                         '/data/rockland_sample/A00060384/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
#                         '/data/rockland_sample/A00060280/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
#                         '/data/rockland_sample/A00059935/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
#                         '/data/rockland_sample/A00059875/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
#                         '/data/rockland_sample/A00059734/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
#                         '/data/rockland_sample/A00059733/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz']

    subject_file_list= ['/data/Projects/anikolai/rockland_downsampled/A00018030/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00027159/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00027167/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00027439/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00027443/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00030980/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00030981/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00031216/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00031219/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00031410/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00031411/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00031578/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00031881/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00032008/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00032817/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00033231/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00033714/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00034073/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00034074/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00034350/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00035291/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00035292/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00035364/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00035377/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00035869/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00035940/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00035941/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00035945/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00037125/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00037368/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00037458/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00037459/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00037483/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00038603/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00038706/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00039075/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00039159/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00039866/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00040342/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00040440/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00040556/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00040798/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00040800/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00040815/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00041503/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00043240/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00043282/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00043494/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00043740/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00043758/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00043788/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00044084/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00044171/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00050743/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00050847/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00051063/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00051603/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00051690/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00051691/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00051758/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00052069/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00052165/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00052183/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00052237/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00052461/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00052613/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00052614/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00053203/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00053320/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00053390/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00053490/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00053744/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00053873/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00054206/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00055693/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00056703/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00057405/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00057480/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00057725/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00057862/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00057863/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00057967/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058004/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058053/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058060/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058061/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058215/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058229/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058516/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058537/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058570/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058685/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00058951/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00059109/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00059325/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00059427/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00059733/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00059734/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00059865/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00059875/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00059935/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00060280/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00060384/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00060429/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00060503/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00060603/3mm_resampled.nii.gz',
    '/data/Projects/anikolai/rockland_downsampled/A00060846/3mm_resampled.nii.gz']

    roi_mask_file=home + '/git_repo/BASC/masks/BG_3mm.nii.gz'
   
    dataset_bootstraps=50
    timeseries_bootstraps=50
    n_clusters=3
    output_size=400
    cross_cluster=True
    bootstrap_list=list(range(0,dataset_bootstraps))
    proc_mem=[10,80]
    
    #roi2_mask_file=home + '/git_repo/BASC/masks/yeo_3mm.nii.gz'

    roi2_mask_file=home + '/git_repo/BASC/masks/yeo2_3mm.nii.gz'
    
    affinity_threshold= [0.5] * 107 #[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    out_dir=   '/data/Projects/anikolai/BASC_outputs/NKITest'
    run=True
    

    basc_test= run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, output_size, bootstrap_list, proc_mem, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)

    print((time.time() - matrixtime))
#%% TESTS UNDER CONSTRUCTION
# NDARRAY TO VOL
def test_ndarray_to_vol():
    import basc
    import nibabel as nb
    
    subject_file = home + '/git_repo/PyBASC/sample_data/sub1/Func_Quarter_Res.nii.gz'
    subject_file = home + '/git_repo/PyBASC/sample_data/test.nii.gz'
    data = nb.load(subject_file).get_data().astype('float32')
    roi_mask_file= home + '/git_repo/PyBASC/masks/LC_Quarter_Res.nii.gz'
    print( 'Data Loaded')

    
    roi_mask_file_nb = nb.load(roi_mask_file)

    roi_mask_nparray = nb.load(roi_mask_file).get_data().astype('float32').astype('bool')
    roi1data = data[roi_mask_nparray]
    
    data_array=roi1data
    sample_file=subject_file
    filename=home + '/git_repo/PyBASC/sample_data/ndarray_to_vol_test.nii.gz'
    
    basc.ndarray_to_vol(data_array, roi_mask_file, roi_mask_file, filename)


def test_co_clustering():
    
    import numpy as np
    import nibabel as nb
    from matplotlib import pyplot as plt
    import sklearn as sk
    from sklearn.datasets import make_biclusters
    from sklearn.datasets import samples_generator as sg
    from sklearn.cluster.bicluster import SpectralCoclustering
    from sklearn.metrics import consensus_score
    
    
    # REAL DATA
    subject_file= '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060280/3mm_bandpassed_demeaned_filtered_antswarp.nii.gz'
    roi_mask_file=home + '/git_repo/basc/masks/BG_3mm.nii.gz'
    roi2_mask_file=home + '/git_repo/basc/masks/yeo2_3mm.nii.gz'

    data = nb.load(subject_file).get_data().astype('float32')
    print( 'Data Loaded')

    print( 'Setting up NIS')
    roi_mask_file_nb = nb.load(roi_mask_file)
    roi2_mask_file_nb= nb.load(roi2_mask_file)

    roi_mask_nparray = nb.load(roi_mask_file).get_data().astype('float32').astype('bool')
    roi2_mask_nparray = nb.load(roi2_mask_file).get_data().astype('float32').astype('bool')


    roi1data = data[roi_mask_nparray]
    roi2data = data[roi2_mask_nparray]
    
    #add code that uploads the roi1data and roi2data, divides by the mean and standard deviation of the timeseries
    roi1data=sk.preprocessing.normalize(roi1data, norm='l2')
    roi2data=sk.preprocessing.normalize(roi2data, norm='l2')
    
    dist_btwn_data_1_2 = np.array(sp.spatial.distance.cdist(roi1data, roi2data, metric = 'correlation'))
    sim_btwn_data_1_2=1-dist_btwn_data_1_2
    sim_btwn_data_1_2[np.isnan(sim_btwn_data_1_2)]=0
    sim_btwn_data_1_2[sim_btwn_data_1_2<0]=0

    sim_btwn_data_1_2=sim_btwn_data_1_2+(np.random.rand(len(sim_btwn_data_1_2),len(sim_btwn_data_1_2[1,:])))/100
    sim_btwn_data_1_2[sim_btwn_data_1_2>1]=1
    
    sum(sum(sim_btwn_data_1_2==np.inf))
    sum(sum(sim_btwn_data_1_2==np.nan))


    model = SpectralCoclustering(n_clusters=5, random_state=0, n_init=100)
    model.fit(sim_btwn_data_1_2)
    
    fit_data = sim_btwn_data_1_2[np.argsort(model.row_labels_)]
    fit_data = fit_data[:, np.argsort(model.column_labels_)]
    
    plt.matshow(fit_data, cmap=plt.cm.Blues)
    plt.title("After biclustering; rearranged to show biclusters")
    
    plt.show()
    
    
    

    #SIMULATION DATA
    import numpy as np
    from matplotlib import pyplot as plt

    from sklearn.datasets import make_biclusters
    from sklearn.datasets import samples_generator as sg
    from sklearn.cluster.bicluster import SpectralCoclustering
    from sklearn.metrics import consensus_score

    #Creating Simulated Data
    data, rows, columns = make_biclusters(
        shape=(300, 100), n_clusters=5, noise=5,
        shuffle=False, random_state=0)
    
    plt.matshow(data, cmap=plt.cm.Blues)
    plt.title("Original dataset")
    
    data, row_idx, col_idx = sg._shuffle(data, random_state=0)
    plt.matshow(data, cmap=plt.cm.Blues)
    plt.title("Shuffled dataset")
    
    
    #Creating Model
    model = SpectralCoclustering(n_clusters=5, random_state=0)
    model.fit(data)
    score = consensus_score(model.biclusters_,
                            (rows[:, row_idx], columns[:, col_idx]))
    
    print("consensus score: {:.3f}".format(score))
    
    fit_data = data[np.argsort(model.row_labels_)]
    fit_data = fit_data[:, np.argsort(model.column_labels_)]
    
    plt.matshow(fit_data, cmap=plt.cm.Blues)
    plt.title("After biclustering; rearranged to show biclusters")
    
    plt.show()
    
    
    
    
    
    
    
    
    
    ####################################################################
    ####################################################################
    from sklearn import cluster
    import scipy as sp
    import time
    from sklearn import cluster, datasets
    import numpy as np
    from matplotlib import pyplot as plt

    from sklearn.datasets import make_biclusters
    from sklearn.datasets import samples_generator as sg
    from sklearn.cluster.bicluster import SpectralCoclustering
    from sklearn.metrics import consensus_score
    
    data1 = generate_simple_blobs(27)
    data2 = generate_simple_blobs(27)
    data2 = data2[0:150,:]
    

    print("Calculating Cross-clustering")
    print("Calculating pairwise distances between areas")
    
    dist_btwn_data_1_2 = np.array(sp.spatial.distance.cdist(roi1data, roi2data, metric = 'correlation'))
    sim_btwn_data_1_2=1-dist_btwn_data_1_2
    sim_btwn_data_1_2[sim_btwn_data_1_2<0]=0
    co_cluster=cluster.SpectralCoclustering()
    co_cluster.fit(sim_btwn_data_1_2)
    score = consensus_score(co_cluster.biclusters_,
                        (rows[:, row_idx], columns[:, col_idx]))

    print("consensus score: {:.3f}".format(score))

    fit_data = data[np.argsort(co_cluster.row_labels_)]
    fit_data = fit_data[:, np.argsort(co_cluster.column_labels_)]

    plt.matshow(fit_data, cmap=plt.cm.Blues)
    plt.title("After biclustering; rearranged to show biclusters")

    plt.show()
