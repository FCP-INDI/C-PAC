#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:35:02 2017

@author: aki.nikolaidis
"""
#adding test change
#import __init__

#__init__.import_all()

NKI_subject_file_list=[ '/data/rockland_sample/A00060603/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00060503/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00060429/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00060384/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00060280/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00059935/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00059875/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00059734/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
                         '/data/rockland_sample/A00059733/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz']

subject_file_list=NKI_subject_file_list

def NKI_Test2():
    #NKI TEST 2
    import os
    import numpy as np
    import nibabel as nb
    import utils
    import basc
    import pandas as pd
    import time
    
#    subject_file_list = ['/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060280/reduced50.nii.gz',
#                         '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060384/reduced50.nii.gz',
#                         '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060429/reduced50.nii.gz',
#                         '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060503/reduced50.nii.gz',
#                         '/Users/aki.nikolaidis/Desktop/NKI_SampleData/A00060603/reduced50.nii.gz',
#                            ]
#    
    
    #
    #                    ['/Users/aki.nikolaidis/BGDev_SampleData/A00060846/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                    '/Users/aki.nikolaidis/BGDev_SampleData/A00060603/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                    '/Users/aki.nikolaidis/BGDev_SampleData/A00060503/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                    '/Users/aki.nikolaidis/BGDev_SampleData/A00060429/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                    '/Users/aki.nikolaidis/BGDev_SampleData/A00060384/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                    '/Users/aki.nikolaidis/BGDev_SampleData/A00060280/bandpassed_demeaned_filtered_antswarp.nii.gz']
    #                     ['/data/rockland_sample/A00060603/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                     '/data/rockland_sample/A00060503/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                     '/data/rockland_sample/A00060429/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                     '/data/rockland_sample/A00060384/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                     '/data/rockland_sample/A00060280/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                     '/data/rockland_sample/A00059935/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                     '/data/rockland_sample/A00059875/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                     '/data/rockland_sample/A00059734/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
    #                     '/data/rockland_sample/A00059733/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz']
    
    
    
    
    
    roi_mask_file='/home/anikolai/C-PAC/CPAC/basc/sampledata/masks/BG.nii.gz'
    roi2_mask_file='/home/anikolai/C-PAC/CPAC/basc/sampledata/masks/yeo_2.nii.gz'
    
    #roi_mask_file= home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/Striatum_2thirdsRes.nii.gz'
    
    dataset_bootstraps=100
    timeseries_bootstraps=100
    n_clusters=4
    cross_cluster=True
    #roi2_mask_file= home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/Yeo_LowRes/yeo_2_2thirdsRes_bin.nii.gz'
    output_size = 2000
    
    ism_list=[]
    ism_file = np.zeros((len(subject_file_list), output_size, output_size))
    
    cbb_block_size=None
    affinity_threshold=0.5
    n_bootstraps=timeseries_bootstraps
    
    
    
        
           
    start = time.time()
    for i in range(len(subject_file_list)):
        data = nb.load(subject_file_list[int(i)]).get_data().astype('float32')
        print 'Data Loaded'
    
       
        print 'Setting up NIS'
        roi_mask_file_nb = nb.load(roi_mask_file)
        roi2_mask_file_nb= nb.load(roi2_mask_file)
    
        roi_mask_nparray = nb.load(roi_mask_file).get_data().astype('float32').astype('bool')
        roi2_mask_nparray = nb.load(roi2_mask_file).get_data().astype('float32').astype('bool')
    
    
    
        roi1data = data[roi_mask_nparray]
        roi2data = data[roi2_mask_nparray]
        print 'compressing data'
        data_dict1 = utils.data_compression(roi1data.T, roi_mask_file_nb, roi_mask_nparray, output_size)
        Y1_compressed = data_dict1['data']
        Y1_compressed = Y1_compressed.T
        Y1_labels = pd.DataFrame(data_dict1['labels'])
        Y1_labels=np.array(Y1_labels)
        print 'Y1 compressed'
        #index=pd.DataFrame(np.arange(1,Y1_labels.shape[0]+1))
        
        print 'compressing y2'
    
        data_dict2 = utils.data_compression(roi2data.T, roi2_mask_file_nb, roi2_mask_nparray, output_size)
        Y2_compressed = data_dict2['data']
        Y2_compressed=Y2_compressed.T
        #Y2_labels = pd.DataFrame(data_dict2['labels'])
        print 'y2 compressed'
        
        print('going into ism')
        print 'Calculating individual stability matrix of:', subject_file_list[int(i)]
        ism = utils.individual_stability_matrix(Y1_compressed, n_bootstraps, n_clusters, Y2_compressed, cross_cluster, cbb_block_size, affinity_threshold)
        
        print('expanding ism')
       # voxel_num=roi1data.shape[0]
        voxel_ism = utils.expand_ism(ism, Y1_labels)
        
    
        #def individual_stability_matrix(Y1, n_bootstraps, k_clusters, Y2=None, cross_cluster=False, cbb_block_size = None, affinity_threshold = 0.5):
        f = '/home/anikolai/C-PAC/CPAC/basc/tests/output/ism_dataset_%i.npy' % i
        ism_list.append(f)
        #np.save(f, ism_dataset[i])
        #ism_file[int(i)] = os.path.join(os.getcwd(), 'individual_stability_matrix.npy')
        np.save(ism_list[i], voxel_ism)
        print 'Saving individual stability matrix %s for %s' % (voxel_ism, subject_file_list[int(i)])
    
    print((time.time() - start))
    GroupAnalysisStart = time.time()
    G, clusters_G, cluster_voxel_scores, gsm_file, clusters_G_file, cluster_voxel_scores_file = basc.group_stability_matrix(ism_list, dataset_bootstraps, n_clusters=n_clusters)
    print((time.time() - GroupAnalysisStart))
    return {'G':G, 'clusters_G':clusters_G, 'cluster_voxel_scores':cluster_voxel_scores, 'gsm_file':gsm_file, 'clusters_G_file':clusters_G_file, 'cluster_voxel_scores_file':cluster_voxel_scores_file}



subject_file_list = ['/data/rockland_sample/A00018030/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00027159/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00027167/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00027439/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00027443/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00030980/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00030981/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00031216/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00031219/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00031410/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00031411/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00031578/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00031881/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00032008/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00032817/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00033231/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00033714/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00034073/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00034074/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00034350/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00035291/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00035292/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00035364/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00035377/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00035869/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00035940/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00035941/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00035945/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00037125/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00037368/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00037458/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00037459/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00037483/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00038603/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00038706/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00039075/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00039159/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00039866/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00040342/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00040440/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00040556/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00040798/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00040800/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00040815/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00041503/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00043240/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00043282/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00043494/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00043740/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00043758/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00043788/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00044084/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00044171/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00050743/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00050847/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00051063/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00051603/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00051690/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00051691/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00051758/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00052069/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00052165/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00052183/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00052237/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00052461/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00052613/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00052614/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00053203/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00053320/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00053390/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00053490/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00053744/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00053873/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00054206/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00055693/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00056703/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00057405/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00057480/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00057725/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00057862/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00057863/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00057967/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058004/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058053/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058060/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058061/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058215/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058229/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058516/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058537/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058570/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058685/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00058951/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00059109/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00059325/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00059427/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00059733/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00059734/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00059865/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00059875/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00059935/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00060280/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00060384/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00060429/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00060503/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00060603/functional_mni/_scan_clg_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz',
'/data/rockland_sample/A00060846/functional_mni/_scan_dsc_2_rest_645/bandpassed_demeaned_filtered_antswarp.nii.gz']