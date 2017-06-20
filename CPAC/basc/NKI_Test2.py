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

def main():
    print 'hey stranger'

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
    
    dataset_bootstraps=5
    timeseries_bootstraps=10
    n_clusters=4
    cross_cluster=True
    #roi2_mask_file= home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/Yeo_LowRes/yeo_2_2thirdsRes_bin.nii.gz'
    output_size = 1000
    
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



if __name__ == '__main__':
    main()
