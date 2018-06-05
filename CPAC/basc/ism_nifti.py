#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 18:12:07 2017

@author: aki.nikolaidis
"""
#roi_mask_file='/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/BilateralStriatumThalamus_3mm.nii.gz'
##ism=np.load('/Users/aki.nikolaidis/PyBASC_outputs/Group_CMTestingFullData/dim_400_4_clusters/workflow_output/basc_workflow_runner/basc/individual_stability_matrices/mapflow/_individual_stability_matrices0/individual_stability_matrix.npy')
#n_clusters=4
#output_size=400
#out_dir= '/Users/aki.nikolaidis/PyBASC_outputs/ISM_Testing_500/dim_' + str(output_size) + '_' + str(n_clusters) + '_clusters'

#Extra code to assess within cluster stability
#clusterism
#Clusterstability=[]
#for k in clustvoxscore:
#    clusterstability.append=


def ism_nifti(roi_mask_file, n_clusters, out_dir):
    import utils
    import basc
    import numpy as np
    import os
#Individual subject ISM to NIFTI and individual

#Inputs Subject ISM, ROIFile, 


    
    
    #for i in range(nSubjects):
    ismdir=out_dir + '/workflow_output/basc_workflow_runner/basc/individual_stability_matrices/mapflow/'
    os.chdir(ismdir)
    os.chdir(out_dir +'/workflow_output/basc_workflow_runner/basc/individual_stability_matrices/mapflow/')
    subdirs_all = [x[1] for x in os.walk(ismdir)]                                                                            
    subdirs=subdirs_all[0]
    
    for subdir in subdirs:
        os.chdir(ismdir + subdir)
        ism=np.load(ismdir + subdir + '/individual_stability_matrix.npy')
        clusters_ism = utils.cluster_timeseries(ism, n_clusters, similarity_metric = 'correlation', affinity_threshold=0.0)
        clusters_ism = clusters_ism+1
        niftifilename = ismdir + subdir +'/ism_clust.nii.gz'
        clusters_ism_file = ismdir + subdir +'/clusters_ism.npy'
        #Saving Individual Level Cluster Solution
        ndarray_to_vol(clusters_ism, roi_mask_file, roi_mask_file, niftifilename)
        np.save(clusters_ism_file, clusters_ism)
        
        
        cluster_ids = np.unique(clusters_ism)
        nClusters = cluster_ids.shape[0]
        nVoxels = clusters_ism.shape[0]
        ism_cluster_voxel_scores = np.zeros((nClusters, nVoxels))
        k_mask=np.zeros((nVoxels, nVoxels))
        ism_cluster_voxel_scores[:,:], k_mask[:,:] = utils.cluster_matrix_average(ism, clusters_ism)
        ism_cluster_voxel_scores=ism_cluster_voxel_scores.astype("uint8")
        
        ind_cluster_stability=[]
        ind_cluster_INSTABILITY=[]
        ind_cluster_stability_Diff=[]
        
        ind_cluster_stability_file = os.path.join(os.getcwd(), 'ind_cluster_stability.npy')
        ind_cluster_INSTABILITY_file = os.path.join(os.getcwd(), 'ind_cluster_INSTABILITY.npy')
        ind_cluster_stability_Diff_file = os.path.join(os.getcwd(), 'ind_cluster_stability_Diff.npy')
        ism_cluster_voxel_scores_file = os.path.join(os.getcwd(), 'ism_cluster_voxel_scores.npy')
        
        os.chdir(ismdir + '/' + subdir)
        
        for k in cluster_ids:
            ind_cluster_stability.append(ism_cluster_voxel_scores[(k-1),clusters_ism==k].mean())
            ind_cluster_INSTABILITY.append(ism_cluster_voxel_scores[(k-1),clusters_ism!=k].mean())
            A, B = basc.ndarray_to_vol(ism_cluster_voxel_scores[k-1,:], roi_mask_file, roi_mask_file, 'ism_single_cluster%i_stability.nii.gz' % k)
        ind_cluster_stability=np.asarray(ind_cluster_stability)
        ind_cluster_INSTABILITY=np.asarray(ind_cluster_INSTABILITY)
        ind_cluster_stability_Diff=ind_cluster_stability-ind_cluster_INSTABILITY
        
        np.save(ind_cluster_stability_file, ind_cluster_stability)
        np.save(ind_cluster_INSTABILITY_file, ind_cluster_INSTABILITY)
        np.save(ind_cluster_stability_Diff_file, ind_cluster_stability_Diff)
        np.save(ism_cluster_voxel_scores_file, ism_cluster_voxel_scores)


    return