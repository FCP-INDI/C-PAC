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


def gsm_nifti(roi_mask_file, n_clusters, out_dir):
    import utils
    import basc
    import numpy as np
    import os
#Individual subject ISM to NIFTI and individual

#Inputs Subject ISM, ROIFile, 


    
    
    #for i in range(nSubjects):
    gsmdir=out_dir + '/workflow_output/basc_workflow_runner/basc/join_group_stability/'
    os.chdir(gsmdir)

    gsm=np.load(gsmdir + '/group_stability_matrix.npy')
    clusters_gsm = utils.cluster_timeseries(gsm, n_clusters, similarity_metric = 'correlation', affinity_threshold=0.0)
    clusters_gsm = clusters_gsm+1
    #niftifilename = gsmdir  +'/gsm_clust.nii.gz'
    #clusters_gsm_file = gsmdir +'/clusters_gsm.npy'
    #Saving Individual Level Cluster Solution
#    ndarray_to_vol(clusters_gsm, roi_mask_file, roi_mask_file, niftifilename)
#    np.save(clusters_gsm_file, clusters_gsm)
    
    
    cluster_ids = np.unique(clusters_gsm)
    nClusters = cluster_ids.shape[0]
    nVoxels = clusters_gsm.shape[0]
    gsm_cluster_voxel_scores = np.zeros((nClusters, nVoxels))
    k_mask=np.zeros((nVoxels, nVoxels))
    gsm_cluster_voxel_scores[:,:], k_mask[:,:] = utils.cluster_matrix_average(gsm, clusters_gsm)
    gsm_cluster_voxel_scores=gsm_cluster_voxel_scores.astype("uint8")
    
    grp_cluster_stability=[]
    grp_cluster_INSTABILITY=[]
    grp_cluster_stability_Diff=[]
    
    grp_cluster_stability_file = os.path.join(os.getcwd(), 'grp_cluster_stability.npy')
    grp_cluster_INSTABILITY_file = os.path.join(os.getcwd(), 'grp_cluster_INSTABILITY.npy')
    grp_cluster_stability_Diff_file = os.path.join(os.getcwd(), 'grp_cluster_stability_Diff.npy')
    gsm_cluster_voxel_scores_file = os.path.join(os.getcwd(), 'gsm_cluster_voxel_scores.npy')
    
    for k in cluster_ids:
        grp_cluster_stability.append(gsm_cluster_voxel_scores[(k-1),clusters_gsm==k].mean())
        grp_cluster_INSTABILITY.append(gsm_cluster_voxel_scores[(k-1),clusters_gsm!=k].mean())
        A, B = basc.ndarray_to_vol(gsm_cluster_voxel_scores[k-1,:], roi_mask_file, roi_mask_file, 'gsm_single_cluster%i_stability.nii.gz' % k)
    grp_cluster_stability=np.asarray(grp_cluster_stability)
    grp_cluster_INSTABILITY=np.asarray(grp_cluster_INSTABILITY)
    grp_cluster_stability_Diff=grp_cluster_stability-grp_cluster_INSTABILITY
    
    np.save(grp_cluster_stability_file, grp_cluster_stability)
    np.save(grp_cluster_INSTABILITY_file, grp_cluster_INSTABILITY)
    np.save(grp_cluster_stability_Diff_file, grp_cluster_stability_Diff)
    np.save(gsm_cluster_voxel_scores_file, gsm_cluster_voxel_scores)


    return