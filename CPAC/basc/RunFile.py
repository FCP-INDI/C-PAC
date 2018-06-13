##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Created on Fri Jul 28 10:44:38 2017
#
#@author: aki.nikolaidis
#"""
#
#import BASC
#from BASC import *
import __init__
import utils
import os
import gc
import numpy as np
import scipy.stats
from os.path import expanduser
from basc_workflow_runner import run_basc_workflow
from basc import save_igcm_nifti, create_group_cluster_maps, ism_nifti, gsm_nifti

home = expanduser("~")
proc_mem= [3,16]


subject_file_list = ['/Users/aki.nikolaidis/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/NKI_Dev/645_Data/TimeTest/645_852_A00029979.nii.gz',
                     '/Users/aki.nikolaidis/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/NKI_Dev/645_Data/TimeTest/645_852_A00031217.nii.gz']


#
roi_mask_file='/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/BilateralStriatumThalamus_3mm.nii.gz'


roi2_mask_file='/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/Yeo_All_7_3mm.nii.gz'



dataset_bootstrap_list=[1,10,30,60,100]
timeseries_bootstrap_list=[1,10,30,60,100]
similarity_metric_list=['correlation'] #['correlation','euclidean','cityblock', 'cosine']
blocklength_list=[1]#[0.5,1,2]
n_clusters_list=[2,4,6,8,10,12,14,16]
output_sizes=[800]#[10,100,600,1200]

cross_cluster=False
affinity_threshold= [0.0] * len(subject_file_list)

ism_gsm_stability=[]
ind_clust_stab_mat=[]
ind_clust_stab_summary=[[1, 2, 3, 4, 5]]
run=True

for (dataset_bootstraps, timeseries_bootstraps) in zip(dataset_bootstrap_list,timeseries_bootstrap_list):
    
    bootstrap_list= list(np.ones(dataset_bootstraps, dtype=int)*dataset_bootstraps)

    for similarity_metric in similarity_metric_list:
        for blocklength in blocklength_list:
                for n_clusters in n_clusters_list:
                    for output_size in output_sizes:
                        #import pdb; pdb.set_trace()
                        out_dir= home + '/PyBASC_outputs/SelfClusterTest5/dim_' + str(output_size) + '_' + str(similarity_metric) + '_' + str(n_clusters) + '_clusters_' +str(timeseries_bootstraps) +'_IndBS_' + str(blocklength) + '_block' + similarity_metric
                        #out_dir= '/Users/aki.nikolaidis/PyBASC_outputs/Testing_Ward/dim_' + str(output_size) + '_' + str(similarity_metric) + '_' + str(n_clusters) + '_clusters_' +str(timeseries_bootstraps) +'_IndBS_' + str(blocklength) + '_block' + similarity_metric
                        PyBASC_test=run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, output_size, bootstrap_list, proc_mem, similarity_metric, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, blocklength=blocklength, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)
#                        del PyBASC_test
#                        gc.collect()
                        #import pdb; pdb.set_trace()
                        ism_gsm_stability.append(np.load(out_dir + '/workflow_output/ism_gsm_corr_file/ism_gsm_corr.npy'))
                        ind_clust_stab_mat = np.load(out_dir + '/workflow_output/ind_group_cluster_stability_set/ind_group_cluster_stability_set.npy')
                        ind_clust_stab_summary=np.concatenate((ind_clust_stab_summary, np.array([[n_clusters, output_size, ind_clust_stab_mat.mean(), scipy.stats.variation(ind_clust_stab_mat).mean(), (ind_clust_stab_mat.mean() - scipy.stats.variation(ind_clust_stab_mat).mean())]])))
                        
                        #Run Group ClusterMaps
                        gsm_file = out_dir + '/workflow_output/basc_workflow_runner/basc/join_group_stability/group_stability_matrix.npy'
                        clusters_G_file = out_dir + '/workflow_output/basc_workflow_runner/basc/join_group_stability/clusters_G.npy'
                        os.chdir(out_dir +'/workflow_output/basc_workflow_runner/basc/join_group_stability/')
                        create_group_cluster_maps(gsm_file,clusters_G_file,roi_mask_file)
                        #Run IGCM on all individual subjects
                        clustvoxscoredir=out_dir + '/workflow_output/basc_workflow_runner/basc/individual_group_clustered_maps/mapflow/'
                        clusters_G_file= out_dir + '/workflow_output/basc_workflow_runner/basc/join_group_stability/clusters_G.npy'
                        #clusters_G = np.load(clusters_G_file)
                        
                        #ism_nifti(roi_mask_file, n_clusters, out_dir)
                        #gsm_nifti(roi_mask_file, n_clusters, out_dir)
                
                
                
                #        subdirs_all = [x[1] for x in os.walk(clustvoxscoredir)]                                                                            
                #        subdirs=subdirs_all[0]
                #        for subdir in subdirs:
                #        
            #            #import pdb; pdb.set_trace()
            #            clustvoxscorefile=clustvoxscoredir + subdir+ '/cluster_voxel_scores.npy'
            #            #clustvoxscores=np.load(clustvoxscorefile)
            #            os.chdir(clustvoxscoredir + '/' + subdir)
            #            save_igcm_nifti(clustvoxscorefile,clusters_G_file, roi_mask_file)
            #        
                    
                    print('saving files: ism_gsm_stability')
                    ism_gsm_stability_file=os.path.join(out_dir, 'ism_gsm_stability_'+ str(n_clusters)+ '.npy')
                    np.save(ism_gsm_stability_file, ism_gsm_stability)
                    ism_gsm_stability=[]
            
                print('saving files: ind_clust_stab_summary')
                ind_clust_stab_summary_file=os.path.join(out_dir, 'ind_clust_stab_summary.npy')
                np.save(ind_clust_stab_summary_file, ind_clust_stab_summary)
    
    #PSEUDO CODE- FOR THE ANALYSIS WITH THE HIGHEST MEAN AND LOWEST CV
    # RUN POST ANALYSIS- INCLUDING FOR LOOPING OVER ALL SUBJECTS FOLDERS AND CREATING NIFTI FILES FOR EVERY CLUSTER
    
    
#os.system('say "Hello, your analysis has completed. Please take a look at the output files, have a nice day, and dont forget to tip your server"')    
        
#clustervoxscores=np.load(outdir,)
