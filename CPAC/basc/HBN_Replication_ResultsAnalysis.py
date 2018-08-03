#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 17:02:14 2018

@author: aki.nikolaidis
"""
import os
import numpy as np
from sklearn.metrics import adjusted_rand_score
import pandas as pd
import nibabel as nb

#TrueMotorClust=np.zeros((4332,4332))
#TrueMotorClust[0:2994,0:2994]=1
#TrueMotorClust[2994:4332,2994:4332]=1
#
#TrueMotorLabels=np.zeros((4332,1))
#TrueMotorLabels[0:2994]=1
#TrueMotorLabels[2994:4332]=2

data_dir='/Users/aki.nikolaidis/PyBASC_outputs'
GS=[1,10,30,60,100]
SNR=['0.05corr_2noise', '0.3corr_3noise','0.3corr_2noise']
Volumes=[100,200,400]

TrueClust=nb.load('/Users/aki.nikolaidis/git_repo/PyBASC/masks/TrueBGClust.nii.gz').get_data().astype('float32').astype('bool')
roi_mask_file='/Users/aki.nikolaidis/git_repo/PyBASC/masks/BG_3mm_TrueClust.nii.gz'
roi_mask_nparray = nb.load(roi_mask_file).get_data().astype('float32').astype('bool')
TrueBG = roi_mask_nparray[roi_mask_nparray]





TrueBG_GSM=adjacency_matrix(TrueBG)
TrueBG_GSM=TrueBG_GSM*1

import pdb; pdb.set_trace()

#SimResults=pd.DataFrame(columns=['AvgAcc', 'Reg1Acc', 'Reg2Acc', 'numsub', 'numvox', 'TRs', 'n_clusters', 'corrstrength', 'bootstraps', 'noiselevel', 'SNR'])
SimResults=pd.DataFrame(columns=['GBS', 'IBS', 'group_label_acc', 'gsm_acc', 'ism_gsm_corrmean', 'ism_gsm_corrstd'])

#'GBS', 'IBS', 'group_label_acc', 'gsm_acc', 'ism_gsm_corrmean', 'ism_gsm_corrstd'])

#GBSTest_0.3corr_2noise_100Vol1GS

for vols in Volumes:
    for SNR_level in SNR: 
        for bootstraps in GS:
            outdir= data_dir+'/GBSTest_'+SNR_level+'_' + str(vols)+ 'Vol' + str(bootstraps) + 'GS/'
            print(outdir)
            import pdb;pdb.set_trace()
            subdirs_all = [x[1] for x in os.walk(outdir)]                                                                            
            subdirs=subdirs_all[0]
            #out_dir= '/Users/aki.nikolaidis/PyBASC_outputs/WWW_BootstrapTest_100GS/dim_' + str(output_size) + '_' + str(similarity_metric) + '_' + str(n_clusters) + '_clusters_' +str(timeseries_bootstraps) +'_IndBS_' + str(blocklength) + '_block' + similarity_metric
            for subdir in subdirs:
                newdir=outdir + subdir
                os.chdir(newdir)
                #Go into each dir and calculate a bunch of things and put them all into a csv file for plotting later.
                print(newdir)
                path=os.path.normpath(newdir)
                specifics=path.split(os.sep)[5]
                dimreduction= specifics.split('_')[1]
                clusternum=specifics.split('_')[3]
                IBS=specifics.split('_')[5]
                GBS=bootstraps
               # import pdb; pdb.set_trace()
        
                group_cluster_labels=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/clusters_G.npy')
                group_label_acc=adjusted_rand_score(TrueBG.ravel(), group_cluster_labels)
                #import pdb;pdb.set_trace()
                gsm=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/group_stability_matrix.npy')
                gsm_acc= np.corrcoef(gsm.ravel(),TrueBG_GSM.ravel())[0][1]
                #import pdb;pdb.set_trace()
                ism_gsm_corr=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/ism_gsm_corr.npy')
                ism_gsm_corrmean= ism_gsm_corr.mean()
                ism_gsm_corrstd= ism_gsm_corr.std()
                
                #import pdb; pdb.set_trace()
                newdata=pd.DataFrame([[GBS, IBS, group_label_acc, gsm_acc, ism_gsm_corrmean, ism_gsm_corrstd]], columns=['GBS', 'IBS', 'group_label_acc', 'gsm_acc', 'ism_gsm_corrmean', 'ism_gsm_corrstd'])
        
                frames=[SimResults, newdata]
                SimResults= pd.concat(frames)
                            
SimResults.to_csv(data_dir+'/NewGBS_IBS4_29_2018.csv')
                            
