#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 15:58:30 2018

@author: aki.nikolaidis
"""

import numpy as np
import pandas as pd
import os

Replications=['10Min','20Min','30Min','40Min']

ref_data='/data2/Projects/BASC/HNU_SSI/PyBASC_outputs/Yeo2_Ref_gsr1_scrub0'
HNUResults=pd.DataFrame(columns=['Replications','clusternum','GBS_IBS', 'group_label_acc', 'gsm_acc', 'ism_gsm_corrmean', 'ism_gsm_corrstd', 'ref_rep_ismgsmcorr'])


#/data2/Projects/BASC/HNU_SSI/PyBASC_outputs/Yeo2_Ref_gsr1_scrub0/dim_800_correlation_8_clusters_60_IndBS_1_blockcorrelation
#collect correlations between ISMs, and ISM/GSM relationships.
data_dir='/data2/Projects/BASC/HNU_SSI/PyBASC_outputs/'

#for Replication in Replications:
#    outdir= data_dir+'Yeo2_' +Replication+ '_gsr1_scrub0' 
#    print(outdir)
#    subdirs_all = [x[1] for x in os.walk(outdir)]                                                                            
#    subdirs=subdirs_all[0]
#    #out_dir= '/Users/aki.nikolaidis/PyBASC_outputs/WWW_BootstrapTest_100GS/dim_' + str(output_size) + '_' + str(similarity_metric) + '_' + str(n_clusters) + '_clusters_' +str(timeseries_bootstraps) +'_IndBS_' + str(blocklength) + '_block' + similarity_metric
#    for subdir in subdirs:
#        refdir=refdata+subdir
#        ref_group_cluster_label=np.load(refdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/clusters_G.npy')
#        ref_gsm=np.load(refdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/group_stability_matrix.npy')
#        ref_ismgsmcorr=np.load(refdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/ism_gsm_corr.npy')
#        
#        
#        newdir=outdir + subdir
#        os.chdir(newdir)
#        
#        #Go into each dir and calculate a bunch of things and put them all into a csv file for plotting later.
#        print(newdir)
#        #newdir='/data2/Projects/BASC/HNU_SSI/PyBASC_outputs/Yeo2_Ref_gsr1_scrub0/dim_800_correlation_8_clusters_60_IndBS_1_blockcorrelation'
#        path=os.path.normpath(newdir)
#        specifics=path.split(os.sep)[7]
#        dimreduction= specifics.split('_')[1]
#        clusternum=specifics.split('_')[3]
#        GBS_IBS=specifics.split('_')[5]
#        #GBS=bootstraps
#       # import pdb; pdb.set_trace()
#
#        group_cluster_labels=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/clusters_G.npy')
#        group_label_acc=adjusted_rand_score(ref_group_cluster_label, group_cluster_labels)
#        #import pdb;pdb.set_trace()
#        gsm=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/group_stability_matrix.npy')
#        gsm_acc= np.corrcoef(gsm.ravel(),ref_gsm.ravel())[0][1]
#        #import pdb;pdb.set_trace()
#        ism_gsm_corr=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/ism_gsm_corr.npy')
#        ism_gsm_corrmean= ism_gsm_corr.mean()
#        ism_gsm_corrstd= ism_gsm_corr.std()
#        ref_rep_ismgsmcorr=np.corrcoef(ism_gsm_corr, ref_ismgsmcorr)[0][1]
#        
#        #import pdb; pdb.set_trace()
#        newdata=pd.DataFrame([[Replications ,clusternum, GBS_IBS, group_label_acc, gsm_acc, ism_gsm_corrmean, ism_gsm_corrstd, ref_rep_ismgsmcorr]], columns=['Replications','clusternum','GBS_IBS', 'group_label_acc', 'gsm_acc', 'ism_gsm_corrmean', 'ism_gsm_corrstd', 'ref_rep_ismgsmcorr'])
#        frames=[HNUResults, newdata]
#        HNUResults= pd.concat(frames)
#
#HNUResults.to_csv(data_dir+'/HNU_5_18_2018.csv')



### INDIVIDUAL LEVEL RESULTS ANALYSIS

HNU_Ind_Results=pd.DataFrame(columns=['Replications','clusternum','GBS_IBS', 'SubjectNum', 'Ref_Rep_ISM_Corr', 'rep_ref_ism_gsm_corr_corr', 'sub_ref_ism_gsm_corr', 'sub_rep_ism_gsm_corr'])


#/data2/Projects/BASC/HNU_SSI/PyBASC_outputs/Yeo2_Ref_gsr1_scrub0/dim_800_correlation_8_clusters_60_IndBS_1_blockcorrelation
#collect correlations between ISMs, and ISM/GSM relationships.
data_dir='/data2/Projects/BASC/HNU_SSI/PyBASC_outputs/'
refdata='/data2/Projects/BASC/HNU_SSI/PyBASC_outputs/Yeo2_Ref_gsr1_scrub0/'

for Replication in Replications:
    outdir= data_dir+'Yeo2_' +Replication+ '_gsr1_scrub0/' 
    print(outdir)
    subdirs_all = [x[1] for x in os.walk(outdir)]                                                                            
    subdirs=subdirs_all[0]
    #out_dir= '/Users/aki.nikolaidis/PyBASC_outputs/WWW_BootstrapTest_100GS/dim_' + str(output_size) + '_' + str(similarity_metric) + '_' + str(n_clusters) + '_clusters_' +str(timeseries_bootstraps) +'_IndBS_' + str(blocklength) + '_block' + similarity_metric
    for subdir in subdirs:
        refdir=refdata+subdir
        #ref_group_cluster_label=np.load(refdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/clusters_G.npy')
        #ref_gsm=np.load(refdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/group_stability_matrix.npy')
        #ref_ismgsmcorr=np.load(refdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/ism_gsm_corr.npy')
        
        
        newdir=outdir + subdir
        os.chdir(newdir)
        
        #Go into each dir and calculate a bunch of things and put them all into a csv file for plotting later.
        print(newdir)
        #newdir='/data2/Projects/BASC/HNU_SSI/PyBASC_outputs/Yeo2_Ref_gsr1_scrub0/dim_800_correlation_8_clusters_60_IndBS_1_blockcorrelation'
        path=os.path.normpath(newdir)
        specifics=path.split(os.sep)[7]
        dimreduction= specifics.split('_')[1]
        clusternum=specifics.split('_')[3]
        GBS_IBS=specifics.split('_')[5]
        #GBS=bootstraps
        #import pdb; pdb.set_trace()
       
        rep_ism_gsm_corr=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/ism_gsm_corr.npy')
        ref_ism_gsm_corr=np.load(refdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/ism_gsm_corr.npy')
       
        rep_ref_ism_gsm_corr_corr= ism_corr= np.corrcoef(rep_ism_gsm_corr.ravel(),ref_ism_gsm_corr.ravel())[0][1]
       
        SubjectNum=0
        subjectpath = newdir+'/workflow_output/basc_workflow_runner/basc/individual_stability_matrices/mapflow/'
        refpath = refdir +'/workflow_output/basc_workflow_runner/basc/individual_stability_matrices/mapflow/'
       
        subjectdirs_all = [x[1] for x in os.walk(subjectpath)]                                                                            
        subjectdirs=subjectdirs_all[0]
        for subject in subjectdirs:
           #NEED TO FIGURE OUT WHETHER POSSIBLE TO USE A SINGLE SUBJECT LEVEL FOR LOOP TO GET THE REFERENCE ONES AS WELL?
           #LOAD REF DATA
           
           ref_ism= np.load(refpath+subject+'/individual_stability_matrix.npy')
           rep_ism= np.load(subjectpath+subject+'/individual_stability_matrix.npy')
           Ref_Rep_ISM_Corr= np.corrcoef(rep_ism.ravel(),ref_ism.ravel())[0][1]
           
           #rep_clust_labels=that
           #LOAD REP DATA
           # CORRELATED REP AND REF DATA
           
           sub_ref_ism_gsm_corr= ref_ism_gsm_corr[SubjectNum]
           sub_rep_ism_gsm_corr= rep_ism_gsm_corr[SubjectNum]
           
           SubjectNum=SubjectNum+1
       

#        group_cluster_labels=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/clusters_G.npy')
#        group_label_acc=adjusted_rand_score(ref_group_cluster_label, group_cluster_labels)
#        #import pdb;pdb.set_trace()
#        gsm=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/group_stability_matrix.npy')
#        gsm_acc= np.corrcoef(gsm.ravel(),ref_gsm.ravel())[0][1]
#        #import pdb;pdb.set_trace()
#        ism_gsm_corr=np.load(newdir+'/workflow_output/basc_workflow_runner/basc/join_group_stability/ism_gsm_corr.npy')
#        ism_gsm_corrmean= ism_gsm_corr.mean()
#        ism_gsm_corrstd= ism_gsm_corr.std()
#        ref_rep_ismgsmcorr=np.corrcoef(ism_gsm_corr, ref_ismgsmcorr)[0][1]
        
        #import pdb; pdb.set_trace()
        newdata=pd.DataFrame([[Replications,clusternum, GBS_IBS, SubjectNum, Ref_Rep_ISM_Corr, rep_ref_ism_gsm_corr_corr, sub_ref_ism_gsm_corr, sub_rep_ism_gsm_corr]], columns=['Replications','clusternum','GBS_IBS', 'SubjectNum', 'Ref_Rep_ISM_Corr', 'rep_ref_ism_gsm_corr_corr', 'sub_ref_ism_gsm_corr', 'sub_rep_ism_gsm_corr'])
        frames=[HNU_Ind_Results, newdata]
        HNU_Ind_Results= pd.concat(frames)

HNU_Ind_Results.to_csv(data_dir+'/HNU_Ind_5_31.csv')











#
#
#for bootstrap in bootstrap_list:
#    for sub in sublist:
#        ref=np.load(path to sub ism)
#        rep1=np.load(path to sub ism 1IBS)
#        rep10=np.load(path to sub ism 10IBS)
#        rep30=np.load(path to sub ism 30IBS)
#        rep60=np.load(path to sub ism 60IBS)
#        rep100=np.load(path to sub ism 100IBS)
#        
#        ref1corr=np.corrcoef(ref.ravel(),rep1.ravel())[0][1]
#        ref10corr=np.corrcoef(ref.ravel(),rep10.ravel())[0][1]
#        ref30corr=np.corrcoef(ref.ravel(),rep30.ravel())[0][1]
#        ref60corr=np.corrcoef(ref.ravel(),rep60.ravel())[0][1]
#        ref100corr=np.corrcoef(ref.ravel(),rep100.ravel())[0][1]
#
#
#                HBNResults=pd.DataFrame([[sub, bootstrap, ref1corr, ref10corr, ref30corr, ref60corr, ref100corr]], columns=['Subject', 'bootstrap', 'ref1corr', 'ref10corr', 'ref30corr', 'ref60corr', 'ref100corr'])
#

