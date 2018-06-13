#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 09:04:14 2017

@author: aki.nikolaidis
"""

#    * Sample 100 samples from a normal distribution, with 1 feature
#    * Compute the mean
#    * Compute the mean by bootstrap procedure
#    * True mean is 0
#    * Compute the error in both situations.
#* Run this simulation 1000 times, if bootstrap is better the distribution on average will be greater than error of Sample mean.
#* Sample 2-
#    * Do the same, but sample with 10% are from a second normal distribution with variance of 100. And non 0 mean. Repeat entire procedure and look for difference between their errors.
#%%
def simbootstrap(data):
    import numpy as np
    k = int(np.ceil(float(data.shape[0])/1))
    r_ind = np.floor(np.random.rand(1,k)*data.shape[0])
    bootstrapsample=data[r_ind.astype('int')]
    return bootstrapsample



#%% FIRST APPROXiMATION- SIMULATION OF ONE DIMENSION
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
resamples=100 #how many times we sample the population
n=100 # how many data points we have in our sample
bootstraps=100 #how many bootstraps samples are collected on this sample

#adding 'noisy' data
intercept=50
variance=2
error_rarity=n/25

#initializing
bootstrapdata=[];
samplemeans=[];
finalbootstrapmeans=[];


for num in range(resamples):
    bootstrapdata=[]
    data=np.random.randn(n) #create data
    ######
    #CREATE NOISE DATA
    m=n/error_rarity #length of 'noisy' data
    morevariancedata= variance * np.random.randn(int(m)) + intercept #create noisy data
    alldata=data.tolist() + morevariancedata.tolist() #add noisy data
    data=np.asarray(alldata)
    ######
    samplemeans.append(data.mean()) #create sample mean and append to sample mean list
    
    for nums in range(bootstraps):
        tempsample= simbootstrap(data).mean() #bootstrap sample
        bootstrapdata.append(tempsample) #append bootstrap to list
    
    bootstrapdata=np.asarray(bootstrapdata).ravel()
    finalbootstrapmeans.append(bootstrapdata.mean()) #take mean of all bootstraps for this sample
        


finalbootstrapmeans=np.asarray(finalbootstrapmeans).ravel()
samplemeans=np.asarray(samplemeans)
error=samplemeans-finalbootstrapmeans

true_sample_error=0-samplemeans
true_bs_error=0-finalbootstrapmeans

#import pdb; pdb.set_trace()

#import matplotlib.pyplot as plt
#plt.hist(finalbootstrapmeans) 
#plt.hist(samplemeans) 


#print(correlation)
print('Average bootstrapmean is ', sum(finalbootstrapmeans)/len(finalbootstrapmeans))
print('Average samplemean is ', sum(samplemeans)/len(samplemeans))
#
print('Bootstrapmean Variance is ', finalbootstrapmeans.var())
print('Samplemean Variance is ', samplemeans.var())



sns.kdeplot(samplemeans, label="samplemeans")
sns.kdeplot(finalbootstrapmeans, label="bootstrapmeans")

#sns.kdeplot(error, label="error")
#sns.kdeplot(true_sample_error, label="true_sample_error")
#sns.kdeplot(true_bs_error, label="true_bs_error")

#ns.kdeplot(samplemeans, bw=.2, label="bw: 0.2")
#sns.kdeplot(samplemeans, bw=2, label="bw: 2")
plt.legend();



#%% MINI SIM OF  SEPARATE SIGNAL SOURCES
import numpy as np
import scipy as sp
import seaborn as sns
import pandas as pd
from sklearn.datasets import make_blobs

n=100

one=np.random.randn(n,1)
two=np.random.randn(n,1)


A=0.4 + 0.8*one
A=A+np.random.rand(n,1) +np.random.rand(n,1)


B=0.4 + 0.8*one + 0.8*two
B=B+np.random.rand(n,1) +np.random.rand(n,1)

C=0.4 + 0.8*two
C=C+np.random.rand(n,1) +np.random.rand(n,1)

A=pd.DataFrame(A, columns=['A'])
B=pd.DataFrame(B, columns=['B'])
C=pd.DataFrame(C, columns=['C'])
one=pd.DataFrame(one, columns=['one'])
two=pd.DataFrame(two, columns=['two'])

new=pd.concat([A,B,C,one,two], axis=1)

#sns.set(style="ticks")
sns.pairplot(new, kind="reg", markers="+", diag_kind="kde")

#1-sp.spatial.distance.cdist(A.T,one.T, metric='correlation')
#plt.scatter(A,two)
#np.corr(one[:,],two[:,])


#%% Null Generator
import numpy as np
import matplotlib.pyplot as plt

truth1500_1=np.zeros((1500,1500))
truth1500_1[0:1000,0:1000]=1
truth1500_1[1000:1500,1000:1500]=1

nullacclist=[]
nullrange=1000

for permutations in range(nullrange):
    null=np.random.rand(2250000)
    nullacc=np.corrcoef(truth1500_1.ravel(),null)[0,1]
    nullacclist.append(nullacc)

plt.hist(nullacclist)

#%%LargerSim

import numpy as np
import scipy as sp
import seaborn as sns
import pandas as pd
#import __init__
from utils import timeseries_bootstrap, adjacency_matrix
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn import cluster, datasets, preprocessing
from sklearn.cluster import FeatureAgglomeration
import seaborn as sns; sns.set(style="ticks", color_codes=True)


Reg1_True=np.load('/Users/aki.nikolaidis/git_repo/PyBASC/Reg1_True.npy')
Reg2_True=np.load('/Users/aki.nikolaidis/git_repo/PyBASC/Reg2_True.npy')


replications=100
numsub_list=[1]
numvox_list=[500]#[6171]
n_list=[200]
n_clusters_list=[2]
corrstrength_list=[0.5]
bootstraps_list=[1,10]#, 100, 200, 400, 600, 800]
noiselevel_list=[2.5]#, 1.5, 2, 2.5]#, 1.4, 1.8, 2.2, 2.6, 3.0]

truth1500_1=np.zeros((1500,1500))
truth1500_1[0:1000,0:1000]=1
truth1500_1[1000:1500,1000:1500]=1

truth1500_2=np.zeros((1500,1500))
truth1500_2[0:500,0:500]=1
truth1500_2[500:1500,500:1500]=1



#create matrix that has all accuracy and parameter information in it. Use Seaborne to plot the effects of different parameters on accuracy of group clustering.
SimResults=pd.DataFrame(columns=['AvgAcc', 'Reg1Acc', 'Reg2Acc', 'numsub', 'numvox', 'TRs', 'n_clusters', 'corrstrength', 'bootstraps', 'noiselevel', 'SNR', 'reps'])




for reps in range(replications):
    for noiselevel in noiselevel_list:
        for bootstraps in bootstraps_list:
            for corrstrength in corrstrength_list:
                for n_clusters in n_clusters_list:
                    for n in n_list:
                        for numvox in numvox_list:
                            for numsub in numsub_list:
                                GSM1=np.zeros((numvox*3, numvox*3))
                                GSM2=np.zeros((numvox*3, numvox*3))
                                pred1_stack=[];
                                pred2_stack=[];
                                for subs in range(numsub):
                                    #CREATING DATA
                                    one=np.random.randn(n,1)
                                    two=np.random.randn(n,1)
                                    region_one=[];
                                    region_two=[];
                                    region_A=[];
                                    region_B=[];
                                    region_C=[];
                                    for vox in range(numvox):
                                        print(vox)
                                        A=0.0 + corrstrength*one +  (corrstrength/10)*two
                                        A=A+(noiselevel*(np.random.rand(n,1)))
                                        A=preprocessing.normalize(A, norm='max', axis=0)
                                        
                                        B=0.0 + (corrstrength/2)*one + (corrstrength/2)*two
                                        B=B+ (noiselevel*(np.random.rand(n,1))) #Could add more noise to each
                                        B=preprocessing.normalize(B, norm='max', axis=0)
                                        
                                        C=0.0 + (corrstrength/10)*one + corrstrength*two
                                        C=C+(noiselevel*(np.random.rand(n,1)))
                                        C=preprocessing.normalize(C, norm='max', axis=0)
                                        
                                        region_one.append(one+np.random.rand(n,1))
                                        region_two.append(two+np.random.rand(n,1))
                                        region_A.append(A+np.random.rand(n,1))
                                        region_B.append(B+np.random.rand(n,1))
                                        region_C.append(C+np.random.rand(n,1))
                                        
                                        
                                           
                                    region_one=np.asarray(region_one)
                                    region_two=np.asarray(region_two)
                                    region_A=np.asarray(region_A)
                                    region_B=np.asarray(region_B)
                                    region_C=np.asarray(region_C)
                                    region_one=np.reshape(region_one, (numvox,n)).T
                                    region_two=np.reshape(region_two, (numvox,n)).T
                                    region_A=np.reshape(region_A, (numvox,n)).T
                                    region_B=np.reshape(region_B, (numvox,n)).T
                                    region_C=np.reshape(region_C, (numvox,n)).T
                                    
                                    
                                    
                                    Regions=np.block([region_A,region_B,region_C])
                                    
                                    N1 = Regions.shape[0]
                                    V1 = Regions.shape[1]
                                    
                                    reg1_ism = np.zeros((V1, V1))
                                    reg2_ism = np.zeros((V1, V1))
                                    
                                    for i in range(bootstraps):
                                        print('bootstraps num is', i)
                                        if bootstraps==1:
                                            bootstraps=0
                                            Regions_b1=Regions
                                            region_one_b1=region_one
                                            region_two_b1=region_two
                                        else:
                                            
                                            print('N1',N1)
                                            print('V1',V1)
                                            print(int(np.sqrt(N1)))
                                            block_size = int(np.sqrt(N1))
                                            #Regions_b1, block_mask = timeseries_bootstrap(Regions, block_size)
                                            tseries=Regions
                                            k = int(np.ceil(float(tseries.shape[0])/block_size))
                                            r_ind = np.floor(np.random.rand(1,k)*tseries.shape[0])
                                            blocks = np.dot(np.arange(0,block_size)[:,np.newaxis], np.ones([1,k]))
                                            #print('timeseries3')
                                        
                                            block_offsets = np.dot(np.ones([block_size,1]), r_ind)
                                            block_mask = (blocks + block_offsets).flatten('F')[:tseries.shape[0]]
                                            block_mask = np.mod(block_mask, tseries.shape[0])
                                            
                                            Regions_b1 =    Regions[block_mask.astype('int'), :]
                                            region_one_b1 = region_one[block_mask.astype('int'), :]
                                            region_two_b1 = region_two[block_mask.astype('int'), :]  
                                        
                                        print('Creating Similarity Data between voxels in ABC and region One')
                                        Dist_to_One= sp.spatial.distance.cdist(Regions_b1.T,region_one_b1.T, metric='correlation') # Try using 1- and euclidean? (Critical point here- correlation first then euclidean distance!!))
                                        Dist_of_ABC_One = sp.spatial.distance.pdist(Dist_to_One, metric = 'euclidean')
                                        Dist_mat_of_ABC_One = sp.spatial.distance.squareform(Dist_of_ABC_One)
                                        Sim_mat_of_ABC_One= 1- preprocessing.normalize(Dist_mat_of_ABC_One, norm='max') # Try max, then try standardizing Dist_of_ABC_One
                                        Sim_mat_of_ABC_One[Sim_mat_of_ABC_One<0]=0
                                        Sim_mat_of_ABC_One[Sim_mat_of_ABC_One>1]=1   
                                        #import pdb;pdb.set_trace()
                                        
                                        print('Creating Similarity Data between voxels in ABC and region Two')
                                        Dist_to_Two=sp.spatial.distance.cdist(Regions_b1.T,region_two_b1.T, metric='correlation') # Try using 1- and euclidean?
                                        Dist_of_ABC_Two =  sp.spatial.distance.pdist(Dist_to_Two, metric = 'euclidean')
                                        Dist_mat_of_ABC_Two = sp.spatial.distance.squareform(Dist_of_ABC_Two)
                                        Sim_mat_of_ABC_Two=1-preprocessing.normalize(Dist_mat_of_ABC_Two, norm='max')
                                        Sim_mat_of_ABC_Two[Sim_mat_of_ABC_Two<0]=0
                                        Sim_mat_of_ABC_Two[Sim_mat_of_ABC_Two>1]=1
                                        #import pdb; pdb.set_trace()
                                        
                                        print('Clustering ABC according to connectivity to region 1')
                                        
                                        ward1 = FeatureAgglomeration(n_clusters=n_clusters, affinity='euclidean', linkage='ward')    
                                        ward1.fit(Sim_mat_of_ABC_One)
                                        y_pred1 = ward1.labels_.astype(np.int)
                                        
    #                                    spectral1 = cluster.SpectralClustering(n_clusters, eigen_solver='arpack', random_state = 5, affinity="precomputed", assign_labels='discretize')
    #                                    spectral2 = cluster.SpectralClustering(n_clusters, eigen_solver='arpack', random_state = 5, affinity="precomputed", assign_labels='discretize')
    #                                    
    #                                    
    #                                    spectral1.fit(Sim_mat_of_ABC_One)
    #                                    y_pred1 = spectral1.labels_.astype(np.int)
                                        
                                        print('Clustering ABC according to connectivity to region 2')
                                        ward2 = FeatureAgglomeration(n_clusters=n_clusters, affinity='euclidean', linkage='ward')    
                                        ward2.fit(Sim_mat_of_ABC_Two)
                                        y_pred2 = ward2.labels_.astype(np.int)
    #                                    spectral2.fit(Sim_mat_of_ABC_Two)
    #                                    y_pred2 = spectral2.labels_.astype(np.int)
                                    
                                        y1_adj=adjacency_matrix(y_pred1)
                                        y2_adj=adjacency_matrix(y_pred2)
                                        
                                        reg1_ism += y1_adj
                                        reg2_ism += y2_adj
                                    #reg1_ism=np.asarray(reg1_ism)
                                    #reg2_ism=np.asarray(reg2_ism)
                                    
                                    #import pdb; pdb.set_trace()
                                    
                                    if bootstraps==0:
                                        #print('No Bootstraps!')
                                        reg1_final=reg1_ism
                                        reg2_final=reg2_ism
                                        reg1_pred=y_pred1
                                        reg2_pred=y_pred2
                                        
                                    else:
                                        #import pdb;pdb.set_trace()
                                    
                                        reg1_ism/=bootstraps
                                        reg2_ism/=bootstraps
                                        
    #                                    spectral1.fit(reg1_ism)
    #                                    reg1_pred = spectral1.labels_.astype(np.int)
                                        ward1.fit(reg1_ism)
                                        reg1_pred = ward1.labels_.astype(np.int)
                                        reg1_final=adjacency_matrix(reg1_pred)
                                        
                                        
    #                                    spectral2.fit(reg2_ism)
    #                                    reg2_pred = spectral2.labels_.astype(np.int)
                                        ward2.fit(reg2_ism)
                                        reg2_pred = ward2.labels_.astype(np.int)
                                        reg2_final=adjacency_matrix(reg2_pred)
                                    
                                    GSM1+=reg1_final
                                    GSM2+=reg2_final
                                    pred1_stack.append(reg1_pred)
                                    pred2_stack.append(reg2_pred)
                                    
                                    print('All Done')
                                       
    #                            Reg1Acc=np.corrcoef(GSM1.ravel(),Reg1_True.ravel())
    #                            Reg1Acc=Reg1Acc[0,1]
    #                            
    #                            Reg2Acc=np.corrcoef(GSM2.ravel(),Reg2_True.ravel())
    #                            Reg2Acc=Reg2Acc[0,1]
                                #Reg2Acc=1-sp.spatial.distance.cdist((GSM2.ravel(),Reg2_True.ravel()), metric='correlation')
                                SNR=corrstrength/noiselevel
                                
                                Reg1Acc=np.corrcoef(truth1500_1.ravel(),reg1_ism.ravel())[0,1]
                                Reg2Acc=np.corrcoef(truth1500_2.ravel(),reg2_ism.ravel())[0,1]
                                
                                AvgAcc=(Reg1Acc+Reg2Acc)/2
                                
                                newdata=pd.DataFrame([[AvgAcc, Reg1Acc, Reg2Acc, numsub, numvox, n, n_clusters, corrstrength, bootstraps, noiselevel, SNR, reps]], columns=['AvgAcc', 'Reg1Acc', 'Reg2Acc', 'numsub', 'numvox', 'TRs', 'n_clusters', 'corrstrength', 'bootstraps', 'noiselevel', 'SNR', 'reps'])
                                #newdata=pd.DataFrame([[numsub, numvox, n, n_clusters, corrstrength, bootstraps, noiselevel]], columns=['numsub', 'numvox', 'TRs', 'n_clusters', 'corrstrength', 'bootstraps', 'noiselevel'])
    
                                frames=[SimResults, newdata]
                                SimResults= pd.concat(frames)
                                SimResults.to_csv('/Users/aki.nikolaidis/git_repo/PyBASC/SimResultsIBS/TestingNoIBS.csv')
                                
    #                            
#plotframe= pd.DataFrame([[SimResults['corrstrength'], SimResults['bootstraps'], SimResults['Reg1Acc'], SimResults['Reg2Acc']]], columns=['corrstrength','bootstraps','Reg1Acc','Reg2Acc'])                            
##plotdata=pd.DataFrame(plotframe)
#plot = sns.pairplot(SimResults)

#plt.imshow(GSM1)
#plt.imshow(GSM2)
#plt.imshow(reg1_ism)
#plt.imshow(reg2_ism)
#
#
#Reg1Acc=np.corrcoef(truth1500_1.ravel(),reg1_ism.ravel())
#Reg2Acc=np.corrcoef(truth1500_2.ravel(),reg2_ism.ravel())
#
#np.corrcoef(truth1500_1.ravel(),GSM1.ravel())
#np.corrcoef(truth1500_2.ravel(),GSM2.ravel())

#FILL COLUMNS OF THE MATRIX WITH THESE VALUES
#olumns: [bootstraps, Signal/Noise ratio, accuracy1, accuracy2]

#plt.imshow(Reg1_True)
#plt.imshow(Reg2_True)
#
#pred1_stack=np.asarray(pred1_stack)
#pred2_stack=np.asarray(pred2_stack)



#for noiselevel in noiselevel_list:
#    for bootstraps in bootstraps_list:
#        for corrstrength in corrstrength_list:
#            for n_clusters in n_clusters_list:
#                for n in n_list:
#                    for numvox in numvox_list:
#                        for numsub in numsub_list
#                            for subs in range(numsub):


#print('Region One Correlation is ', Reg1Acc)
#print('Region Two Correlation is ', Reg2Acc)


#%% JACCARD INDEX PERMUTATION TESTING
from sklearn.metrics import adjusted_rand_score, jaccard_similarity_score, adjusted_mutual_info_score

permutations=100

np.random.permutation(pred1_stack)
permscore=[];
permscore2=[];
pred_fullstack=np.vstack((pred1_stack,pred2_stack))
pred_fullstack_perm=np.random.permutation(pred_fullstack)
test_stat=adjusted_rand_score(pred1_stack[0,:],pred2_stack[0,:])
test_stat2=adjusted_mutual_info_score(pred1_stack[0,:],pred2_stack[0,:])
#plt.imshow(pred_fullstack_perm)

for perm in range(permutations):
    pred_fullstack_perm=np.random.permutation(pred_fullstack)
    
    permscore.append(adjusted_rand_score(pred_fullstack_perm[0,:],pred_fullstack_perm[1,:]))
    #permscore2.append(adjusted_mutual_info_score(pred_fullstack_perm[0,:], pred_fullstack_perm[1,:]))
    
permscore=np.asarray(permscore)
plt.hist(permscore)
print(test_stat)


permscore2=np.asarray(permscore2)
plt.hist(permscore2)
print(test_stat2)

#One way to do this would be to 
#1 - Take the group level GSM and calculate difference to other.
#2 - For perm in permutations:
#        Permute over ISMs
#        Calculate new GSMs
#        Test Jaccard Distance between groups
#        Save.
#3-  Test difference with original.

#permutation testing

#%% CREATE DATA FOR SIMULATION ON REAL DATA

import numpy as np
import scipy as sp
import seaborn as sns
import pandas as pd
import basc
import utils
import os
#import __init__
from utils import timeseries_bootstrap, adjacency_matrix
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn import cluster, datasets, preprocessing
import seaborn as sns; sns.set(style="ticks", color_codes=True)


#Reg1_True=np.load('/Users/aki.nikolaidis/git_repo/PyBASC/Reg1_True.npy')
#Reg2_True=np.load('/Users/aki.nikolaidis/git_repo/PyBASC/Reg2_True.npy')


#NOTES
rightvox = 1284 #Region A voxels, should be mostly motor signal 
thalvox = 1775# Region B voxels- should be a compostiion of signal from both= like region B 
leftvox = 1273# Region C voxels, should be mostly visual signal. 
fullvox= 4332


TrueMotorClust=np.zeros((4332,4332))
TrueMotorClust[0:3059,0:3059]=1
TrueMotorClust[3059:4332,3059:4332]=1

TrueVisualClust=np.zeros((4332,4332))
TrueVisualClust[0:1284,0:1284]=1
TrueVisualClust[1284:4332,1284:4332]=1

TrueBGClust=np.zeros((4332,4332))
TrueBGClust[0:1284,0:1284]=1
TrueBGClust[1284:3059,1284:3059]=1
TrueBGClust[3059:4332,3059:4332]=1

#truth1500_1[0:1000,0:1000]=1
#truth1500_1[1000:1500,1000:1500]=1


numsub=10
motorvox=11656
visualvox=11941

n_list=[200]
n_clusters_list=[2]
corrstrength_list=[0.05,0.3]
bootstraps_list=[1]
noiselevel_list=[2, 3]

for subs in range(numsub):
    print(subs)

for n in n_list:
    for noiselevel in noiselevel_list:
        for corrstrength in corrstrength_list:
            for subs in range(numsub):
                #CREATING DATA
                
    
                SimBG_right=[];
                SimBG_thal=[];
                SimBG_left=[];
                SimMotor=[];
                SimVisual=[];
    
                SimMotor_signal=np.random.randn(n,1)
                SimVisual_signal=np.random.randn(n,1)
                
                interaction=SimMotor_signal*SimVisual_signal # THink about adding this to the model
    
    #            for vox in range(BGvox):
    #                print('BG', vox)
    #                
    #                if (vox < 2000):
    #                    A=0.0 + corrstrength*SimMotor_signal + (corrstrength/10)*SimVisual_signal + (0.1*interaction)
    #                    A=A+(noiselevel*(np.random.rand(n,1)))
    #                    SimBG.append(A+np.random.rand(n,1))
    #                if (2000 <= vox < 4000):
    #                    B=0.0 + (corrstrength/2)*SimMotor_signal + (corrstrength/2)*SimVisual_signal + interaction
    #                    B=B+ (noiselevel*(np.random.rand(n,1))) #Could add more noise to each
    #                    SimBG.append(B+np.random.rand(n,1))
    #                if (vox>=4000):
    #                    C=0.0 + (corrstrength/10)*SimMotor_signal + corrstrength*SimVisual_signal + (0.1 * interaction)
    #                    C=C+(noiselevel*(np.random.rand(n,1)))
    #                    SimBG.append(C+np.random.rand(n,1))
                
                for vox in range(rightvox):
                    A=0.0 + corrstrength*SimMotor_signal# + (corrstrength/10)*SimVisual_signal + (0.1*interaction)
                    A=A+(noiselevel*(np.random.rand(n,1)))
                    A=preprocessing.normalize(A, norm='max', axis=0)
                    SimBG_right.append(A)
            
                for vox in range(thalvox):
                    B=0.0 + (corrstrength)*SimMotor_signal + (corrstrength)*SimVisual_signal# + 0.5*interaction
                    B=B+ (noiselevel*(np.random.rand(n,1))) #Could add more noise to each
                    B=preprocessing.normalize(B, norm='max', axis=0)
                    SimBG_thal.append(B)
                        
                for vox in range(leftvox):
                    C=0.0 + corrstrength*SimVisual_signal #+ (corrstrength/10)*SimMotor_signal  + (0.1*interaction)
                    C=C+(noiselevel*(np.random.rand(n,1)))
                    C=preprocessing.normalize(C, norm='max', axis=0)
                    SimBG_left.append(C)            
                
                for vox in range(motorvox):
                    print('Motor', vox)
                    Motor=0.0 + SimMotor_signal + (noiselevel*np.random.rand(n,1)) #+ (0.1*SimVisual_signal)
                    Motor=preprocessing.normalize(Motor, norm='max', axis=0)
                    SimMotor.append(Motor)
                    
                for vox in range(visualvox):
                    print('Visual', vox)
                    Visual= 0.0 + SimVisual_signal + (noiselevel*np.random.rand(n,1)) #+ (0.1*SimMotor_signal)
                    Visual=preprocessing.normalize(Visual, norm='max', axis=0)
                    SimVisual.append(Visual)
                
                SimBG_right=np.asarray(SimBG_right)
                SimBG_thal=np.asarray(SimBG_thal)
                SimBG_left=np.asarray(SimBG_left)
                SimMotor=np.asarray(SimMotor)
                SimVisual=np.asarray(SimVisual)
                #region_A=np.reshape(region_A, (numvox,n)).T
                #SimBG=np.reshape(SimBG,(BGvox,n))
                #import pdb; pdb.set_trace()
                SimBG_right=np.reshape(SimBG_right, (rightvox,n))
                SimBG_thal=np.reshape(SimBG_thal, (thalvox,n))
                SimBG_left=np.reshape(SimBG_left, (leftvox,n))
                
                SimMotor=np.reshape(SimMotor,(motorvox,n))
                SimVisual=np.reshape(SimVisual,(visualvox,n))

                data_array_right = SimBG_right
                roi_mask_file_right = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Right_Caud_Put_Pall_bin_thresh_3mm.nii.gz'
                sample_file_right = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Right_Caud_Put_Pall_bin_thresh_3mm.nii.gz'
                filename_right = '/Users/aki.nikolaidis/git_repo/PyBASC/SimData4/SimBG_right_' + 'sub_'+ str(subs) +'corrstrength_' + str(corrstrength) + 'noise_' + str(noiselevel) + '.nii.gz'
                
                data_array_thal = SimBG_thal
                roi_mask_file_thal = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Bilateral_Thalamus3mm.nii.gz'
                sample_file_thal = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Bilateral_Thalamus3mm.nii.gz'
                filename_thal = '/Users/aki.nikolaidis/git_repo/PyBASC/SimData4/SimBG_Thal_' + 'sub_'+ str(subs) +'corrstrength_' + str(corrstrength) + 'noise_' + str(noiselevel) + '.nii.gz'
                
                data_array_left = SimBG_left
                roi_mask_file_left = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Left_Caud_Put_Pall_bin_thresh_3mm.nii.gz'
                sample_file_left = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Left_Caud_Put_Pall_bin_thresh_3mm.nii.gz'
                filename_left = '/Users/aki.nikolaidis/git_repo/PyBASC/SimData4/SimBG_left_' + 'sub_'+ str(subs) +'corrstrength_' + str(corrstrength) + 'noise_' + str(noiselevel) + '.nii.gz'
                
                #roi_mask_file='/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/BilateralStriatumThalamus_3mm.nii.gz'
                #roi2_mask_file='/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/Yeo_All_7_3mm.nii.gz'
    
                data_array2 = SimVisual
                roi_mask_file2 = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/Yeo_1_3mm.nii.gz'
                sample_file2 = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/Yeo_1_3mm.nii.gz'
                filename2 = '/Users/aki.nikolaidis/git_repo/PyBASC/SimData4/SimVisual_' + 'sub_'+ str(subs) +'corrstrength_' + str(corrstrength) + 'noise_' + str(noiselevel) + '.nii.gz'
                
                data_array3 = SimMotor
                roi_mask_file3 = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/Yeo_2_3mm.nii.gz'
                sample_file3 = '/Users/aki.nikolaidis/git_repo/PyBASC/masks/Yeo7_3mmMasks/Yeo_2_3mm.nii.gz'
                filename3 = '/Users/aki.nikolaidis/git_repo/PyBASC/SimData4/SimMotor_' + 'sub_'+ str(subs) +'corrstrength_' + str(corrstrength) + 'noise_' + str(noiselevel) + '.nii.gz'
                
              
                
                #write Regions to nifti file
                basc.ndarray_to_vol(data_array_right, roi_mask_file_right, sample_file_right, filename_right)
                basc.ndarray_to_vol(data_array_thal, roi_mask_file_thal, sample_file_thal, filename_thal)
                basc.ndarray_to_vol(data_array_left, roi_mask_file_left, sample_file_left, filename_left)
                #write Region One to nifti file
                basc.ndarray_to_vol(data_array2, roi_mask_file2, sample_file2, filename2)
                #
                basc.ndarray_to_vol(data_array3, roi_mask_file3, sample_file3, filename3)
                
                niftiadditionfile='fslmaths ' + filename_right + ' -add ' + filename_thal + ' -add ' + filename_left + ' -add ' + filename2 + ' -add ' + filename3 + ' /Users/aki.nikolaidis/git_repo/PyBASC/SimData4/sub_' + str(subs) +'corr_' + str(corrstrength) + '_noise_' + str(noiselevel) + '_TRs_' + str(n) + '.nii.gz'
    
                os.system(niftiadditionfile)
                
                os.system('rm /Users/aki.nikolaidis/git_repo/PyBASC/SimData4/Sim*')
                
                

#%%    
#plt.imshow(Sim_to_One)
#plt.imshow(Sim_of_ABC_One)
#plt.imshow(Sim_mat_of_ABC_One)
#plt.imshow(reg1_final)

Reg1Corr=1-sp.spatial.distance.pdist((GSM1.ravel(),Reg1_True.ravel()), metric='correlation')
Reg2Corr=1-sp.spatial.distance.pdist((GSM2.ravel(),Reg2_True.ravel()), metric='correlation')

print('Region One Correlation is ', Reg1Corr)
print('Region Two Correlation is ', Reg2Corr)
#%%


plt.imshow(reg1_final)
#%%
plt.imshow(reg2_final)

#%%
plt.imshow(GSM1)

#%%
plt.imshow(GSM2)

1-sp.spatial.distance.pdist((reg2_final.ravel(),Reg2_True.ravel()), metric='correlation')
#
#
#A=1-sp.spatial.distance.cdist(region_A,region_one, metric='correlation')

#%%

region_one=np.asarray(region_one)
region_onedf=pd.DataFrame(region_one)

A=0.4 + 0.8*one
A=A+np.random.rand(n,1) +np.random.rand(n,1)


B=0.4 + 0.8*one + 0.8*two
B=B+np.random.rand(n,1) +np.random.rand(n,1)

C=0.4 + 0.8*two
C=C+np.random.rand(n,1) +np.random.rand(n,1)

A=pd.DataFrame(A, columns=['A'])
B=pd.DataFrame(B, columns=['B'])
C=pd.DataFrame(C, columns=['C'])
one=pd.DataFrame(one, columns=['one'])
two=pd.DataFrame(two, columns=['two'])

new=pd.concat([A,B,C,one,two], axis=1)

#sns.set(style="ticks")
sns.pairplot(new, kind="reg", markers="+", diag_kind="kde")




#%% Plotting GMMs SciKitLearn

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from sklearn import mixture

n_samples = 300

# generate random sample, two components
np.random.seed(0)

# generate spherical data centered on (20, 20)
shifted_gaussian = np.random.randn(n_samples, 2) + np.array([20, 20])

# generate zero centered stretched Gaussian data
#C = np.array([[0, -0.7], [3.5, .7]])
#stretched_gaussian = np.dot(np.random.randn(n_samples, 2), C)

C = np.array([[-0.5, -4],
              [0.5, 4]])
stretched_gaussian = np.dot(np.random.randn(n_samples, 2), C)

C2 = np.array([[0.6, -2],
              [-0.6, 2]])
stretched_gaussian2 = np.dot(np.random.randn(n_samples, 2), C2)

C3 = np.array([[5, 2],
              [9, 2]])
stretched_gaussian3 = np.dot(np.random.randn(n_samples, 2), C3)

# concatenate the two datasets into the final training set
#X_train = np.vstack([shifted_gaussian, stretched_gaussian])
X_train = np.vstack([stretched_gaussian, stretched_gaussian2, stretched_gaussian3])


# fit a Gaussian Mixture Model with two components
clf = mixture.GaussianMixture(n_components=2, covariance_type='full')
clf.fit(X_train)

# display predicted scores by the model as a contour plot
x = np.linspace(-20., 30.)
y = np.linspace(-20., 40.)
X, Y = np.meshgrid(x, y)
XX = np.array([X.ravel(), Y.ravel()]).T
Z = -clf.score_samples(XX)
Z = Z.reshape(X.shape)

CS = plt.contour(X, Y, Z, norm=LogNorm(vmin=1.0, vmax=1000.0),
                 levels=np.logspace(0, 3, 10))
CB = plt.colorbar(CS, shrink=0.8, extend='both')
plt.scatter(X_train[:, 0], X_train[:, 1], .8)

plt.title('Negative log-likelihood predicted by a GMM')
plt.axis('tight')
plt.show()

#%%More GMM Sample Code



#%%JOVO Matlab Code
#clear
#
#n=100;
#
#T=1000;
#
#B=1000;
#
# 
#
#xbar=zeros(T,1);
#
#xbbar=zeros(T,1);
#
#for i=1:T
#
#    x=randn(n,1);
#
#    xbar(i)=mean(x);
#
#    
#
#    xboot=zeros(B,1);
#
#    for b=1:B
#
#        idx=ceil(n*rand(n*0.8,1));
#
#        xboot(b)=mean(x(idx));
#        
#        k = int(np.ceil(float(data.shape[0])/1))
#        r_ind = np.floor(np.random.rand(1,k)*data.shape[0])
#        bootstrapsample=data[r_ind.astype('int')]
#
#    end
#
#    xbbar(i)=mean(xboot);
#
#end
#
# 
#
#%%
#
# 
#
#clc
#
#err(i)=xbar(i)-xbbar(i);
#
#[mean(xbar), mean(xbbar), mean(xbar-xbbar)]
#
#[std(xbar),std(xbbar)]


#%%
#    #print('Calculating Timeseries Bootstrap')
#    k = int(np.ceil(float(tseries.shape[0])/block_size))
#
#    r_ind = np.floor(np.random.rand(1,k)*tseries.shape[0])
#    blocks = np.dot(np.arange(0,block_size)[:,np.newaxis], np.ones([1,k]))
#    #print('timeseries3')
#
#    block_offsets = np.dot(np.ones([block_size,1]), r_ind)
#    block_mask = (blocks + block_offsets).flatten('F')[:tseries.shape[0]]
#    block_mask = np.mod(block_mask, tseries.shape[0])
#    #print('block_offsets shape0', block_offsets.shape[0])
#    #print('block_offsets shape1', block_offsets.shape[1])
#    #print('block_mask shape0', block_mask.shape[0])
#    #print('block_mask shape1', block_mask.shape[1])
#    #print('tseries shape0', tseries[block_mask.astype('int'), :].shape[0])
#    #print('tseries shape1', tseries[block_mask.astype('int'), :].shape[1])
#    #print('Finished: ', (time.time() - bootstraptime), ' seconds')
#    return tseries[block_mask.astype('int'), :]

##SimTest
#
#import numpy as np
#import pandas as pd
#n = 100
#df = pd.DataFrame()
#np.random.seed(1)
#df['x1'] = np.random.randn(n)
#df['x2'] = np.random.randn(n)
#df['x3'] = np.random.randn(n)
#df['x4'] = np.random.randn(n)
#df['y'] = 10 + -100*df['x1'] +  75*df['x3'] + np.random.randn(n)
#
#import statsmodels.formula.api as smf
#results = smf.ols('y ~ x1 + x2 + x3 + x4 ', data=df).fit()
#print(results.summary())
#
#
##* Create simulation data-
#
#
#import numpy as np
#
#sample=np.random.randn(n)
#print('sample mean is', sample.mean())
#data=np.random.randn(n)
#simbootstrap(data,bootstraps)
