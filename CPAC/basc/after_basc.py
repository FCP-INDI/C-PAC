#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 14:48:45 2018

Inputs: List of All ISM

@author: aki.nikolaidis
"""


#package imports
import os
from sklearn.cluster.bicluster import SpectralCoclustering
import sklearn as sk
from sklearn import cluster, datasets, preprocessing
import scipy as sp
import time 
from sklearn.cluster import FeatureAgglomeration
from sklearn.feature_extraction import image
from sklearn.metrics import adjusted_rand_score

all_ind_clusterlabels=[]
out_dir='/Users/aki.nikolaidis/PyBASC_outputs/Self_Sim_WWS/dim_500_correlation_2_clusters_100_IndBS_1_blockcorrelation'
ismdir=out_dir + '/workflow_output/basc_workflow_runner/basc/individual_stability_matrices/mapflow/'
os.chdir(ismdir)
subdirs_all = [x[1] for x in os.walk(ismdir)]                                                                            
subdirs=subdirs_all[0]
n_clusters=2
    #roi_mask_nparray = nb.load(roi_mask_file).get_data().astype('float32').astype('bool')
ward = FeatureAgglomeration(n_clusters=n_clusters, affinity='euclidean', linkage='ward')    


for subdir in subdirs:
    os.chdir(ismdir + subdir)
    temp_ism_file = os.path.join(os.getcwd(), 'individual_stability_matrix.npy')
    temp_ism=np.load(temp_ism_file)
    print("Calculating Hierarchical Clustering")
    ward.fit(temp_ism)
    y_pred = ward.labels_.astype(np.int)
    #APPLY CLUSTERING
    all_ind_clusterlabels.append(y_pred)
    
    

#for ism in ismlist:
    #temp=np.load(ism)
    #all_ind_clusterlabels.append(temp)
    

score_similarities= np.zeros((len(all_ind_clusterlabels),len(all_ind_clusterlabels)))

for i in range(len(all_ind_clusterlabels)):
    for j in range(len(all_ind_clusterlabels)):
        score_similarities[i,j]=adjusted_rand_score(all_ind_clusterlabels[i],all_ind_clusterlabels[j])

ypred2=[]
for i in range(1,5):
    print(i)
    ward2 = FeatureAgglomeration(n_clusters=i, affinity='euclidean', linkage='ward')    
    ward2.fit(score_similarities)
    #ytemp=
    ypred2.append(ward2.labels_.astype(np.int))   
    
    
    
    
#CODE TO COMPARE Cluster ASSIGNMENTS ACROSS NETWORKS
#LOAD ALL CLUSTER ASSIGNMENTS INTO MATRIX WITH LABEL ON EACH COLUMN OF RESPECTIVE SOURCE
#CALCULATE SIMILARITIES ACROSS THEM WITH RAND INDEX
#CLUSTER THE RESULTING SIMILARITY MATRIX
#PLOT THE HEAT MAP AND THE DENDROGRAM CLUSTERING    
    

import os
import numpy as np
import utils
import basc
gsm=np.load('/Users/aki.nikolaidis/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/YeoNetworksGSMs/SelfGSM.npy')
clusters_G=np.load('/Users/aki.nikolaidis/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/YeoNetworksGSMs/Selfclusters_G.npy')

nClusters=5
nVoxels=6171
cluster_voxel_scoresSelf = np.zeros((nClusters, nVoxels))
k_mask=np.zeros((nVoxels, nVoxels))
cluster_voxel_scoresSelf[:,:], k_mask[:,:] = utils.cluster_matrix_average(gsm, clusters_G)
    
    
#Plotly Code
import plotly
import plotly.plotly as py
from plotly.tools import FigureFactory as FF
plotly.tools.set_credentials_file(username='AkiNikolaidis', api_key='WnYhFXtIHn2FHv1gbTLH')

import numpy as np

X = np.random.rand(15, 15)
dendro = FF.create_dendrogram(X)
dendro['layout'].update({'width':800, 'height':500})
py.plot(dendro, filename='simple_dendrogram')

import plotly.plotly as py
from plotly.tools import FigureFactory as FF

import numpy as np

X = np.random.rand(10, 10)
names = ['Jack', 'Oxana', 'John', 'Chelsea', 'Mark', 'Alice', 'Charlie', 'Rob', 'Lisa', 'Lily']
fig = FF.create_dendrogram(X, orientation='left', labels=names)
fig['layout'].update({'width':800, 'height':800})
A=py.iplot(fig, filename='dendrogram_with_labels')






#PLOTTING PLOTLY DOUBLE DENDROGRAM OF ALL 35 CLUSTER STABILITY

plotly.tools.set_credentials_file(username='AkiNikolaidis', api_key='WnYhFXtIHn2FHv1gbTLH')
import plotly.plotly as py
from plotly.graph_objs import *
from plotly.tools import FigureFactory as FF

import numpy as np
from scipy.spatial.distance import pdist, squareform



# get data
#data = np.genfromtxt("http://files.figshare.com/2133304/ExpRawData_E_TABM_84_A_AFFY_44.tab",
#                     names=True,usecols=tuple(range(1,30)),dtype=float, delimiter="\t")
#data_array = data.view((np.float, len(data.dtype.names)))
#data_array = data_array.transpose()
#labels = data.dtype.names

data_array = Clust_Voxel_Scores_Total#.view((np.float, len(data.dtype.names)))
#data_array = cluster_stab.transpose()
labels=labels2=('Yeo1_1','Yeo1_2','Yeo1_3','Yeo1_4','Yeo1_5','Yeo2_1','Yeo2_2','Yeo2_3','Yeo2_4','Yeo2_5','Yeo3_1','Yeo3_2','Yeo3_3','Yeo3_4','Yeo3_5','Yeo4_1','Yeo4_2','Yeo4_3','Yeo4_4','Yeo4_5','Yeo5_1','Yeo5_2','Yeo5_3','Yeo5_4','Yeo5_5','Yeo6_1','Yeo6_2','Yeo6_3','Yeo6_4','Yeo6_5','Yeo7_1','Yeo7_2','Yeo7_3','Yeo7_4','Yeo7_5','Self_1','Self_2','Self_3','Self_4','Self_5')


# Initialize figure by creating upper dendrogram
figure = FF.create_dendrogram(data_array, orientation='bottom', labels=labels)
for i in range(len(figure['data'])):
    figure['data'][i]['yaxis'] = 'y2'

# Create Side Dendrogram
dendro_side = FF.create_dendrogram(data_array, orientation='right')
for i in range(len(dendro_side['data'])):
    dendro_side['data'][i]['xaxis'] = 'x2'

# Add Side Dendrogram Data to Figure
figure['data'].extend(dendro_side['data'])

# Create Heatmap
dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
dendro_leaves = list(map(int, dendro_leaves))
data_dist = pdist(data_array, metric='correlation')
heat_data = squareform(data_dist)
heat_data = heat_data[dendro_leaves,:]
heat_data = heat_data[:,dendro_leaves]

heatmap = Data([
    Heatmap(
        x = dendro_leaves, 
        y = dendro_leaves,
        z = heat_data,    
        colorscale = 'YIGnBu'
    )
])

heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

# Add Heatmap Data to Figure
figure['data'].extend(Data(heatmap))

# Edit Layout
figure['layout'].update({'width':800, 'height':800,
                         'showlegend':False, 'hovermode': 'closest',
                         })
# Edit xaxis
figure['layout']['xaxis'].update({'domain': [.15, 1],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'ticks':""})
# Edit xaxis2
figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                   'mirror': False,
                                   'showgrid': False,
                                   'showline': False,
                                   'zeroline': False,
                                   'showticklabels': False,
                                   'ticks':""}})

# Edit yaxis
figure['layout']['yaxis'].update({'domain': [0, .85],
                                  'mirror': False,
                                  'showgrid': False,
                                  'showline': False,
                                  'zeroline': False,
                                  'showticklabels': False,
                                  'ticks': ""})
# Edit yaxis2
figure['layout'].update({'yaxis2':{'domain':[.825, .975],
                                   'mirror': False,
                                   'showgrid': False,
                                   'showline': False,
                                   'zeroline': False,
                                   'showticklabels': False,
                                   'ticks':""}})

# Plot!
py.iplot(figure, filename='dendrogram_with_heatmap')

py.plot(figure)

dendro_side['layout']['xaxis']

help(FF.create_dendrogram)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#Code to Borrow    
#ISM NIFTI CODE

        
#
#    
##SKLEARN CODE
#import numpy as np
#from matplotlib import pyplot as plt
#
#from sklearn.datasets import make_biclusters
#from sklearn.datasets import samples_generator as sg
#from sklearn.cluster.bicluster import SpectralCoclustering
#from sklearn.metrics import consensus_score
#
#data, rows, columns = make_biclusters(
#    shape=(300, 300), n_clusters=5, noise=5,
#    shuffle=False, random_state=0)
#
#plt.matshow(data, cmap=plt.cm.Blues)
#plt.title("Original dataset")
#
#data, row_idx, col_idx = sg._shuffle(data, random_state=0)
#plt.matshow(data, cmap=plt.cm.Blues)
#plt.title("Shuffled dataset")
#
#model = SpectralCoclustering(n_clusters=5, random_state=0)
#model.fit(data)
#score = consensus_score(model.biclusters_,
#                        (rows[:, row_idx], columns[:, col_idx]))
#
#print("consensus score: {:.3f}".format(score))
#
#fit_data = data[np.argsort(model.row_labels_)]
#fit_data = fit_data[:, np.argsort(model.column_labels_)]
#
#plt.matshow(fit_data, cmap=plt.cm.Blues)
#plt.title("After biclustering; rearranged to show biclusters")
#
#plt.show()