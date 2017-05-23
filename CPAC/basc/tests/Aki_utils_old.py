#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 19:09:11 2017

@author: aki.nikolaidis
"""
#%%
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


matplotlib.style.use('ggplot')
home = expanduser("~")

#%%
#Data Preparation

#load in dataset/func file- these are exchangable. (make sure it's in MNI152 space)
#dataset = datasets.fetch_adhd(n_subjects=2)
func_file = home + '/Desktop/residual_antswarp.nii.gz'
func = nb.load(func_file)


#load masks
bg_file = home + '/Dropbox/1_Projects/1_Research/2_BASC/Data/BasalGanglia_MNI2mm/BG.nii.gz'
yeo_7_lib_file = home + '/Dropbox/1_Projects/1_Research/2_BASC/Data/Yeo_JNeurophysiol11_MNI152/FSL_Reorient2Std_Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz'

#put masks into func space
bg = resample_img(bg_file, target_affine=func.affine,target_shape=func.shape[:3], interpolation='nearest')
bg = bg.get_data()

yeo_7_lib = resample_img(yeo_7_lib_file, target_affine=func.affine,target_shape=func.shape[:3], interpolation='nearest')
yeo_7_lib= yeo_7_lib.get_data()

#%%
print('check your mask and func data have the same shape')
print('')
print('func shape:', func.shape)
print('')
print('mask shapes:')
print('bg shape', bg.shape)
print('yeo_7_lib', yeo_7_lib.shape)

#%%
#Get the vector size of the func data
new_shape = func.shape[0]*func.shape[1]*func.shape[2]
new_shape

#put func data into a numpy array
func = func.get_data()

#string out the vector of functional data
func = np.resize(func, (new_shape,func.shape[-1] ))

#string out the vectors of mask data
bg_mask = np.resize(bg, (new_shape,1))
yeo_7_lib = np.resize(yeo_7_lib, (new_shape,1))


#%%
#Extract masks from functional data
#Action- add remainder networks
yeo_motor = func[np.where( yeo_7_lib == 2 )[0], :]
yeo_vis = func[np.where( yeo_7_lib == 1 )[0], :]
bg_func = func[np.where( bg_mask == 1 )[0], :]

print 'motor', type(yeo_motor), yeo_motor.shape, 'vis', yeo_vis.shape, 'bg', bg_func.shape
#%%

#Testing out Sample data on sclustering algorithm
#Actions for today to figure out- use the scikit learn clustering code instead of the ncut
#Figure out how to compute distances between voxels of A based on connectivity to B and vice versa.
#       #This may be the cdist function

#Figure out how to create required distance matrices
#transform to similarity metrics- similarity = np.exp(-beta * distance / distance.std())
#calculate distances- Y = cdist(XA, XB, 'euclidean')
#reshape save the y_pred as an image to check how it looks
sampledata1=generate_blobs()
np.random.seed(30)
offset = np.random.randn(30)
x1 = np.random.randn(200,30) + 2*offset
x2 = np.random.randn(100,30) + 44*np.random.randn(30)
x3 = np.random.randn(400,30)
sampledata2 = np.vstack((x1,x2,x3))
sampledata3 = sampledata2[0:400,:]


#%%
bg_func_mini1 = bg_func[0:1000,:]
bg_func_mini2 = bg_func[1001:2000,:]
bg_func_mini3 = bg_func[1500:3400,:]

dist=sp.spatial.distance.cdist(bg_func_mini3, bg_func_mini1, 'euclidean')
dist_pd=pd.DataFrame(dist)

y_predict = cluster_timeseries(dist_pd, 3, similarity_metric = 'correlation')


#%%
#plt.imshow(dist_pd, cmap='hot', interpolation='nearest')

fig, ax = plt.subplots()

im = ax.imshow(dist_pd, cmap='hot', interpolation='nearest')
fig.colorbar(im)

ax.axis('tight')

plt.show()


#blobs = generate_blobs_3d()



#%% TESTING BASC.PY FUNCTIONS

subject_file= home + '/Desktop/residual_antswarp.nii.gz'
roi_mask_file= home + '/Dropbox/1_Projects/1_Research/2_BASC/Data/BasalGanglia_MNI2mm/Right_Caudate.nii.gz'
n_bootstraps=10
k_clusters=2
cross_cluster=True
roi2_mask_file= home + '/Dropbox/1_Projects/1_Research/2_BASC/Data/BasalGanglia_MNI2mm/Left_Caudate.nii.gz'
cbb_block_size=None
affinity_threshold=0.5

nifti_individual_stability(subject_file, roi_mask_file, n_bootstraps, k_clusters, cross_cluster=True, roi2_mask_file=roi2_mask_file, cbb_block_size=None, affinity_threshold=0.5)

plt.imshow(ism, cmap='hot', interpolation='nearest')



#%%
X= bg_func

    # normalize dataset for easier parameter selection
X = StandardScaler().fit_transform(X)

    # estimate bandwidth for mean shift
bandwidth = cluster.estimate_bandwidth(X, quantile=0.3)

    # connectivity matrix for structured Ward
connectivity = kneighbors_graph(X, n_neighbors=10, include_self=False)
    # make connectivity symmetric
connectivity = 0.5 * (connectivity + connectivity.T)
spectral = cluster.SpectralClustering(n_clusters=2,eigen_solver='arpack', affinity="nearest_neighbors")


#%%
t0 = time.time()
spectral.fit(X)
t1 = time.time()
if hasattr(spectral, 'labels_'):
    y_pred = spectral.labels_.astype(np.int)
else:
    y_pred = spectral.predict(X)



#%%
#Reshape y_pred
#Save as img
#91, 109, 91

#the challenge is to get the y_pred reshaped into the shape of the original whole brain.


y_pred_mask = np.reshape(y_pred(x,y))

img = nib.Nifti1Image(data, np.eye(4))
img.get_data_dtype() == np.dtype(np.int16)
nib.save(img, os.path.join('build','test4d.nii.gz'))






#%#

    blobs = generate_blobs()
    y_predict = cluster_timeseries(blobs, 3, similarity_metric = 'correlation')

#%%
tseries1=yeo_motor.T



tseries2=tseries1.sum(axis=1)

df=pd.DataFrame(tseries2)
df_norm = df.apply((df - df.mean()) / (df.max() - df.min()))

ism = individual_stability_matrix(tseries1, 5, 10)

#%%
tseries2pd.plot()
show()

plt.imshow(yeo_motor, cmap='hot', interpolation='nearest')
plt.show()









#%%
a = np.random.random((16, 16))
plt.imshow(A, cmap='hot', interpolation='nearest')
plt.show()

plt.imshow(tseries, cmap='hot', interpolation='nearest')
plt.show()

plt.imshow(S, cmap='hot', interpolation='nearest')
plt.show()