#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 19:09:11 2017

@author: aki.nikolaidis
"""
#%%
import matplotlib
import numpy as np
import nibabel as nb
import pandas as pd
import nilearn.image as image

from nilearn import datasets
from nilearn.input_data import NiftiMasker
from nilearn.plotting import plot_roi, show
from nilearn.image.image import mean_img
from nilearn.image import resample_img
from matplotlib import pyplot as plt

matplotlib.style.use('ggplot')

#%%
#Data Preparation

#load in dataset/func file- these are exchangable. (make sure it's in MNI152 space)
dataset = datasets.fetch_adhd(n_subjects=2)
func_file = '/Users/aki.nikolaidis/nilearn_data/residual_antswarp.nii.gz'
func = nb.load(func_file)


#load masks
bg_file = '/Users/aki.nikolaidis/Dropbox/1_Projects/1_Research/2_BASC/Data/BasalGanglia_MNI2mm/BG.nii.gz'
yeo_7_lib_file = '/Users/aki.nikolaidis/Dropbox/1_Projects/1_Research/2_BASC/Data/Yeo_JNeurophysiol11_MNI152/FSL_Reorient2Std_Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz'

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
bg = np.resize(bg, (new_shape,1))
yeo_7_lib = np.resize(yeo_7_lib, (new_shape,1))


#%%
#Extract masks from functional data
#Action- add remainder networks
yeo_motor = func[np.where( yeo_7_lib == 2 )[0], :]
yeo_vis = func[np.where( yeo_7_lib == 1 )[0], :]
bg_func = func[np.where( bg == 1 )[0], :]

print 'motor', type(yeo_motor), yeo_motor.shape, 'vis', yeo_vis.shape, 'bg', bg_func.shape
#%%
tseries1=yeo_motor.T

#Hey Stranger....

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