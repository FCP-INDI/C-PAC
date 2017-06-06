p# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 15:02:25 2017

@author: lbclab
"""
#Testing Yeo Resampling issues

import numpy as np
import nibabel as nb
import nilearn
from nilearn.image import resample_img
import matplotlib.pyplot as plt

template = load_mni152_template()
func_file = home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/Test_Data/residual_antswarp.nii.gz'

bg_file = home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/BasalGanglia_MNI2mm/BG.nii.gz'
yeo_7_lib_file = home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/Yeo_JNeurophysiol11_MNI152/FSL_Reorient2Std_Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz'

yeo_7_lib = nb.load(yeo_7_lib_file)

new_yeo2 = resample_img(yeo_7_lib, target_affine=template.affine, target_shape=template.shape[:3])
nb.save(new_yeo2, home + '/Dropbox/1_Projects/1_Research/2_CMI_BG_DEV/1_BASC/Data/new_yeo2.nii.gz')






#############UNNECESSARY###
img = resample_to_img(img1, template)
img = nilearn.image.resample_img(img1, target_affine=func.affine,target_shape=func.shape[:3], interpolation='nearest')
img = nilearn.image.resample_img(img1, target_shape=func.shape[:3], target_affine=np.diag((-2, 2, 2, 1)), interpolation='nearest')
img = nilearn.image.resample_img(img1, target_shape=func.shape[:3], target_affine=None, interpolation='nearest')



###EXAMPLE CODE

img_in_mm_space = resample_img(img, target_affine=np.eye(4),
                               target_shape=(512, 512, 1))

target_affine_3x3 = np.eye(3) * 2
target_affine_4x4 = np.eye(4) * 2
target_affine_4x4[3, 3] = 1.
img_3d_affine = resample_img(img, target_affine=target_affine_3x3)
img_4d_affine = resample_img(img, target_affine=target_affine_4x4)
target_affine_mm_space_offset_changed = np.eye(4)
target_affine_mm_space_offset_changed[:3, 3] = \
    img_3d_affine.affine[:3, 3]

img_3d_affine_in_mm_space = resample_img(
    img_3d_affine,
    target_affine=target_affine_mm_space_offset_changed,
    target_shape=(np.array(img_3d_affine.shape) * 2).astype(int))

img_4d_affine_in_mm_space = resample_img(
    img_4d_affine,
    target_affine=np.eye(4),
    target_shape=(np.array(img_4d_affine.shape) * 2).astype(int))