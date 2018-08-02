#!/usr/bin/env python
# -*- coding: utf-8 -*-

from CPAC.connectome.cross_validation import CVedInterface
from nilearn.decomposition import DictLearning
import nibabel as nib
import numpy as np


class DictLearningInterface(CVedInterface):

    def fit(self, X, y=None):
        X = [nib.Nifti1Image(x, np.eye(len(x.shape))) for x in X]
        dict_learn = DictLearning()
        dict_learn.fit(X)
        return dict_learn

    def transform(self, model, X, y=None):
        X = [nib.Nifti1Image(x, np.eye(len(x.shape))) for x in X]
        return [x.T for x in model.transform(X)], y
