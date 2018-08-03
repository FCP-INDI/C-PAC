#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import nibabel as nib
from nilearn.decomposition import DictLearning as _DictLearning


class Roi:

    def __init__(self, **parameters):
        self.parameters = parameters

    def fit(self, X, y=None):
        raise NotImplementedError

    def transform(self, model, X, y=None):
        raise NotImplementedError


class DictLearning(Roi):

    def fit(self, X, y=None):
        X = [nib.Nifti1Image(x, np.eye(len(x.shape))) for x in X]
        self.dict_learn = _DictLearning()
        self.dict_learn.fit(X)

    def transform(self, X, y=None):
        X = [nib.Nifti1Image(x, np.eye(len(x.shape))) for x in X]
        return [x.T for x in self.dict_learn.transform(X)], y
