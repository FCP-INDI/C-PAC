#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import nibabel as nib
from nilearn.decomposition import DictLearning as _DictLearning
from nilearn.decomposition import CanICA
from nilearn.regions import Parcellations

class ROI:

    def __init__(self, **parameters):
        self.parameters = parameters

    def fit(self, X, y=None):
        X = [nib.Nifti1Image(x, np.eye(len(x.shape))) for x in X]
        return self.fit_images(X, y)
    
    def fit_images(self, X, y=None):
        raise NotImplementedError

    def transform(self, model, X, y=None):
        X = [nib.Nifti1Image(x, np.eye(len(x.shape))) for x in X]
        return self.transform_images(X, y)

    def transform_images(self, model, X, y=None):
        raise NotImplementedError


class DictLearning(ROI):

    def fit_images(self, X, y=None):
        self.dict_learn = _DictLearning()
        self.dict_learn.fit(X)

    def transform_images(self, X, y=None):
        return [x.T for x in self.dict_learn.transform(X)], y


class KMeans(ROI):

    def fit_images(self, X, y=None):
        self.parc = Parcellations(method='kmeans', n_parcels=self.components)
        self.parc.fit(X)

    def transform_images(self, X, y=None):
        return [x.T for x in self.parc.transform(X)], y


class Ward(ROI):

    def fit_images(self, X, y=None):
        self.parc = Parcellations(method='ward', n_parcels=self.components)
        self.parc.fit(X)

    def transform_images(self, X, y=None):
        return [x.T for x in self.parc.transform(X)], y


class GroupICA(ROI):

    def fit_images(self, X, y=None):
        self.ica = CanICA(n_components=self.components)
        self.ica.fit(X)

    def transform_images(self, X, y=None):
        return [x.T for x in self.ica.transform(X)], y
