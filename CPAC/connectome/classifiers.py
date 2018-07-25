#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from CPAC.connectome.cross_validation import CVedInterface
from sklearn.svm import SVC


class SVM(CVedInterface):

    def fit(self, X, y=None):
        clf = SVC()
        clf.fit(X, y) 
        return clf

    def transform(self, model, X, y=None):
        return model.predict(X), y
