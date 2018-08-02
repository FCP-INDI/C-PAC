#!/usr/bin/env python
# -*- coding: utf-8 -*-

from CPAC.connectome.cross_validation import CVedInterface
from sklearn.svm import SVC


class SVMInterface(CVedInterface):

    def fit(self, X, y=None):
        clf = SVC()
        clf.fit(X, y) 
        return clf

    def transform(self, model, X, y=None):
        return list(model.predict(X)), y
