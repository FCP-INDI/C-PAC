#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sklearn.svm import SVC as _SVC


class Classifier:

    def fit(self, X, y=None):
        raise NotImplementedError

    def transform(self, X, y=None):
        raise NotImplementedError


class SVC(Classifier):

    def __init__(self, **parameters):
        self.parameters = parameters

    def fit(self, X, y=None):
        self.clf = _SVC()
        self.clf.fit(X, y) 

    def transform(self, X, y=None):
        return list(self.clf.predict(X)), y
