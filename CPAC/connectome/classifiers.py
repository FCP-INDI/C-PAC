#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import (
    Ridge as _Ridge,
    LogisticRegression as _LogisticRegression
)
from sklearn import feature_selection


class Classifier:

    def fit(self, X, y=None):
        raise NotImplementedError

    def transform(self, X, y=None):
        raise NotImplementedError


class SVC(Classifier):

    def __init__(self, penalty='l2', anova_percentile=False):
        self.penalty = penalty
        self.anova_percentile = anova_percentile if anova_percentile > 0 else False

    def fit(self, X, y=None):

        pipeline = []

        if self.anova_percentile:
            feature_selection = feature_selection.SelectPercentile(
                feature_selection.f_classif,
                percentile=self.anova_percentile
            )
            pipeline += [('anova', feature_selection)]
        
        pipeline += [('svc', LinearSVC(penalty=self.penalty))]

        self.clf = Pipeline(pipeline)
        self.clf.fit(X, y) 

    def transform(self, X, y=None):
        return list(self.clf.predict(X)), y


class Ridge(Classifier):

    def fit(self, X, y=None):
        self.clf = _Ridge()
        self.clf.fit(X, y) 

    def transform(self, X, y=None):
        return list(self.clf.predict(X)), y


class LogisticRegression(Classifier):

    def __init__(self, penalty='l2'):
        self.penalty = penalty

    def fit(self, X, y=None):
        self.clf = _LogisticRegression(penalty=self.penalty)
        self.clf.fit(X, y) 

    def transform(self, X, y=None):
        return list(self.clf.predict(X)), y


class RandomForest(Classifier):

    def fit(self, X, y=None):
        self.clf = RandomForestClassifier()
        self.clf.fit(X, y) 

    def transform(self, X, y=None):
        return list(self.clf.predict(X)), y


class KNearestNeighbors(Classifier):

    def fit(self, X, y=None):
        self.clf = KNeighborsClassifier()
        self.clf.fit(X, y) 

    def transform(self, X, y=None):
        return list(self.clf.predict(X)), y


class NaiveBayes(Classifier):

    def fit(self, X, y=None):
        self.clf = GaussianNB()
        self.clf.fit(X, y) 

    def transform(self, X, y=None):
        return list(self.clf.predict(X)), y
    