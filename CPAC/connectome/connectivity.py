#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from nilearn.connectome import ConnectivityMeasure
from CPAC.connectome.cross_validation import CVedInterface


def _normalize(X):
    return [x.reshape((np.prod(x.shape[0:-1]), -1)).T for x in X]


class FunctionalConnectivity(CVedInterface):

    def fit(self, X, y=None):
        X = _normalize(X)

        metric = self.metric
        if metric == 'partial':
            metric = 'partial correlation'

        correlation_measure = ConnectivityMeasure(kind=metric,
                                                  vectorize=True,
                                                  discard_diagonal=True)
        correlation_measure.fit(X)
        return correlation_measure

    def transform(self, model, X, y=None):
        X = _normalize(X)
        return list(model.transform(X)), y


class CorrelationConnectivityInterface(FunctionalConnectivity):
    metric = 'correlation'


class PartialCorrelationConnectivityInterface(FunctionalConnectivity):
    metric = 'partial correlation'


class TangentConnectivityInterface(FunctionalConnectivity):
    metric = 'tangent'