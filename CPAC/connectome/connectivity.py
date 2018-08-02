#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from nipype.interfaces.base import traits
from nilearn.connectome import ConnectivityMeasure
from CPAC.connectome.cross_validation import CVedInterface, CVedInputSpec, CVedOutputSpec


class FunctionalConnectivityInputSpec(CVedInputSpec):
    metric = traits.Enum('correlation', 'partial correlation', 'partial', 'tangent',
                         desc='The measure to be used to compute the connectivity matrix',
                         mandatory=True)
    vectorize = traits.Bool(desc='If True, connectivity matrices are reshaped into 1D '
                                 'arrays and only their flattened lower triangular '
                                 'parts are returned.', usedefault=True)


class FunctionalConnectivity(CVedInterface):
    input_spec = FunctionalConnectivityInputSpec

    def fit(self, X, y=None):
        X = [x.reshape((np.prod(x.shape[0:-1]), -1)).T for x in X]

        metric = self.inputs.metric
        if metric == 'partial':
            metric = 'partial correlation'

        correlation_measure = ConnectivityMeasure(kind=metric,
                                                  vectorize=self.inputs.vectorize,
                                                  discard_diagonal=True)
        correlation_measure.fit(X)
        return correlation_measure

    def transform(self, model, X, y=None):
        X = [x.reshape((np.prod(x.shape[0:-1]), -1)).T for x in X]
        return list(model.transform(X)), y
