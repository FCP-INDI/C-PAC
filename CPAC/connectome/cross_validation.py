#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import copy
import numpy as np
import nibabel as nb
from collections import OrderedDict

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from sklearn.preprocessing import LabelEncoder
from nilearn.connectome import ConnectivityMeasure

from CPAC.connectome.rois import (
    DictLearning,
    KMeans,
    Ward,
    GroupICA
)

from CPAC.connectome.classifiers import (
    SVC,
    Ridge,
    LogisticRegression,
    RandomForest,
    KNearestNeighbors,
    NaiveBayes,
)

from nipype import logging
logger = logging.getLogger('interface')


def _load(x):
    if isinstance(x, six.string_types):
        if x.endswith('.npy'):
            return np.load(x)
        return nb.load(x).get_data()
    else:
        return x


def _normalize(X):
    return [x.reshape((np.prod(x.shape[0:-1]), -1)).T for x in X]


class CVInputSpec(TraitedSpec):
    folds = traits.Int(mandatory=True)

    X = traits.List(traits.Any, mandatory=True)
    y = traits.List(traits.Any, mandatory=False)

    stratification = traits.List(trait=traits.Any,
                                 mandatory=False)


class CVOutputSpec(TraitedSpec):

    label_encoder = traits.Instance(klass=LabelEncoder)

    fold = traits.List(trait=traits.Int)
    config = traits.List(trait=traits.Any)

    X_train = traits.List(traits.List(trait=traits.Any))
    y_train = traits.List(traits.List(trait=traits.Int))
    X_valid = traits.List(traits.List(trait=traits.Any))
    y_valid = traits.List(traits.List(trait=traits.Int))


class CVInterface(BaseInterface):
    input_spec = CVInputSpec
    output_spec = CVOutputSpec

    def _run_interface(self, runtime):
        from sklearn.model_selection import KFold
        from sklearn.preprocessing import LabelEncoder

        X = self.inputs.X
        y = self.inputs.y

        self._le = None
        if hasattr(y, "__getitem__") and hasattr(y, "__iter__"):

            assert len(y) == len(X)

            self._le = LabelEncoder()
            self._le.fit(y)

        kf = KFold(n_splits=self.inputs.folds, shuffle=False)
        self._folds = list(kf.split(range(len(X))))

        return runtime

    def _list_outputs(self):

        # Select multi indexes from a list
        from operator import itemgetter 
        g = lambda l, i: list(map(l.__getitem__, i))

        X = self.inputs.X
        y = self.inputs.y

        outputs = self._outputs().get()

        if hasattr(y, "__getitem__") and hasattr(y, "__iter__"):
            y = self._le.transform(y)
            outputs['y_train'] = [y[d[0]].tolist() for d in self._folds]
            outputs['y_valid'] = [y[d[1]].tolist() for d in self._folds]
            outputs['label_encoder'] = self._le

        outputs['X_train'] = [g(X, d[0]) for d in self._folds]
        outputs['X_valid'] = [g(X, d[1]) for d in self._folds]
        outputs['fold'] = list(range(self.inputs.folds))
        outputs['config'] = [OrderedDict() for d in self._folds]
        return outputs


class CVedInputSpec(TraitedSpec):
    X_train = traits.List(traits.Any, mandatory=True)
    y_train = traits.Array(mandatory=False)
    X_valid = traits.List(traits.Any, mandatory=True)
    y_valid = traits.Array(mandatory=False)

    fold = traits.Int(mandatory=True)
    config = traits.Any(mandatory=False)
    node_config = traits.Any(mandatory=False)


class CVedOutputSpec(TraitedSpec):
    X_train = traits.List(traits.Any, mandatory=True)
    y_train = traits.Array(mandatory=False)
    X_valid = traits.List(traits.Any, mandatory=True)
    y_valid = traits.Array(mandatory=False)

    fold = traits.Int(mandatory=True)
    config = traits.Any(mandatory=False)


class CVedInterface(BaseInterface):
    input_spec = CVedInputSpec
    output_spec = CVedOutputSpec

    _node_outputs = None

    def fit(self, X, y=None):
        raise NotImplementedError

    def transform(self, model, X, y=None):
        raise NotImplementedError

    def _run_interface(self, runtime):
        self._node_outputs = self._outputs().get()

        logger.info("Loading data on %s %s (fold %d)" % (self.__class__.__name__, self.inputs.node_config['type'], self.inputs.fold))
        X_train = [_load(x) for x in self.inputs.X_train]
        y_train = self.inputs.y_train
        X_valid = [_load(x) for x in self.inputs.X_valid]
        y_valid = self.inputs.y_valid

        self._model = self.fit(X_train, y_train)
    
        self._node_outputs['X_train'], self._node_outputs['y_train'] = \
            self.transform(self._model, X_train, y_train)

        self._node_outputs['X_valid'], self._node_outputs['y_valid'] = \
            self.transform(self._model, X_valid, y_valid)

        return runtime

    def _list_outputs(self):

        logger.info("Outputing %s (fold %d)" % (self.__class__.__name__, self.inputs.fold))
        
        outputs = self._node_outputs \
                  if self._node_outputs is not None \
                  else self._outputs().get()

        outputs['fold'] = self.inputs.fold
        outputs['config'] = self.inputs.config
        outputs['config'].update({
            "fold": self.inputs.fold,
            self.config_key: copy.deepcopy(self.inputs.node_config)
        })

        del self._node_outputs

        return outputs


class ROI(CVedInterface):

    config_key = "roi"

    def factory(self, config):
        methods = {
            'dict_learning': DictLearning,
            'k_means': KMeans,
            'kmeans': KMeans,
            'ward': Ward,
            'group_ica': GroupICA,
        }
        config['parameters'] = {} if not config['parameters'] else config['parameters']
        return methods[config['type']](**config['parameters'])

    def fit(self, X, y=None):
        config = copy.deepcopy(self.inputs.node_config)
        model = self.factory(config)
        model.fit(X, y)
        return model

    def transform(self, model, X, y=None):
        return model.transform(X, y)


class Connectivity(CVedInterface):

    config_key = "connectivity"

    def fit(self, X, y=None):
        X = _normalize(X)

        metric = self.inputs.node_config['type']
        if metric == 'partial' or metric == 'partial_correlation':
            metric = 'partial correlation'

        correlation_measure = ConnectivityMeasure(kind=metric,
                                                  vectorize=True,
                                                  discard_diagonal=True)
        correlation_measure.fit(X)
        return correlation_measure

    def transform(self, model, X, y=None):
        X = _normalize(X)
        return list(model.transform(X)), y


class Classifier(CVedInterface):

    config_key = "classifier"

    def factory(self, config):
        methods = {
            'svc': SVC,
            'svm': SVC,
            'ridge': Ridge,
            'logistic': LogisticRegression,
            'random_forest': RandomForest,
            'knn': KNearestNeighbors,
            'knearestneighbors': KNearestNeighbors,
            'k_nearest_neighbors': KNearestNeighbors,
            'naive_bayes': NaiveBayes,
        }
        config['parameters'] = {} if not config['parameters'] else config['parameters']
        return methods[config['type']](**config['parameters'])

    def fit(self, X, y=None):
        config = copy.deepcopy(self.inputs.node_config)
        model = self.factory(config)
        model.fit(X, y)
        return model

    def transform(self, model, X, y=None):
        return model.transform(X, y)


class JoinerInputSpec(TraitedSpec):
    X_train = traits.Any(mandatory=True)
    y_train = traits.Any(mandatory=True)
    X_valid = traits.Any(mandatory=True)
    y_valid = traits.Any(mandatory=True)

    fold = traits.Any(mandatory=True)
    config = traits.Any(mandatory=True)

def detect_depth_on_folds(folds):
    start = 0
    while not isinstance(folds, (int)):
        start += 1
        folds = folds[0]
    return start

def flatten(l):
    return [item for sublist in l for item in sublist]

class JoinerInterface(BaseInterface):
    input_spec = JoinerInputSpec
    output_spec = JoinerInputSpec

    def _run_interface(self, runtime):
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()

        X_train = self.inputs.X_train
        y_train = self.inputs.y_train
        X_valid = self.inputs.X_valid
        y_valid = self.inputs.y_valid
        fold = self.inputs.fold
        config = self.inputs.config

        for _ in range(detect_depth_on_folds(fold) - 1):
            X_train = flatten(X_train)
            y_train = flatten(y_train)
            X_valid = flatten(X_valid)
            y_valid = flatten(y_valid)
            fold = flatten(fold)
            config = flatten(config)

        outputs['X_train'] = X_train
        outputs['y_train'] = y_train
        outputs['X_valid'] = X_valid
        outputs['y_valid'] = y_valid
        outputs['fold']    = fold
        outputs['config']  = config
 
        return outputs
