#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import numpy as np
import nibabel as nb
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from sklearn.preprocessing import LabelEncoder


class CVInputSpec(BaseInterfaceInputSpec):
    folds = traits.Int(mandatory=True)

    X = traits.List(traits.Any, mandatory=True)
    y = traits.List(traits.Any, mandatory=False)

    stratification = traits.List(trait=traits.Any,
                                 mandatory=False)


class CVOutputSpec(TraitedSpec):

    label_encoder = traits.Instance(klass=LabelEncoder)

    fold = traits.List(trait=traits.Int)

    train_X = traits.List(traits.List(trait=traits.Any))
    train_y = traits.List(traits.List(trait=traits.Int))

    valid_X = traits.List(traits.List(trait=traits.Any))
    valid_y = traits.List(traits.List(trait=traits.Int))


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
        g = lambda l, i: list(itemgetter(*i)(l))

        X = self.inputs.X
        y = self.inputs.y

        outputs = self._outputs().get()

        if hasattr(y, "__getitem__") and hasattr(y, "__iter__"):
            y = self._le.transform(y)
            outputs['train_y'] = [g(y, d[0]) for d in self._folds]
            outputs['valid_y'] = [g(y, d[1]) for d in self._folds]
            outputs['label_encoder'] = self._le

        outputs['train_X'] = [g(X, d[0]) for d in self._folds]
        outputs['valid_X'] = [g(X, d[1]) for d in self._folds]
        outputs['fold'] = list(range(self.inputs.folds))

        return outputs


class CVedInputSpec(BaseInterfaceInputSpec):
    X = traits.List(traits.Any, mandatory=True)
    y = traits.Array(mandatory=False)

    fold = traits.Int(mandatory=True)
    model = traits.Any(mandatory=False)


class CVedOutputSpec(TraitedSpec):
    X = traits.List(traits.Any, mandatory=True)
    y = traits.Array(mandatory=False)

    fold = traits.Int(mandatory=True)
    model = traits.Any(mandatory=False)


class CVedInterface(BaseInterface):
    input_spec = CVedInputSpec
    output_spec = CVedOutputSpec

    _node_outputs = None

    def fit(self, X, y=None):
        raise NotImplementedError

    def transform(self, X, y=None):
        raise NotImplementedError

    def _run_interface(self, runtime):
        from nipype import logging
        logger = logging.getLogger('interface')

        self._model = self.inputs.model
        self._node_outputs = self._outputs().get()

        logger.info("Loading data on %s (fold %d)" % (self.__class__.__name__, self.inputs.fold))
        X = [_load(x) for x in self.inputs.X]
        y = self.inputs.y

        if not self._model:
            logger.info("Fitting %s (fold %d)" % (self.__class__.__name__, self.inputs.fold))
            self._model = self.fit(X, y)
        else:
            logger.info("Validating %s (fold %d)" % (self.__class__.__name__, self.inputs.fold))

        self._node_outputs['X'], self._node_outputs['y'] = \
            self.transform(self._model, X, y)

        return runtime

    def _list_outputs(self):
        outputs = self._node_outputs \
                  if self._node_outputs is not None \
                  else self._node_outputs().get()

        outputs['model'] = self._model
        outputs['fold'] = self.inputs.fold

        del self._node_outputs

        return outputs


def _load(x):
    if isinstance(x, six.string_types):
        if x.endswith('.npy'):
            return np.load(x)
        return nb.load(x).get_data()
    else:
        return x
