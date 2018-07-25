#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from sklearn.preprocessing import LabelEncoder


class CVInputSpec(BaseInterfaceInputSpec):
    folds = traits.Int(mandatory=True)

    X = traits.List(trait=traits.String,
                           mandatory=True)

    y = traits.List(trait=traits.String,
                    mandatory=False)

    stratification = traits.List(trait=traits.String,
                                 mandatory=False)


class CVOutputSpec(TraitedSpec):

    label_encoder = traits.Instance(klass=LabelEncoder)

    fold = traits.List(trait=traits.Int)

    train_X = traits.List(trait=traits.List(trait=traits.String))
    train_y = traits.List(trait=traits.List(trait=traits.Int))

    valid_X = traits.List(trait=traits.List(trait=traits.String))
    valid_y = traits.List(trait=traits.List(trait=traits.Int))


class CVInterface(BaseInterface):
    input_spec = CVInputSpec
    output_spec = CVOutputSpec

    def _run_interface(self, runtime):
        from sklearn.model_selection import KFold
        from sklearn.preprocessing import LabelEncoder

        X = self.inputs.X
        y = self.inputs.y

        self._le = None
        if y:
            self._le = LabelEncoder()
            self._le.fit(y)

        kf = KFold(n_splits=self.inputs.folds)
        self._folds = list(kf.split(X))

        return runtime

    def _list_outputs(self):

        from operator import itemgetter 

        g = lambda l, i: list(itemgetter(*i)(l))
        X = self.inputs.X
        y = self.inputs.y

        outputs = self._outputs().get()

        if y:
            y = self._le.transform(y)
            outputs['train_y'] = [g(y, d[0]) for d in self._folds]
            outputs['valid_y'] = [g(y, d[1]) for d in self._folds]
            outputs['label_encoder'] = self._le
        outputs['train_X'] = [g(X, d[0]) for d in self._folds]
        outputs['valid_X'] = [g(X, d[1]) for d in self._folds]
        outputs['fold'] = list(range(self.inputs.folds))

        return outputs


class CVedInputSpec(BaseInterfaceInputSpec):
    X = traits.List(trait=traits.File(exists=True),
                    mandatory=True)
    y = traits.List(trait=traits.Int,
                    mandatory=False)

    model = trait=traits.File(exists=True, mandatory=False)


class CVedOutputSpec(TraitedSpec):
    X = traits.List(trait=traits.File(exists=True),
                    mandatory=True)
    y = traits.List(trait=traits.Int,
                    mandatory=False)

    model = trait=traits.File(exists=True, mandatory=False)


class CVedInterface(BaseInterface):
    input_spec = CVedInputSpec
    output_spec = CVedOutputSpec

    def fit(self, X, y=None):
        raise NotImplementedError

    def transform(self, X, y=None):
        raise NotImplementedError

    def _run_interface(self, runtime):
        self._model = self.inputs.model

        if not model:
            self._model = self.fit(self.inputs.X, self.inputs.y)

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['X'], outputs['y'] = self._model.transform(self.inputs.X, self.inputs.y)
        outputs['model'] = self._model
        return outputs