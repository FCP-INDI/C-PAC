#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import numpy as np
import nibabel as nb
import pandas as pd
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

from sklearn import metrics
from sklearn.preprocessing import LabelEncoder


def serialize_parameters(parameters):
    return '; '.join(
        '%s: %s' % (k, parameters[k])
        for k in sorted(parameters.keys())
    )

class ReportInputSpec(TraitedSpec):
    X_train = traits.Any(mandatory=True)
    y_train = traits.Any(mandatory=True)
    X_valid = traits.Any(mandatory=True)
    y_valid = traits.Any(mandatory=True)

    fold = traits.Any(mandatory=True)
    config = traits.Any(mandatory=True)

    label_encoder = traits.Instance(klass=LabelEncoder)


class ReportOutputSpec(TraitedSpec):
    data = traits.Instance(klass=pd.DataFrame)


class ReportInterface(BaseInterface):
    input_spec = ReportInputSpec
    output_spec = ReportOutputSpec

    def _run_interface(self, runtime):
        from nipype import logging
        logger = logging.getLogger('interface')

        combination = zip(self.inputs.X_train, self.inputs.y_train,
                          self.inputs.X_valid, self.inputs.y_valid,
                          self.inputs.fold, self.inputs.config)

        rows = []

        for X_train, y_train, X_valid, y_valid, fold, config in combination:
            row = {}

            if self.inputs.label_encoder:
                le = self.inputs.label_encoder
                if len(le.classes_) == 2:
                    row['train tn'], \
                    row['train fp'], \
                    row['train fn'], \
                    row['train tp'] = metrics.confusion_matrix(np.round(y_train), np.round(X_train)).ravel()

                    row['valid tn'], \
                    row['valid fp'], \
                    row['valid fn'], \
                    row['valid tp'] = metrics.confusion_matrix(np.round(y_valid), np.round(X_valid)).ravel()

            row['train roc auc'] = metrics.roc_auc_score(y_train, X_train)
            row['valid roc auc'] = metrics.roc_auc_score(y_valid, X_valid)
                
            row['roi'] = config['roi']['type']
            row['roi parameters'] = \
                '' if not config['roi']['parameters'] \
                else serialize_parameters(config['roi']['parameters'][k])

            row['connectivity'] = config['connectivity']['type']
            row['connectivity parameters'] = \
                '' if not config['connectivity']['parameters'] \
                else serialize_parameters(config['connectivity']['parameters'][k])
                

            row['classifier'] = config['classifier']['type']
            row['classifier parameters'] = \
                '' if not config['classifier']['parameters'] \
                else serialize_parameters(config['classifier']['parameters'][k])

            rows += [row]

        
        self._df = pd.DataFrame(rows)
                
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['data'] = self._df

        self._df.to_csv('./report.csv')

        return outputs

