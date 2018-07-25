#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

class ClassifierInputSpec(BaseInterfaceInputSpec):
    subjects = traits.List(trait=traits.File(exists=True, desc='Model'),
                           desc='',
                           mandatory=True)
    model = traits.File(exists=True, desc='Model',
                        mandatory=False)
    vectorize = traits.Bool(desc='If True, connectivity matrices are reshaped into 1D '
                                 'arrays and only their flattened lower triangular '
                                 'parts are returned.', usedefault=True)


class ClassifierOutputSpec(TraitedSpec):
    connectivity = File(exists=True, desc="The connectivity matrix or vector")


class Classifier(BaseInterface):
    input_spec = ClassifierInputSpec
    output_spec = ClassifierOutputSpec

    def _run_interface(self, runtime):
        img = nb.load(self.inputs.functional)
        data = np.array(img.get_data())
        data = data.reshape((np.prod(data.shape[0:-1]), -1)).T

        # Normalize input for nilearn component
        metric = self.inputs.metric
        if metric == 'partial':
            metric = 'partial correlation'

        # Normalize input for file name
        metric_name = self.inputs.metric
        if metric_name == 'partial correlation':
            metric_name = 'partial'

        correlation_measure = ConnectivityMeasure(kind=metric,
                                                  vectorize=self.inputs.vectorize,
                                                  discard_diagonal=True)
        correlation_matrix = correlation_measure.fit_transform([data])[0]

        correlation_file = os.path.join(
            os.getcwd(), 'correlation_%s.npy' % metric_name)
        np.save(correlation_file, correlation_matrix)

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()

        # Normalize input for file name
        metric_name = self.inputs.metric
        if metric_name == 'partial correlation':
            metric_name = 'partial'

        correlation_file = os.path.join(
            os.getcwd(), 'correlation_%s.npy' % metric_name)
        outputs['connectivity'] = correlation_file
        return outputs
