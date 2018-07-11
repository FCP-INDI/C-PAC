#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile

import logging
logger = logging.getLogger()
logger.setLevel(logging.WARNING)

from nipype import config
config.enable_debug_mode()

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, File, traits, TraitedSpec


class TestConnectomeInputSpec(BaseInterfaceInputSpec):
    connectivity = File(exists=True)


class TestConnectome(BaseInterface):
    input_spec = TestConnectomeInputSpec

    def _run_interface(self, runtime):
        import numpy as np
        img = np.load(self.inputs.connectivity)

        assert img[0].mean() < 0.8
        assert img[:, 0].mean() < 0.8

        return runtime


def test_wf():

    from CPAC.connectome import create_connectivity

    import numpy as np
    import nibabel as nb

    _, image_path = tempfile.mkstemp(suffix='.nii.gz')
    
    data = np.ones((10, 10, 10, 100), dtype=np.float64)

    # Add same noise for each observation
    data += np.random.uniform(-0.01, 0.01, data.shape[-1])

    # Add into 0,0,0 a different noise
    data[0, 0, 0] += np.random.uniform(-0.01, 0.01, data.shape[-1])

    img = nb.Nifti1Image(data, np.eye(4))
    img.to_filename(image_path)

    wf = create_connectivity()
    wf.inputs.inputspec.functional = image_path
    wf.inputs.inputspec.metric = 'correlation'

    conn = pe.Node(TestConnectome(),
                   name='test_connectivity')

    outputspec = wf.get_node('outputspec')

    wf.connect(outputspec, 'connectivity',
               conn, 'connectivity')

    runtime = wf.run(plugin='Linear')

    os.remove(image_path)
    

class TestVectorShapeInputSpec(BaseInterfaceInputSpec):
    connectivity = File(exists=True)
    expected_shape = traits.Tuple()


class TestVectorShape(BaseInterface):
    input_spec = TestVectorShapeInputSpec

    def _run_interface(self, runtime):
        import numpy as np
        img = np.load(self.inputs.connectivity)

        assert img.shape == self.inputs.expected_shape

        return runtime


def test_wf_vector():

    from CPAC.connectome import create_connectivity

    import numpy as np
    import nibabel as nb

    _, image_path = tempfile.mkstemp(suffix='.nii.gz')
    
    data = np.random.uniform(-0.1, 0.1, (10, 10, 10, 100))

    img = nb.Nifti1Image(data, np.eye(4))
    img.to_filename(image_path)

    wf = create_connectivity()
    wf.inputs.inputspec.functional = image_path
    wf.inputs.inputspec.metric = 'correlation'
    wf.inputs.inputspec.vectorize = True

    variables = np.prod(data.shape[0:3])

    # Count only strict lower triangle
    expected_shape = variables * (variables - 1) / 2

    conn = pe.Node(TestVectorShape(),
                   name='test_connectivity')

    conn.inputs.expected_shape = (expected_shape, )

    outputspec = wf.get_node('outputspec')
    wf.connect(outputspec, 'connectivity',
               conn, 'connectivity')

    runtime = wf.run(plugin='Linear')

    os.remove(image_path)