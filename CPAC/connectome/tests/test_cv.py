#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import tempfile
import shutil

import logging
logging.basicConfig(level=logging.DEBUG)

from nipype import config
config.enable_debug_mode()

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, TraitedSpec


def test_wf():

    from CPAC.connectome.cross_validation import CVInterface
    from CPAC.connectome.classifiers import SVMInterface
    from CPAC.connectome.rois import DictLearningInterface

    import numpy as np
    import nibabel as nb
    from nilearn import datasets

    wf = pe.Workflow(name='test')

    adhd_dataset = datasets.fetch_adhd(n_subjects=5)
    func_filenames = adhd_dataset.func
    
    cv = pe.Node(CVInterface(), name='cv')
    cv.inputs.folds = 2
    # cv.inputs.X = list(nb.load(f).get_data() for f in func_filenames)
    # cv.inputs.X = list(np.ones((2, 2, 2, 3)) for f in func_filenames)
    cv.inputs.X = func_filenames
    cv.inputs.y = [1, 2, 2, 1, 2]
    # cv.synchronize = True

    dictlearn = DictLearningInterface()
    dictlearn_train = pe.MapNode(interface=dictlearn, name='dictlearn_train', iterfield=['X', 'y', 'fold'])
    dictlearn_valid = pe.MapNode(interface=dictlearn, name='dictlearn_valid', iterfield=['X', 'y', 'fold', 'model'])

    svm = SVMInterface()
    svm_train = pe.MapNode(interface=svm, name='svm_train', iterfield=['X', 'y', 'fold'])
    svm_valid = pe.MapNode(interface=svm, name='svm_valid', iterfield=['X', 'y', 'fold', 'model'])

    wf.connect(cv, 'train_X', dictlearn_train, 'X')
    wf.connect(cv, 'train_y', dictlearn_train, 'y')
    wf.connect(cv, 'fold', dictlearn_train, 'fold')

    wf.connect(cv, 'valid_X', dictlearn_valid, 'X')
    wf.connect(cv, 'valid_y', dictlearn_valid, 'y')
    wf.connect(cv, 'fold', dictlearn_valid, 'fold')

    wf.connect(dictlearn_train, 'model', dictlearn_valid, 'model')

    wf.connect(dictlearn_train, 'X', svm_train, 'X')
    wf.connect(dictlearn_train, 'y', svm_train, 'y')
    wf.connect(dictlearn_train, 'fold', svm_train, 'fold')

    wf.connect(dictlearn_valid, 'X', svm_valid, 'X')
    wf.connect(dictlearn_valid, 'y', svm_valid, 'y')
    wf.connect(dictlearn_valid, 'fold', svm_valid, 'fold')

    wf.connect(svm_train, 'model', svm_valid, 'model')

    wf.base_dir = '/tmp/cv'

    wf.write_graph(graph2use='exec')

    shutil.rmtree('/tmp/cv')

    runtime = wf.run(plugin='Linear')

    