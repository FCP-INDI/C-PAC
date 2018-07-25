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
    BaseInterfaceInputSpec, traits, TraitedSpec


def test_wf():

    from CPAC.connectome.cross_validation import CVInterface

    import numpy as np
    import nibabel as nb

    wf = pe.Workflow(name='test')

    cv = pe.Node(CVInterface(), name='cv')
    cv.inputs.folds = 5
    cv.inputs.X = [
        'a1', 'a2', 'a3',
        'b1', 'b2', 'b3',
        'c1', 'c2', 'c3',
        'd1', 'd2', 'd3',
        'e1', 'e2', 'e3',
    ]
    cv.inputs.y = ['asd'] * len(cv.inputs.X)
    cv.synchronize = True

    def std_func(X, y, fold, model=None):
        print('Std ' + ('Train ' if not model else 'Test ') + str(fold))
        print(X, y)
        return (X, y, fold, model if model else True)

    std = util.Function(
        input_names=['X', 'y', 'fold', 'model'],
        output_names=['X', 'y', 'fold', 'model'],
        function=std_func
    )

    def svm_func(X, y, fold, model=None):
        print('SVM ' + ('Train ' if not model else 'Test ') + str(fold))
        print(X, y)
        return (X, y, fold, model if model else True)

    svm = util.Function(
        input_names=['X', 'y', 'fold', 'model'],
        output_names=['X', 'y', 'fold', 'model'],
        function=svm_func
    )

    std_train = pe.MapNode(interface=std, name='std_train', iterfield=['X', 'y', 'fold'])
    std_valid = pe.MapNode(interface=std, name='std_valid', iterfield=['X', 'y', 'fold', 'model'])

    svm_train = pe.MapNode(interface=svm, name='svm_train', iterfield=['X', 'y', 'fold'])
    svm_valid = pe.MapNode(interface=svm, name='svm_valid', iterfield=['X', 'y', 'fold', 'model'])

    wf.connect(cv, 'train_X', std_train, 'X')
    wf.connect(cv, 'train_y', std_train, 'y')
    wf.connect(cv, 'fold', std_train, 'fold')

    wf.connect(cv, 'valid_X', std_valid, 'X')
    wf.connect(cv, 'valid_y', std_valid, 'y')
    wf.connect(cv, 'fold', std_valid, 'fold')

    wf.connect(std_train, 'model', std_valid, 'model')

    wf.connect(std_train, 'X', svm_train, 'X')
    wf.connect(std_train, 'y', svm_train, 'y')
    wf.connect(std_train, 'fold', svm_train, 'fold')

    wf.connect(std_valid, 'X', svm_valid, 'X')
    wf.connect(std_valid, 'y', svm_valid, 'y')
    wf.connect(std_valid, 'fold', svm_valid, 'fold')

    wf.connect(svm_train, 'model', svm_valid, 'model')

    wf.base_dir = '/tmp/cv'

    runtime = wf.run(plugin='Linear')

    