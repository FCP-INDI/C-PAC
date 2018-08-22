#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec

import os
import copy
import numpy as np
import nibabel as nb


def create_isfc(configuration, name='isfc'):
    """
    Inter-Subject Functional Correlation
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
        
    Returns
    -------
    workflow : nipype.pipeline.engine.Workflow
        ISFC workflow.
        
    Notes
    -----
    
    Workflow Inputs::
    
        
    Workflow Outputs::

    
    References
    ----------
    .. [1] Simony, E., Honey, C. J., Chen, J., Lositsky, O., Yeshurun,
           Y., Wiesel, A., & Hasson, U. (2016). Dynamic reconfiguration of the
           default mode network during narrative comprehension.
           Nature Communications, 7(May 2015), 1-13.
           https://doi.org/10.1038/ncomms12141
    
    """

    # inputspec = pe.Node(util.IdentityInterface(fields=[]),
    #                     name='inputspec')

    # outputspec = pe.Node(util.IdentityInterface(fields=[]),
    #                      name='outputspec')

    wf = pe.Workflow(name=name)

    from nilearn import datasets
    adhd_dataset = datasets.fetch_adhd(n_subjects=5)
    func_filenames = adhd_dataset.func
    
    wf.connect([
    ])

    return wf