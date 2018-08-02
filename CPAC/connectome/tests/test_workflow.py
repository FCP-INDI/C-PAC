#!/usr/bin/env python
# -*- coding: utf-8 -*-

import yaml

import logging
logger = logging.getLogger()
logger.setLevel(logging.WARNING)

from nipype import config
config.enable_debug_mode()


def test_wf():

    from CPAC.connectome import create_connectome_study

    config = yaml.load("""

connectome:

  roi:
    - type: dict_learning
      parameters:
        components: 40
    - type: dict_learning
      parameters:
        components: 60
    - type: dict_learning
      parameters:
        components: 80

  connectivity: 
    - type: partial
    - type: tangent
    - type: correlation

  classifier:

    - type: svc
      parameters:
        penalty: l2
        anova_percentile: 10

    - type: svc
      parameters:
        penalty: l1
    
    """)

    wf = create_connectome_study(config['connectome'])
