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

  folds: 2

  roi:
    - type: aal
      parameters:

    # - type: harvard_oxford
    #   parameters:

    # - type: power
    #   parameters:

    # - type: basc
    #   parameters:

    # - type: ward
    #   parameters: 40

    # - type: ward
    #   parameters: 60

    # - type: ward
    #   parameters: 80

    # - type: ward
    #   parameters: 100

    # - type: ward
    #   parameters: 120

    # - type: ward
    #   parameters: 150

    # - type: ward
    #   parameters: 200

    # - type: ward
    #   parameters: 300

    # - type: k_means
    #   parameters: 40

    # - type: k_means
    #   parameters: 60

    # - type: k_means
    #   parameters: 80

    # - type: k_means
    #   parameters: 100

    # - type: k_means
    #   parameters: 120

    # - type: k_means
    #   parameters: 150

    # - type: k_means
    #   parameters: 200

    # - type: k_means
    #   parameters: 300

    # - type: group_ica
    #   parameters:
    #     components: 40
    # - type: group_ica
    #   parameters:
    #     components: 60
    # - type: group_ica
    #   parameters:
    #     components: 80
    # - type: group_ica
    #   parameters:
    #     components: 100
    # - type: group_ica
    #   parameters:
    #     components: 120

    # - type: dict_learning
    #   parameters:
    #     components: 40
    # - type: dict_learning
    #   parameters:
    #     components: 60
    # - type: dict_learning
    #   parameters:
    #     components: 80
    # - type: dict_learning
    #   parameters:
    #     components: 100
    # - type: dict_learning
    #   parameters:
    #     components: 120

  connectivity: 
    - type: partial
      parameters:

    # - type: tangent
    #   parameters:

    # - type: correlation
    #   parameters:

  classifier:

    - type: knn
      parameters:

    # - type: random_forest
    #   parameters:

    # - type: naive_bayes
    #   parameters:

    # - type: ridge
    #   parameters:

    # - type: svc
    #   parameters:
    #     penalty: l1
    #     anova_percentile: 10
          
    # - type: svc
    #   parameters:
    #     penalty: l2
    #     anova_percentile: 10

    # - type: svc
    #   parameters:
    #     penalty: l1

    # - type: svc
    #   parameters:
    #     penalty: l2

    # - type: logistic
    #   parameters:
    #     penalty: l1

    # - type: logistic
    #   parameters:
    #     penalty: l2

    """)

    wf = create_connectome_study(config['connectome'])
    wf.write_graph(simple_form=False, graph2use='exec', format='png')
    wf.run()