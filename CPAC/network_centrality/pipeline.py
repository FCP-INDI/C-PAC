import os
from nipype import logging
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl

from CPAC.utils.interfaces.function import Function
from CPAC.network_centrality.network_centrality import create_centrality_wf
from CPAC.network_centrality.utils import merge_lists, check_centrality_params

logger = logging.getLogger('workflow')


def create_network_centrality_workflow(workflow, c, strategies):

    if not any((
        True in c.degWeightOptions,
        True in c.eigWeightOptions,
        True in c.lfcdWeightOptions
    )):
        return strategies


    for num_strat, strat in enumerate(strategies[:]):

        # Resample the functional mni to the centrality mask resolution
        resample_functional_to_template = pe.Node(
            interface=fsl.FLIRT(),
            name='resample_functional_to_template_%d' % num_strat

        )
        resample_functional_to_template.inputs.set(
            interp='trilinear',
            in_matrix_file=c.identityMatrix,
            apply_xfm=True
        )

        node, out_file = strat['functional_to_standard']

        workflow.connect(node, out_file,
                         resample_functional_to_template, 'in_file')

        workflow.connect(c.templateSpecificationFile, 'local_path',
                         resample_functional_to_template, 'reference')

        merge_node = pe.Node(Function(input_names=['deg_list',
                                                   'eig_list',
                                                   'lfcd_list'],
                                      output_names=['merged_list'],
                                      function=merge_lists,
                                      as_module=True),
                             name='merge_node_%d' % num_strat)

        if True in c.degWeightOptions:
            connect_centrality_workflow(
                workflow, c, strat, num_strat,
                resample_functional_to_template,
                c.templateSpecificationFile, merge_node,
                'degree',
                c.degCorrelationThresholdOption,
                c.degCorrelationThreshold
            )
            
        if True in c.eigWeightOptions:
            connect_centrality_workflow(
                workflow, c, strat, num_strat,
                resample_functional_to_template,
                c.templateSpecificationFile, merge_node,
                'eigenvector',
                c.eigCorrelationThresholdOption,
                c.eigCorrelationThreshold
            )
                
        if True in c.lfcdWeightOptions:
            connect_centrality_workflow(
                workflow, c, strat, num_strat,
                resample_functional_to_template,
                c.templateSpecificationFile, merge_node,
                'lfcd',
                c.lfcdCorrelationThresholdOption,
                c.lfcdCorrelationThreshold
            )

        if 0 in c.runNetworkCentrality:
            strategies += [strat.fork()]

        strat.update_resource_pool({
            'centrality': (merge_node, 'merged_list')
        })

    return strategies


# Function to connect the into pipeline
def connect_centrality_workflow(workflow, c, strat, num_strat,
                                resample_functional_to_template,
                                template,
                                merge_node,
                                method_option, threshold_option,
                                threshold):

    # Set method_options variables
    if method_option == 'degree':
        out_list = 'deg_list'
    elif method_option == 'eigenvector':
        out_list = 'eig_list'
    elif method_option == 'lfcd':
        out_list = 'lfcd_list'

    # Init workflow name and resource limits
    wf_name = 'afni_centrality_%d_%s' % (num_strat, method_option)
    num_threads = c.maxCoresPerParticipant
    memory = c.memoryAllocatedForDegreeCentrality

    # Format method and threshold options properly and check for
    # errors
    method_option, threshold_option = \
        check_centrality_params(method_option,
                                threshold_option,
                                threshold)

    # Change sparsity thresholding to % to work with afni
    if threshold_option == 'sparsity':
        threshold = threshold * 100

    afni_centrality_wf = \
        create_centrality_wf(wf_name, method_option,
                             threshold_option,
                             threshold, num_threads, memory)

    workflow.connect(resample_functional_to_template, 'out_file',
                     afni_centrality_wf, 'inputspec.in_file')

    workflow.connect(template, 'local_path',
                     afni_centrality_wf, 'inputspec.template')

    workflow.connect(afni_centrality_wf,
                     'outputspec.outfile_list',
                     merge_node,
                     out_list)
