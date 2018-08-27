import os
from nipype import logging
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl

from CPAC.utils import function
from CPAC.network_centrality.utils import merge_lists
from CPAC.network_centrality import (
    create_resting_state_graphs,
    get_cent_zscore
)

logger = logging.getLogger('workflow')

# TODO ASH redo?
# Check for the existence of AFNI 3dDegreeCentrality/LFCD binaries
import subprocess
try:
    ret_code = subprocess.check_call(['which', '3dDegreeCentrality'],
                                     stdout=open(os.devnull, 'wb'))
    if ret_code == 0:
        afni_centrality_found = True
except subprocess.CalledProcessError as exc:
    afni_centrality_found = False

try:
    ret_code = subprocess.check_call(['which', '3dLFCD'],
                                     stdout=open(os.devnull, 'wb'))
    if ret_code == 0:
        afni_lfcd_found = True
except subprocess.CalledProcessError as exc:
    afni_lfcd_found = False


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

        # Get nipype  node and out file of the func mni img
        node, out_file = strat.get_node_from_resource_pool(
            'functional_to_standard'
        )

        # Resample the input functional file to template(roi/mask)
        workflow.connect(node, out_file,
                         resample_functional_to_template, 'in_file')

        resample_functional_to_template.inputs.reference = \
            c.templateSpecificationFile

        # Init merge node for appending method output lists to one another
        merge_node = pe.Node(function.Function(input_names=['deg_list',
                                                            'eig_list',
                                                            'lfcd_list'],
                                               output_names=['merged_list'],
                                               function=merge_lists,
                                               as_module=True),
                             name='merge_node_%d' % num_strat)

        # Function to connect the CPAC centrality python workflow
        # into pipeline

        # Degree/eigen check
        if afni_centrality_found:
            if True in c.degWeightOptions:
                connect_afni_centrality_wf(
                    workflow, c, strat, num_strat,
                    resample_functional_to_template, merge_node,
                    'degree',
                    c.degCorrelationThresholdOption,
                    c.degCorrelationThreshold
                )
            if True in c.eigWeightOptions:
                connect_afni_centrality_wf(
                    workflow, c, strat, num_strat,
                    resample_functional_to_template, merge_node,
                    'eigenvector',
                    c.eigCorrelationThresholdOption,
                    c.eigCorrelationThreshold
                )
        # Otherwise run the CPAC python workflow
        else:
            # If we're calculating degree centrality
            if True in c.degWeightOptions:
                connect_centrality_workflow(
                    workflow, c, strat, num_strat,
                    resample_functional_to_template, merge_node,
                    'degree',
                    c.degCorrelationThresholdOption,
                    c.degCorrelationThreshold,
                    c.degWeightOptions,
                    'deg_list'
                )
            # If we're calculating eigenvector centrality
            if True in c.eigWeightOptions:
                connect_centrality_workflow(
                    workflow, c, strat, num_strat,
                    resample_functional_to_template, merge_node,
                    'eigenvector',
                    c.eigCorrelationThresholdOption,
                    c.eigCorrelationThreshold,
                    c.eigWeightOptions,
                    'eig_list'
                )
        # LFCD check
        if afni_lfcd_found:
            # If we're calculating lFCD
            if True in c.lfcdWeightOptions:
                connect_afni_centrality_wf(
                    workflow, c, strat, num_strat,
                    resample_functional_to_template, merge_node,
                    'lfcd',
                    c.lfcdCorrelationThresholdOption,
                    c.lfcdCorrelationThreshold
                )
        # Otherwise run the CPAC python workflow
        else:
            # If we're calculating lFCD
            if True in c.lfcdWeightOptions:
                connect_centrality_workflow(
                    workflow, c, strat, num_strat,
                    resample_functional_to_template, merge_node,
                    'lfcd',
                    c.lfcdCorrelationThresholdOption,
                    c.lfcdCorrelationThreshold,
                    c.lfcdWeightOptions,
                    'lfcd_list'
                )

        # Update resource pool with centrality outputs
        strat.update_resource_pool({
            'centrality': (merge_node, 'merged_list')
        })

        if 0 in c.runNetworkCentrality:
            strat = strat.fork()
            strategies += [strat]

    return strategies


def connect_centrality_workflow(workflow, c, strat, num_strat,
                                resample_functional_to_template, merge_node,
                                methodOption, thresholdOption,
                                threshold, weightOptions, mList):

    # Create centrality workflow
    network_centrality = \
        create_resting_state_graphs(
            wf_name='network_centrality_%d_%s' % (num_strat, methodOption),
            allocated_memory=c.memoryAllocatedForDegreeCentrality
        )

    # Connect resampled (to template/mask resolution)
    # functional_mni to inputspec
    workflow.connect(resample_functional_to_template, 'out_file',
                     network_centrality, 'inputspec.in_file')

    network_centrality.inputs.inputspec.set(
        # Subject mask/parcellation image
        template=c.templateSpecificationFile,
        # Give which method we're doing
        method_option=methodOption,
        # Type of threshold
        threshold_option=thresholdOption,
        # Connect threshold value (float)
        threshold=threshold
    )

    # Merge output with others via merge_node connection
    workflow.connect(network_centrality, 'outputspec.centrality_outputs',
                     merge_node, mList)

    # Append this as a strategy
    strat.append_name(network_centrality.name)


# Function to connect the afni 3dDegreeCentrality workflow
# into pipeline
def connect_afni_centrality_wf(workflow, c, strat, num_strat,
                               resample_functional_to_template, merge_node,
                               method_option, threshold_option,
                               threshold):

    # Import packages
    from CPAC.network_centrality.afni_network_centrality \
        import create_afni_centrality_wf
    import CPAC.network_centrality.utils as cent_utils

    # Init variables
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
        cent_utils.check_centrality_params(method_option,
                                           threshold_option,
                                           threshold)

    # Change sparsity thresholding to % to work with afni
    if threshold_option == 'sparsity':
        threshold = threshold * 100

    # Init the workflow
    afni_centrality_wf = \
        create_afni_centrality_wf(wf_name, method_option,
                                  threshold_option,
                                  threshold, num_threads, memory)

    # Connect pipeline resources to workflow
    workflow.connect(resample_functional_to_template, 'out_file',
                     afni_centrality_wf, 'inputspec.in_file')

    # Mask
    afni_centrality_wf.inputs.inputspec.template = \
        c.templateSpecificationFile

    # Connect outputs to merge node
    workflow.connect(afni_centrality_wf,
                     'outputspec.outfile_list',
                     merge_node,
                     out_list)
