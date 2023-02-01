# -*- coding: utf-8 -*-

import glob
from copy import deepcopy
from CPAC.pipeline import nipype_pipeline_engine as pe
from nipype.interfaces.utility import Function
from nipype.pipeline.engine.utils import generate_expanded_graph
import networkx as nx

from indi_aws import fetch_creds

from CPAC.utils.datasource import (
    create_check_for_s3_node,
)

def expand_workflow(wf):
    return generate_expanded_graph(deepcopy(wf._create_flat_graph()))

def is_datasink(n):
    return type(n).__name__ == 'Node' and type(n.interface).__name__ == 'DataSink'

def compute_datasink_dirs(graph, datasink, output_dir=None, container=None):
    directories = {}
    for inp in graph.in_edges(datasink):
        src, _ = inp
        for edge in graph.get_edge_data(*inp)['connect']:
            _, derivative_name = edge

            datasink_output_dir = datasink.interface.inputs.base_directory
            if output_dir is not None:
                datasink_output_dir = output_dir

            datasink_container = datasink.interface.inputs.container
            if container is not None:
                datasink_container = container

            # Look if there is an output in this datasink directory

            iterables = datasink.parameterization
            path = '/'.join(['', derivative_name] + iterables)
            path = datasink.interface._substitute(path)[1:]
            path = '/'.join([datasink_output_dir, datasink_container, path])

            directories[(src, derivative_name)] = path

    return directories

def list_files(path, s3_creds_path=None):
    if path.startswith('s3://'):
        pieces = path[5:].split('/')
        bucket_name, path = pieces[0], '/'.join(pieces[1:])
        bucket = fetch_creds.return_bucket(s3_creds_path, bucket_name)
        return [
            's3://%s/%s' % (bucket, obj['Key'])
            for obj in bucket.objects.filter(Prefix=path)
        ]
    else:
        return list(glob.glob(path + '/*'))

def the_trimmer(wf, output_dir=None, container=None, s3_creds_path=None):
    """
    The trimmer: trimming your workflow based on its datasinks.

    The assumption here is that all your desired outputs will be placed in an
    output directory by using a DataSink node.

    By analysing a specific output directory, and comparing what is in it with
    the DataSinks of a workflow, it is possible to audit which Datasinks have
    already outputted to the output directory. If the DataSink has already
    fulfilled its function, we infer that previous nodes also had, since they
    are prerequisites for the DataSink to run. This is the simplest case in which
    we can prune nodes (i.e. do not execute them), optimizing the execution
    time.

    A brief syntax note:
    [node] is a node
    → is a connection, disconsidering which field/attribute
    →(field)→ indicates a connection considering the field
    ✓ is a datasink with an existing file
    ❌ is a datasing witn an missing file

    E.g.

    [node1] → [node2] → [node3] → [datasink to file.txt ✓]

    since file.txt already exist, there is no need to execute the [node3].
    Since [node2] only outputs to [node3], and [node3] will not be executed,
    we can prune [node2]. Same for [node1]. In this case, our workflow will
    not have any more nodes to run.

    There are more complex cases:

    1) A node outputs for several nodes, and some of their results are not
    in the output directory.

    [node1] → [node2] → [node3] → [datasink to file1.txt ✔]
            ↳ [node4] → [datasink to file2.txt ❌]

    for this case, we cannot prune [node1], since its output is used in
    another branch, for [node4], that is not cached. After trimming,
    the remaining workflow is:

    [node1] 
            ↳ [node4] → [datasink to file2.txt ❌]

    2) The node has several outputs, and an uncached branch down the
    graph requires one of its outputs.

    [registration] →(warped image)→ [datasink to warped.nii.gz ✔]
                   ↳(transforms)→ [apply transforms] → [datasink to func_warped.nii.gz ❌]
                   [functional] ↗

    given func_warped.nii.gz is not cached, we need to perform "apply transforms", that
    requires the transforms from the [registration] node. In this case, even that warped.nii.gz
    is cached, we will reexecute the [registration] again to get the transforms. After trimming,
    the remaining workflow is:

    [registration] 
                   ↳(transforms)→ [apply transforms] → [datasink to func_warped.nii.gz ❌]
                   [functional] ↗


    For this implementation, we disregard MapNodes, as their outputs is harder to check.

    Iterables are considered in the implementation by expanding the original workflow
    into what is called an execution graph, creating a node for each iterable value.

    Parameters
    ----------
    wf : Workflow
        A Nipype workflow to be pruned.

    output_dir : Path
        The directory in which the outputs are stored. If not provided, value is inferred
        from the DataSink nodes.

    container : Path
        The subdirectory from the output_dir in which the output are stored. If not provided,
        value is inferred from the DataSink nodes.
    
    s3_creds_path : Path
        Path to S3 credentials, in case output_dir is in a S3 bucket.
    
    Returns
    -------
    wf_new : Workflow
        Prunned workflow

    (replacement_mapping, deletions): (Dict, List)
        
        replacement_mapping contains the nodes replaces with input nodes, pointing to
        files from the output_dir

        deletions contains the nodes removed from the workflow, as they do not need to be
        executed
    
    """

    # Expand graph, to flatten out sub-workflows and iterables
    execgraph = expand_workflow(wf)

    replacements = {}
    deletions = []
    
    # Check out for datasinks (i.e. the ones who throws things at the output dir)
    datasinks = [
        n for n in execgraph.nodes()
        if is_datasink(n)
    ]
    for datasink in datasinks:

        for (src, derivative_name), path in \
            compute_datasink_dirs(execgraph,
                                  datasink,
                                  output_dir=output_dir,
                                  container=container).items():

            files = list_files(path, s3_creds_path=s3_creds_path)
            if len(files) == 1:  # Ignore multi-file nodes
                if src not in replacements:
                    replacements[src] = {}

                replacements[src][src_field] = files[0]
        
        # if the replacements have all the fields from the datasink, datasink 
        # can be deleted (we do not want to output again the same file :))
        if all(
            any(
                field in replacements.get(src, {})
                for field, _ in execgraph.get_edge_data(src, dst)['connect']
            )
            for src, dst in execgraph.in_edges(datasink)
        ):
            deletions += [datasink]

    # Remove from replacement list the nodes that gives other output
    # for other nodes, since it seems like not all fields are cached
    for node, cached_fields in replacements.items():
        for edge in execgraph.out_edges(node):
            if any(
                src_field not in cached_fields
                for src_field, _ in execgraph.get_edge_data(*edge)['connect']
            ):
                del replacements[node]
                break
            
    # Delete them! It also removes the edges, and recursively delete nodes
    # before rationalizing about replacements
    for node in reversed(nx.topological_sort(execgraph)):
        if node in deletions:
            execgraph.remove_node(node)
        if is_datasink(node):
            continue
        
        if len(execgraph.out_edges(node)) == 0:
            deletions += [node]
            if node in replacements:
                del replacements[node]
            execgraph.remove_node(node)
                    
    # And now we replace the cached with a data input node, from the 
    # output directory.
    replacement_mapping = {}
    for replacement, cached_files in replacements.items():
        out_edges_data = execgraph.edge[replacement]
        
        # Get this cached node, and replace all the out-connections
        # from this node with a data input node
        out_edges = execgraph.successors(replacement)
        if out_edges:
            for to_node in out_edges:
                for from_field, to_field in out_edges_data[to_node]['connect']:
                    
                    # Reuse the data input node for this field
                    if replacement not in replacement_mapping:
                        replacement_mapping[replacement] = {}
                    if from_field not in replacement_mapping[replacement]:
                
                        new_node = create_check_for_s3_node(
                            name='%s_%s_triminput' % (replacement.name, from_field),
                            file_path=cached_files[from_field],
                            img_type='other',
                            creds_path=s3_creds_path,
                            dl_dir=None
                        )
                        new_node._hierarchy = deepcopy(replacement._hierarchy)

                        execgraph.add_node(new_node)
                        replacement_mapping[replacement][from_field] = new_node
                
                    # Connect the new data input node to the node
                    # it was connected
                    execgraph.add_edge(
                        replacement_mapping[replacement][from_field],
                        to_node,
                        connect=[('local_path', to_field)]
                    )

        execgraph.remove_node(replacement)
        
    # Second round of backtrack deletion, affected by replacements
    for node in reversed(nx.topological_sort(execgraph)):
        if is_datasink(node):
            continue

        if len(execgraph.out_edges(node)) == 0:
            deletions += [node]
            if node in replacements:
                del replacements[node]
            execgraph.remove_node(node)
        
    wf_new = wf.clone(wf.name + '_trimmed')
    wf_new.name = wf.name
    wf_new._graph = execgraph
    
    return wf_new, (replacement_mapping, deletions)
