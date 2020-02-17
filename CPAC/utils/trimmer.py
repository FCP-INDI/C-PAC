import glob
from copy import deepcopy
import nipype.pipeline.engine as pe
from nipype.interfaces.utility import Function
from nipype.pipeline.engine.utils import generate_expanded_graph

from indi_aws import fetch_creds

from CPAC.utils.datasource import (
    create_check_for_s3_node,
)


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

    # Expand graph, to flatten out sub-workflows and iterables
    execgraph = generate_expanded_graph(deepcopy(wf._create_flat_graph()))

    replacements = {}
    deletions = []
    
    # Check out for datasinks (i.e. the ones who throws things at the output dir)
    execnodes = [
        n for n in execgraph.nodes()
        if type(n).__name__ == 'Node' and type(n.interface).__name__ == 'DataSink'
    ]
    for datasink in execnodes:

        # For each input node (DataSink may have several, but C-PAC usually uses one)
        for inp in execgraph.in_edges(datasink):

            src, _ = inp

            # ... and it can have several fields per node
            for edge in execgraph.get_edge_data(*inp)['connect']:

                src_field, derivative_name = edge

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

                # TODO support S3
                files = list_files(path, s3_creds_path=None)
                if len(files) == 1:  # Ignore multi-file nodes
                    if src not in replacements:
                        replacements[src] = {}

                    replacements[src][src_field] = files[0]
        
        # if the replacements have all the fields from the datasink, datasink can be deleted
        # (we do not want to output again the same file :))
        if all(
            any(
                field in replacements.get(src, {})
                for field, _ in execgraph.get_edge_data(src, dst)['connect']
            )
            for src, dst in execgraph.in_edges(datasink)
        ):
            deletions += [datasink]

    # Delete nodes that gives other output for other nodes, since it seems like not all fields are cached
    for node, cached_fields in replacements.items():
        for edge in execgraph.out_edges(node):
            if any(
                src_field not in cached_fields
                for src_field, _ in execgraph.get_edge_data(*edge)['connect']
            ):
                del replacements[node]
                break

    def recurse_predecessors(execgraph, cached_node):

        # for each predecessor of a cached node
        for node, _ in execgraph.in_edges(cached_node):

            # if it only points to the cached node
            if any(outnode != cached_node for _, outnode in execgraph.out_edges(node)):
                continue

            # we dont need to run it
            yield node
            # ... and we can check if there are other nodes to remove
            for n in recurse_predecessors(execgraph, node):
                yield n

    # Recurse nodes that generates the replacement inputs and delete them
    # whenever possible i.e. there is no other connections down the line
    for node, cached_fields in replacements.items():
        deletions += recurse_predecessors(execgraph, node)
    
    # Delete them! It also removes the edges
    for deletion in deletions:
        execgraph.remove_node(deletion)
        
    replacement_mapping = {}
    
    _input_nodes = []
    for replacement, cached_files in replacements.items():
        
        out_edges = execgraph.successors(replacement)
        if out_edges:
            out_edges_data = execgraph.edge[replacement]

            for to_node in out_edges:
                
                for from_field, to_field in out_edges_data[to_node]['connect']:
                    
                    if replacement not in replacement_mapping:
                        replacement_mapping[replacement] = {}
                    
                    # Reuse input fields
                    if from_field not in replacement_mapping[replacement]:
                
                        new_node = create_check_for_s3_node(
                            name='%s_%s_input' % (replacement.name, from_field),
                            file_path=cached_files[from_field],
                            img_type='other',
                            creds_path=s3_creds_path,
                            dl_dir=None
                        )
                        new_node._hierarchy = deepcopy(replacement._hierarchy)

                        execgraph.add_node(new_node)
                        replacement_mapping[replacement][from_field] = new_node
                
                        _input_nodes += [new_node]

                    execgraph.add_edge(
                        replacement_mapping[replacement][from_field],
                        to_node,
                        connect=[('file', to_field)]
                    )

        execgraph.remove_node(replacement)
        
    # Double check to assess if input nodes are really required
    for node in _input_nodes:
        out_edges = execgraph.successors(node)
        if not out_edges:
            execgraph.remove_node(node)
        
    wf_new = wf.clone(wf.name + '_trimmed')
    wf_new._graph = execgraph
    
    return wf_new, (replacement_mapping, deletions)
