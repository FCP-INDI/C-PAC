
def log_nodes_cb(node, status):
    """Function to record node run statistics to a log file as json
    dictionaries

    Parameters
    ----------
    node : nipype.pipeline.engine.Node
        the node being logged
    status : string
        acceptable values are 'start', 'end'; otherwise it is
        considered and error

    Returns
    -------
    None
        this function does not return any values, it logs the node
        status info to the callback logger
    """

    if status != 'end':
        return

    # Import packages
    import json
    import logging
    import nipype.pipeline.engine.nodes as nodes

    logger = logging.getLogger('callback')

    if isinstance(node, nodes.MapNode):
        return

    runtime = node.result.runtime

    status_dict = {
        'name': node.name,
        'id': str(node),
        'start': getattr(runtime, 'startTime'),
        'finish': getattr(runtime, 'endTime'),
        'duration': getattr(runtime, 'duration'),
        'runtime_threads': getattr(runtime, 'cpu_percent', 'N/A'),
        'runtime_memory_gb': getattr(runtime, 'mem_peak_gb', 'N/A'),
        'estimated_memory_gb': node.mem_gb,
        'num_threads': node.n_procs,
    }

    if status_dict['start'] is None or status_dict['finish'] is None:
        status_dict['error'] = True

    logger.debug(json.dumps(status_dict))