import os
import glob
import json
import time
import logging
import datetime
import threading
import queue
import asyncio
import websockets

import networkx as nx
import nipype.pipeline.engine as pe
import nipype.pipeline.engine.nodes as nodes

is_monitoring = threading.Lock()
monitor_nodes_queue = queue.Queue()
monitor_wait_lock = threading.Lock()

monitor_message = lambda data: f'{{"time": {time.time()}, "message": {data}}}'

def recurse_nodes(workflow, prefix=''):
    """
    Function to traverse the workflow, yielding its nodes

    Parameters
    ----------
    workflow : nipype.pipeline.engine.Workflow
        the workflow

    prefix : string
        custom prefix to prepend to the node name

    """
    for node in nx.topological_sort(workflow._graph):
        if isinstance(node, pe.Workflow):
            for subnode in recurse_nodes(node, prefix + workflow.name + '.'):
                yield subnode
        else:
            yield {
                "id": prefix + workflow.name + '.' + node.name,
                "hash": node.inputs.get_hashval()[1],
            }


def log_nodes_initial(workflow):
    """
    Function to record all the existing nodes in the workflow

    Parameters
    ----------
    workflow : nipype.pipeline.engine.Workflow
        the workflow

    """
    logger = logging.getLogger('callback')
    for node in recurse_nodes(workflow):
        data = json.dumps(node)
        logger.debug(data)
        if is_monitoring.locked():
            monitor_nodes_queue.put(monitor_message(data))


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

    if isinstance(node, nodes.MapNode):
        return

    logger = logging.getLogger('callback')

    try:
        runtime = node.result.runtime
    except:
        runtime = None

    status_dict = {
        'id': str(node),
        'hash': node.inputs.get_hashval()[1],
        'start': getattr(runtime, 'startTime') if runtime else None,
        'end': getattr(runtime, 'endTime') if runtime else None,
        'runtime_threads': getattr(runtime, 'cpu_percent', 'N/A') if runtime else 'N/A',
        'runtime_memory_gb': getattr(runtime, 'mem_peak_gb', 'N/A') if runtime else 'N/A',
        'estimated_memory_gb': node.mem_gb,
        'num_threads': node.n_procs,
    }

    if status_dict['start'] is None or status_dict['end'] is None:
        status_dict['error'] = True

    data = json.dumps(status_dict)
    logger.debug(data)

    if is_monitoring.locked():
        monitor_nodes_queue.put(monitor_message(data))


async def send_logs(websocket, path):
    """
    Function executed when a new websocket connection is opened.

    It will iterate through the nodes reporting queue, and send it to
    the websocked.
    
    Parameters
    ----------
    websocket : Websocket
        The websocket object

    path : string
        Path for the websocket (the accepted is ws://0.0.0.0:8080/log)
    
    """
    if path != '/log':
        return

    monitor_wait_lock.release()
    while True:
        item = monitor_nodes_queue.get()
        await websocket.send(item)
        monitor_nodes_queue.task_done()

def ws(loop, host='0.0.0.0', port=8080):
    """
    Function to start the monitoring server.
    
    Parameters
    ----------
    loop : asyncio.Loop
        The loop used by the server

    host : string
        Address or IP that the websocket server binds to.

    port : integer
        Port for the websocket server.
    
    """
    asyncio.set_event_loop(loop)
    start_server = websockets.serve(send_logs, "0.0.0.0", 8080)
    loop.run_until_complete(start_server)
    loop.run_forever()


def monitor_server(host='0.0.0.0', port=8080, wait=False):
    """
    Function to start the monitoring server thread.
    It sets the `monitor_wait_lock` so the pipeline waits for the websocket
    to connect. This way, no reported data is lost.
    
    Parameters
    ----------
    host : string
        Address or IP that the websocket server binds to.

    port : integer
        Port for the websocket server.

    wait : boolean
        Only sets the lock if required.
    
    Returns
    -------
    server_thread : threading.Thread
        The thread the server is running

    loop: asyncio.Loop
        The loop used by the server
    
    """
    if wait:
        monitor_wait_lock.acquire()
    
    is_monitoring.acquire()

    loop = asyncio.new_event_loop()
    server_thread = threading.Thread(target=ws, args=[loop, host, port])
    server_thread.isDaemon = True
    server_thread.start()
    return server_thread, loop
