import os
import glob
import json
import datetime
import threading
import SocketServer


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
        'hash': node.inputs.get_hashval()[1],
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


class LoggingRequestHandler(SocketServer.BaseRequestHandler):

    def handle(self):

        tree = {}

        logs = glob.glob(
            os.path.join(
                self.server.logging_dir,
                "pipeline_" + self.server.pipeline_name,
                "*"
            )
        )

        for log in logs:
            subject = log.split('/')[-1]
            tree[subject] = {}
            
            callback_file = os.path.join(log, "callback.log")

            if not os.path.exists(callback_file):
                continue

            with open(callback_file, 'rb') as lf:
                for l in lf.readlines():
                    l = l.strip()
                    try:
                        node = json.loads(l)
                        if node["id"] not in tree[subject]:
                            tree[subject][node["id"]] = {
                                "hash": node["hash"]
                            }
                            if "start" in node and "finish" in node:
                                tree[subject][node["id"]]["start"] = node["start"]
                                tree[subject][node["id"]]["finish"] = node["finish"]

                        else:
                            if "start" in node and "finish" in node:
                                if tree[subject][node["id"]]["hash"] == node["hash"]:
                                    tree[subject][node["id"]]["cached"] = {
                                        "start": node["start"],
                                        "finish": node["finish"],
                                    }

                                # pipeline was changed, and we have a new hash
                                else:
                                    tree[subject][node["id"]]["start"] = node["start"]
                                    tree[subject][node["id"]]["finish"] = node["finish"]

                    except:
                        break

        self.request.sendall(json.dumps(tree))


class LoggingHTTPServer(SocketServer.ThreadingTCPServer, object):
    
    def __init__(self, pipeline_name, logging_dir='', host='', port=8080, request=LoggingRequestHandler):
        super(LoggingHTTPServer, self).__init__((host, port), request)

        if not logging_dir:
            logging_dir = os.getcwd()

        self.logging_dir = logging_dir
        self.pipeline_name = pipeline_name


def monitor_server(pipeline_name, logging_dir, host='0.0.0.0', port=8080):
    httpd = LoggingHTTPServer(pipeline_name, logging_dir, host, port, LoggingRequestHandler)

    server_thread = threading.Thread(target=httpd.serve_forever)
    server_thread.isDaemon = True
    server_thread.start()

    return server_thread