'''Common functions for overriding Nipype plugins'''
import faulthandler
import sys
from nipype.pipeline.plugins.base import format_exception, logger

logfile = logger.handlers[0].baseFilename if (
    len(logger.handlers) and hasattr(logger.handlers[0].baseFilename)
) else sys.stderr


# Run node
def run_node(node, updatehash, taskid):
    """Function to execute node.run(), catch and log any errors and
    return the result dictionary
    Parameters
    ----------
    node : nipype Node instance
        the node to run
    updatehash : boolean
        flag for updating hash
    taskid : int
        an identifier for this task
    Returns
    -------
    result : dictionary
        dictionary containing the node runtime results and stats
    """

    # Init variables
    result = dict(result=None, traceback=None, taskid=taskid)

    faulthandler.dump_traceback_later(timeout=10, file=logfile, exit=True)
    # Try and execute the node via node.run()
    try:
        result["result"] = node.run(updatehash=updatehash)
    except:  # noqa: E722, intendedly catch all here
        result["traceback"] = format_exception(*sys.exc_info())
        result["result"] = node.result

    # Return the result dictionary
    return result
