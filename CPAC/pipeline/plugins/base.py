'''Common functions for overriding Nipype plugins'''
import faulthandler
import sys
from nipype.pipeline.plugins.base import format_exception, \
    logger  # noqa F401


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

    # Try and execute the node via node.run()
    if not node.skip_timeout:
        # 86400 seconds is 24 hours
        faulthandler.dump_traceback_later(timeout=86400, exit=True)
    try:
        result["result"] = node.run(updatehash=updatehash)
    except:  # noqa: E722, intendedly catch all here
        result["traceback"] = format_exception(*sys.exc_info())
        result["result"] = node.result

    # Return the result dictionary
    return result
