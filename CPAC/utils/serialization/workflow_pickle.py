import pickle
from typing import Any

from .core import workflow_container


def save_workflow_pickle(filename: str, workflow: Any) -> None:
    """
    Serialize and save workflow object to a file.

    Parameters
    ----------
    filename : Filename to save to.
    workflow : Workflow object.
    """

    obj = workflow_container(workflow)
    with open(filename, 'wb') as handle:
        pickle.dump(workflow_container(obj), file=handle, protocol=pickle.HIGHEST_PROTOCOL)
