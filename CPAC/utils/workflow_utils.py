import pickle
from typing import Any

from nipype import __version__ as version_nipype

from CPAC import __version__ as version_cpac

VERSION_WORKFLOW = 1
"""Placeholder for version checking."""


def dump_workflow(filename: str, workflow: Any) -> None:
    """
    Serialize and save workflow object to a file.

    Parameters
    ----------
    filename : Filename to save to.
    workflow : Workflow object.
    """

    obj = {
        'version': {
            'workflow': VERSION_WORKFLOW,
            'cpac': version_cpac,
            'nipype': version_nipype
        },
        'workflow': workflow
    }
    with open(filename, 'wb') as handle:
        pickle.dump(obj, file=handle, protocol=pickle.HIGHEST_PROTOCOL)
