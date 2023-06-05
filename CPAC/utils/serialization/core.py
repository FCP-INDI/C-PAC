from typing import Any

from nipype import __version__ as version_nipype
from nipype.pipeline.engine import Workflow, Node

from CPAC import __version__ as version_cpac

WorkflowRaw = Workflow
"""Alias for internal use. CPAC derives: CPAC.pipeline.nipype_pipeline_engine.engine.Workflow"""

NodeRaw = Node
"""Alias for internal use. CPAC derives: CPAC.pipeline.nipype_pipeline_engine.engine.Node"""

VERSION_WORKFLOW = 1
"""Placeholder for version checking. Increase when file structure changes drastically."""


def workflow_container(workflow: Any, meta: Any = None) -> dict:
    """
    Construct a container dictionary with some version information to allow version checks when reading.

    Parameters
    ----------
    workflow : Workflow object.
    meta : Meta information.

    Returns
    -------

    """
    return {
        'version': {
            'workflow': VERSION_WORKFLOW,
            'cpac': version_cpac,
            'nipype': version_nipype
        },
        'meta': {} if meta is None else meta,
        'workflow': workflow,
    }
