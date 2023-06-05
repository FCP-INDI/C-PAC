import dataclasses
import json
from dataclasses import dataclass
from datetime import datetime, date
from os import PathLike
from typing import List, Dict, Union

import networkx

from .core import workflow_container, WorkflowRaw, NodeRaw


def _object_as_strdict(obj: object) -> Dict[str, str]:
    """Extracts and converts all fields of an object to a {str: str} dict."""
    return {str(k): str(v) for k, v in obj.__dict__.items()}


def _workflow_get_graph(wf: WorkflowRaw) -> networkx.DiGraph:
    """Get graph from workflow instance."""
    return wf._graph  # noqa


@dataclass
class EdgeData:
    """Data class for serializing workflow edges."""

    fullname_origin: str
    fullname_target: str

    @classmethod
    def from_obj(cls, obj):
        return cls(
            fullname_origin=obj[0].fullname,
            fullname_target=obj[1].fullname
        )


@dataclass
class NodeData:
    """Data class for serializing workflow nodes."""

    name: str
    fullname: str
    type: str
    repr: str

    inputs: dict
    outputs: dict

    nodes: List['NodeData'] = dataclasses.field(default_factory=lambda: [])
    edges: List['EdgeData'] = dataclasses.field(default_factory=lambda: [])

    @classmethod
    def from_obj(cls, obj):
        node_data = cls(
            name=obj.name,  # There is name, fullname and itername
            fullname=obj.fullname,
            type='node',
            repr=str(obj),
            inputs=_object_as_strdict(obj.inputs),
            outputs=_object_as_strdict(obj.outputs),
            nodes=[],
            edges=[]
        )

        if isinstance(obj, NodeRaw):
            return node_data

        if isinstance(obj, WorkflowRaw):
            node_data.type = 'workflow'
            graph = _workflow_get_graph(obj)

            for child_node in graph.nodes:
                node_data_child = cls.from_obj(child_node)
                node_data.nodes.append(node_data_child)

            for child_edgle in graph.edges:
                edge_data_child = EdgeData.from_obj(child_edgle)
                node_data.edges.append(edge_data_child)

            return node_data

        raise TypeError(f'Unknown Node type found in graph: {type(obj)}')


@dataclass
class WorkflowJSONMeta:
    """Data class for meta information."""
    pipeline_name: str
    stage: str
    time: datetime = dataclasses.field(default_factory=lambda: datetime.now().astimezone())

    def filename(self) -> str:
        """Generate filename from fields"""
        timestamp = self.time.strftime("%Y-%m-%d_%H-%M-%S")
        return f'workflow_{self.stage}_{timestamp}_{self.pipeline_name}.json'


class WorkflowJSONEncoder(json.JSONEncoder):
    """Custom JSON encoder. Unpacks dataclasses to dicts."""

    def default(self, o):
        if dataclasses.is_dataclass(o):
            return dataclasses.asdict(o)
        if isinstance(o, (datetime, date)):
            return o.isoformat()
        return super().default(o)


def save_workflow_json(
        filename: Union[str, PathLike],
        workflow: WorkflowRaw,
        meta: WorkflowJSONMeta
) -> None:
    """
    Serialize and save workflow object to a file.

    Parameters
    ----------
    filename : Filename to save to.
    workflow : Workflow object.
    meta : Meta information.
    """

    node_data = NodeData.from_obj(workflow)
    obj = workflow_container(workflow=node_data, meta=meta)
    with open(filename, 'w', encoding='utf-8') as file:
        json.dump(obj, file, indent=2, cls=WorkflowJSONEncoder)
