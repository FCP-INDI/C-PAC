import dataclasses
import json
import pathlib as pl
from dataclasses import dataclass
from datetime import datetime, date
from os import PathLike
from typing import List, Dict, Union, Optional

import networkx
from traits.has_traits import HasTraits
from traits.trait_base import _Undefined  # noqa

from .core import workflow_container, WorkflowRaw, NodeRaw
from .nipype_report_parser import read_report_rst, ReportSection
from .. import Configuration


def _object_as_strdict(obj: object) -> Dict[str, str]:
    """Extracts and converts all fields of an object to a {str: str} dict."""
    if isinstance(obj, dict):
        obj_dict = obj
    elif hasattr(obj, '__dict__'):
        obj_dict = obj.__dict__
    else:
        return {'value': str(obj)}
    return {str(k): str(v) for k, v in obj_dict.items()}


def _serialize_inout(obj: object):
    if isinstance(obj, dict):
        if 'json_data' in obj:
            obj_clone = obj.copy()
            obj_clone['json_data'] = '[truncated]'
            obj = obj_clone
        return {str(k): _serialize_inout(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_serialize_inout(i) for i in obj]
    if isinstance(obj, tuple):
        return (_serialize_inout(i) for i in obj)
    elif isinstance(obj, HasTraits):
        return _serialize_inout(obj.trait_get())
    elif isinstance(obj, (int, float, str)):
        return obj
    elif isinstance(obj, _Undefined):
        return '[Undefined]'
    elif isinstance(obj, Configuration):
        return '[C-PAC config]'
    elif obj is None:
        return None
    else:
        print(type(obj))
        return str(obj)


def _workflow_get_graph(wf: WorkflowRaw) -> networkx.DiGraph:
    """Get graph from workflow instance."""
    return wf._graph  # noqa


def _node_get_wd_path(node: NodeRaw) -> pl.Path:
    """
    Get working directory path for a given node.
    """
    return pl.Path(node.output_dir())  # noqa


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
class NodeResultData:
    result_inputs: dict
    result_outputs: dict
    runtime_info: dict
    wd_path: str
    read_success: bool


@dataclass
class NodeData:
    """Data class for serializing workflow nodes."""

    name: str
    fullname: str
    type: str
    repr: str

    inputs: dict
    outputs: dict

    result: Optional[NodeResultData] = None

    nodes: List['NodeData'] = dataclasses.field(default_factory=lambda: [])
    edges: List['EdgeData'] = dataclasses.field(default_factory=lambda: [])

    @classmethod
    def from_obj(cls, obj, serialze_postex: bool = True):
        node_data = cls(
            name=obj.name,  # There is name, fullname and itername
            fullname=obj.fullname,
            type='node',
            repr=str(obj),
            inputs=_object_as_strdict(_serialize_inout(obj.inputs)),
            outputs=_object_as_strdict(_serialize_inout(obj.outputs))
        )

        if isinstance(obj, NodeRaw):
            if serialze_postex:
                node_data.result = NodeResultData(
                    result_inputs={},
                    result_outputs={},
                    runtime_info={},
                    wd_path=_node_get_wd_path(obj).as_posix(),
                    read_success=True
                )
                report_path = _node_get_wd_path(obj) / '_report' / 'report.rst'
                if report_path.exists():
                    report_data = read_report_rst(report_path)
                    node_data.result_inputs = report_data.get(ReportSection.EXECUTION_INPUTS, {})
                    node_data.result_outputs = report_data.get(ReportSection.EXECUTION_OUTPUTS, {})
                    node_data.runtime_info = report_data.get(ReportSection.EXECUTION_INFO, {})
                else:
                    node_data.read_success = False
                    print(f'Report RST does not exist: "{report_path}"')
                    print(f' > Node.base_dir: "{obj.base_dir}"')
                    print(f' > Node.output_dir: "{obj.output_dir()}"')

            return node_data

        if isinstance(obj, WorkflowRaw):
            node_data.type = 'workflow'
            graph = _workflow_get_graph(obj)
        elif isinstance(obj, networkx.DiGraph):
            graph = obj
        else:
            raise TypeError(f'Unknown Node type found in graph: {type(obj)}')

        for child_node in graph.nodes:
            node_data_child = cls.from_obj(child_node, serialze_postex=serialze_postex)
            node_data.nodes.append(node_data_child)

        for child_edge in graph.edges:
            edge_data_child = EdgeData.from_obj(child_edge)
            node_data.edges.append(edge_data_child)

        return node_data


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

    node_data = NodeData.from_obj(workflow, serialze_postex=meta.stage == 'post')
    obj = workflow_container(workflow=node_data, meta=meta)
    with open(filename, 'w', encoding='utf-8') as file:
        json.dump(obj, file, indent=2, cls=WorkflowJSONEncoder)
