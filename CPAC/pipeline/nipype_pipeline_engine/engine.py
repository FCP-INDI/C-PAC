# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Supports just-in-time dynamic memory allocation
#     * Skips doctests that require files that we haven't copied over
#     * Applies a random seed
#     * Supports overriding memory estimates via a log file and a buffer
#     * Adds quotation marks around strings in dotfiles

# ORIGINAL WORK'S ATTRIBUTION NOTICE:
#     Copyright (c) 2009-2016, Nipype developers

#     Licensed under the Apache License, Version 2.0 (the "License");
#     you may not use this file except in compliance with the License.
#     You may obtain a copy of the License at

#         http://www.apache.org/licenses/LICENSE-2.0

#     Unless required by applicable law or agreed to in writing, software
#     distributed under the License is distributed on an "AS IS" BASIS,
#     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#     See the License for the specific language governing permissions and
#     limitations under the License.

#     Prior to release 0.12, Nipype was licensed under a BSD license.

# Modifications Copyright (C) 2022-2023 C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
'''Module to import Nipype Pipeline engine and override some Classes.
See https://fcp-indi.github.io/docs/developer/nodes
for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.pipeline.engine.html
for Nipype's documentation.'''  # noqa: E501  # pylint: disable=line-too-long
import os
import re
from copy import deepcopy
from inspect import Parameter, Signature, signature
from typing import ClassVar, Optional, Union
from nibabel import load
from nipype.interfaces.utility import Function
from nipype.pipeline import engine as pe
from nipype.pipeline.engine.utils import (
    _create_dot_graph,
    format_dot,
    generate_expanded_graph,
    get_print_name,
    load_resultfile as _load_resultfile,
    _replacefunk,
    _run_dot
)
from nipype.utils.filemanip import fname_presuffix
from nipype.utils.functions import getsource
from numpy import prod
from traits.trait_base import Undefined
from traits.trait_handlers import TraitListObject
from CPAC.utils.monitoring.custom_logging import getLogger
from CPAC.utils.typing import DICT

# set global default mem_gb
DEFAULT_MEM_GB = 2.0
UNDEFINED_SIZE = (42, 42, 42, 1200)

logger = getLogger("nipype.workflow")


def _check_mem_x_path(mem_x_path):
    '''Function to check if a supplied multiplier path exists.

    Parameters
    ----------
    mem_x_path : str, iterable, Undefined or None

    Returns
    -------
    bool
    '''
    mem_x_path = _grab_first_path(mem_x_path)
    try:
        return mem_x_path is not Undefined and os.path.exists(
            mem_x_path)
    except (TypeError, ValueError):
        return False


def _doctest_skiplines(docstring, lines_to_skip):
    '''
    Function to add '  # doctest: +SKIP' to the end of docstring lines
    to skip in imported docstrings.

    Parameters
    ----------
    docstring : str

    lines_to_skip : set or list

    Returns
    -------
    docstring : str

    Examples
    --------
    >>> _doctest_skiplines('skip this line', {'skip this line'})
    'skip this line  # doctest: +SKIP'
    '''
    if (
        not isinstance(lines_to_skip, set) and
        not isinstance(lines_to_skip, list)
    ):
        raise TypeError(
            '_doctest_skiplines: `lines_to_skip` must be a set or list.')

    return '\n'.join([
        f'{line}  # doctest: +SKIP' if line in lines_to_skip else line
        for line in docstring.split('\n')
    ])


def _grab_first_path(mem_x_path):
    '''Function to grab the first path if multiple paths for given
    multiplier input

    Parameters
    ----------
    mem_x_path : str, iterable, Undefined or None

    Returns
    -------
    str, Undefined or None
    '''
    if isinstance(mem_x_path, (list, TraitListObject, tuple)):
        mem_x_path = mem_x_path[0] if len(mem_x_path) else Undefined
    return mem_x_path


class Node(pe.Node):
    # pylint: disable=empty-docstring,too-many-instance-attributes
    __doc__ = _doctest_skiplines(
        pe.Node.__doc__,
        {"    >>> realign.inputs.in_files = 'functional.nii'"}
    )

    def __init__(
            self,
            *args,
            mem_gb: Optional[float] = DEFAULT_MEM_GB,
            throttle: Optional[bool] = False,
            **kwargs
        ) -> None:
        # pylint: disable=import-outside-toplevel
        from CPAC.pipeline.random_state import random_seed
        super().__init__(*args, mem_gb=mem_gb, **kwargs)
        self.logger = getLogger("nipype.workflow")
        self.seed = random_seed()
        self.seed_applied = False
        self.input_data_shape = Undefined
        self._debug = False
        if throttle:
            self.throttle = True
        self.verbose_logger = None
        self._mem_x = {}
        if 'mem_x' in kwargs and isinstance(
            kwargs['mem_x'], (tuple, list)
        ):
            if len(kwargs['mem_x']) == 3:
                (
                    self._mem_x['multiplier'],
                    self._mem_x['file'],
                    self._mem_x['mode']
                ) = kwargs['mem_x']
            else:
                self._mem_x['mode'] = 'xyzt'
                if len(kwargs['mem_x']) == 2:
                    (
                        self._mem_x['multiplier'],
                        self._mem_x['file']
                    ) = kwargs['mem_x']
                else:
                    self._mem_x['multiplier'] = kwargs['mem_x']
                    self._mem_x['file'] = None
        else:
            delattr(self, '_mem_x')
        setattr(self, 'skip_timeout', False)

    orig_sig_params = list(signature(pe.Node).parameters.items())

    __init__.__signature__ = Signature(parameters=[
        p[1] if p[0] != 'mem_gb' else (
            'mem_gb',
            Parameter('mem_gb', Parameter.POSITIONAL_OR_KEYWORD,
                      default=DEFAULT_MEM_GB)
        )[1] for p in orig_sig_params[:-1]] + [
            Parameter('mem_x', Parameter.KEYWORD_ONLY, default=None),
            Parameter("throttle", Parameter.KEYWORD_ONLY, default=False),
            orig_sig_params[-1][1]
        ])

    __init__.__doc__ = re.sub(r'(?<!\s):', ' :', '\n'.join([
        pe.Node.__init__.__doc__.rstrip(),
        '''
        mem_gb : int or float
            Estimate (in GB) of constant memory to allocate for this
            node.

        mem_x : 2-tuple or 3-tuple
            (``multiplier``, ``input_file``)
            (int or float, str)

            (``multiplier``, ``input_file``, ``mode``)
            (int or float, str, str)

            **Note**
            This parameter (``mem_x``) is likely to change in a future
            release as we incorporate more factors into memory
            estimates.

            See also: `⚡️ Setting data- and operation-dependent memory-estimates <https://github.com/FCP-INDI/C-PAC/issues/1509>`_
                GitHub epic of issues related to improving Node
                memory estimates based on the data and operations
                involved.

            Multiplier for memory allocation such that ``multiplier``
            times ``mode`` of 4-D file at ``input_file`` plus
            ``self._mem_gb`` equals the total memory allocation for
            the node. ``input_file`` can be a Node input string or
            an actual path.

            ``mode`` can be any one of
            * 'xyzt' (spatial * temporal) (default if not specified)
            * 'xyz' (spatial)
            * 't' (temporal)

        throttle : bool, optional
            Assume this Node will use all available memory if no observation run is
            provided.''']))  # noqa: E501  # pylint: disable=line-too-long

    def _add_flags(self, flags):
        r'''
        Parameters
        ----------
        flags : list or tuple
            If a list, add ``flags`` to ``self.inputs.flags`` or
            ``self.inputs.args``

            If a tuple, remove ``flags[1]`` from and add ``flags[0]``
            to ``self.inputs.flags`` or ``self.inputs.args``
        '''
        def prep_flags(attr):
            to_remove = []
            if isinstance(flags, tuple):
                to_remove += flags[1]
                new_flags = flags[0]
            else:
                new_flags = flags
            old_flags = getattr(self.inputs, attr)
            if isinstance(old_flags, str):
                to_remove.sort(key=lambda x: -x.count(' '))
                for flag in to_remove:
                    if f' {flag} ' in old_flags:
                        old_flags = old_flags.replace(f' {flag}', '')
                old_flags = [old_flags]
            if isinstance(old_flags, list):
                new_flags = [flag for flag in old_flags if
                             flag not in to_remove] + new_flags
            if attr == 'args':
                new_flags = ' '.join(new_flags)
                while '  ' in new_flags:
                    new_flags = new_flags.replace('  ', ' ')
            return new_flags
        if hasattr(self.inputs, 'flags'):
            self.inputs.flags = prep_flags('flags')
        else:
            self.inputs.args = prep_flags('args')

    def _apply_mem_x(self, multiplicand=None):
        '''Method to calculate and memoize a Node's estimated memory
        footprint.

        Parameters
        ----------
        multiplicand : str or int or float or list thereof or
                       3-or-4-tuple or None
            Any of
            * path to file(s) with shape to multiply by multiplier
            * multiplicand
            * shape of image to consider with mode

        Returns
        -------
        number
            estimated memory usage (GB)
        '''
        def parse_multiplicand(multiplicand):
            '''
            Returns an numeric value for a multiplicand if
            multipland is a string or None.

            Parameters
            ----------
            muliplicand : any

            Returns
            -------
            int or float
            '''
            if self._debug:
                self.verbose_logger.debug('%s multiplicand: %s', self.name,
                                          multiplicand)
            if isinstance(multiplicand, list):
                return max([parse_multiplicand(part) for part in multiplicand])
            if isinstance(multiplicand, (int, float)):
                return multiplicand
            if (
                isinstance(multiplicand, tuple) and
                3 <= len(multiplicand) <= 4 and
                all(isinstance(i, (int, float)) for i in multiplicand)
            ):
                return get_data_size(
                    multiplicand,
                    getattr(self, '_mem_x', {}).get('mode'))
            if _check_mem_x_path(multiplicand):
                return get_data_size(
                    _grab_first_path(multiplicand),
                    getattr(self, '_mem_x', {}).get('mode'))
            return 1

        if hasattr(self, '_mem_x'):
            if self._debug:
                self.verbose_logger.debug('%s._mem_x: %s', self.name,
                                          self._mem_x)
            if multiplicand is None:
                multiplicand = self._mem_x_file()
            setattr(self, '_mem_gb', (
                self._mem_gb +
                self._mem_x.get('multiplier', 0) *
                parse_multiplicand(multiplicand)))
            try:
                if self._mem_gb > 1000:
                    self.logger.warning(
                        '%s is estimated to use %.3f GB (%s).',
                        self.name,
                        self._mem_gb,
                        getattr(self, '_mem_x')
                    )
            except FileNotFoundError:
                pass
            del self._mem_x
        if self._debug:
            self.verbose_logger.debug('%s._mem_gb: %s', self.name,
                                      self._mem_gb)
        return self._mem_gb

    def _apply_random_seed(self):
        '''Apply flags for the first matched interface'''
        # pylint: disable=import-outside-toplevel
        from CPAC.pipeline.random_state import random_seed_flags
        if isinstance(self.interface, Function):
            for rsf, flags in random_seed_flags()['functions'].items():
                if self.interface.inputs.function_str == getsource(rsf):
                    self.interface.inputs.function_str = flags(
                        self.interface.inputs.function_str)
                    self.seed_applied = True
                    return
        for rsf, flags in random_seed_flags()['interfaces'].items():
            if isinstance(self.interface, rsf):
                self._add_flags(flags)
                self.seed_applied = True
                return

    @property
    def mem_gb(self):
        """Get estimated memory (GB)"""
        if hasattr(self._interface, "estimated_memory_gb"):
            self._mem_gb = self._interface.estimated_memory_gb
            self.logger.warning(
                'Setting "estimated_memory_gb" on Interfaces has been '
                "deprecated as of nipype 1.0, please use Node.mem_gb."
            )
        if hasattr(self, '_mem_x'):
            if self._mem_x['file'] is None:
                return self._apply_mem_x()
            try:
                mem_x_path = getattr(self.inputs, self._mem_x['file'])
            except AttributeError as attribute_error:
                raise AttributeError(
                    f'{attribute_error.args[0]} in Node \'{self.name}\''
                ) from attribute_error
            if _check_mem_x_path(mem_x_path):
                # constant + mem_x[0] * t
                return self._apply_mem_x()
            raise FileNotFoundError(2, 'The memory estimate for Node '
                                    f"'{self.name}' depends on the input "
                                    f"'{self._mem_x['file']}' but "
                                    'no such file or directory', mem_x_path)
        return self._mem_gb

    @property
    def mem_x(self):
        """Get dict of 'multiplier' (memory multiplier), 'file' (input file)
        and multiplier mode (spatial * temporal, spatial only or
        temporal only). Returns ``None`` if already consumed or not set."""
        return getattr(self, '_mem_x', None)

    def _mem_x_file(self):
        return getattr(self.inputs, getattr(self, '_mem_x', {}).get('file'))

    def override_mem_gb(self, new_mem_gb):
        """Override the Node's memory estimate with a new value.

        Parameters
        ----------
        new_mem_gb : int or float
            new memory estimate in GB
        """
        if hasattr(self, '_mem_x'):
            delattr(self, '_mem_x')
        setattr(self, '_mem_gb', new_mem_gb)

    def run(self, updatehash=False):
        self.__doc__ = getattr(super(), '__doc__', '')
        if self.seed is not None:
            self._apply_random_seed()
            if self.seed_applied:
                random_state_logger = getLogger('random')
                random_state_logger.info('%s\t%s', '# (Atropos constant)' if
                                         'atropos' in self.name else
                                         str(self.seed), self.name)
        return super().run(updatehash)


class MapNode(Node, pe.MapNode):
    # pylint: disable=empty-docstring
    __doc__ = _doctest_skiplines(
        pe.MapNode.__doc__,
        {"    ...                           'functional3.nii']"}
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self.name.endswith('_'):
            self.name = f'{self.name}_'

    _parameters: ClassVar[DICT[str, Parameter]] = {}
    _custom_params: ClassVar[DICT[str, Union[bool, float]]] = {
        "mem_gb": DEFAULT_MEM_GB,
        "throttle": False,
    }
    for param, default in _custom_params.items():
        for p in signature(pe.Node).parameters.items():
            if p[0] in _custom_params:
                _parameters[p[0]] = Parameter(
                    param, Parameter.POSITIONAL_OR_KEYWORD, default=default
                )
            else:
                _parameters[p[0]] = p[1]
    __init__.__signature__ = Signature(parameters=list(_parameters.values()))
    del _custom_params, _parameters


class Workflow(pe.Workflow):
    """Controls the setup and execution of a pipeline of processes."""

    def __init__(self, name, base_dir=None, debug=False):
        """Create a workflow object.
        Parameters
        ----------
        name : alphanumeric string
            unique identifier for the workflow
        base_dir : string, optional
            path to workflow storage
        debug : boolean, optional
            enable verbose debug-level logging
        """
        import networkx as nx

        super().__init__(name, base_dir)
        self._debug = debug
        self.verbose_logger = getLogger('engine') if debug else None
        self._graph = nx.DiGraph()

        self._nodes_cache = set()
        self._nested_workflows_cache = set()

    def _configure_exec_nodes(self, graph):
        """Ensure that each node knows where to get inputs from"""
        for node in graph.nodes():
            node._debug = self._debug  # pylint: disable=protected-access
            node.verbose_logger = self.verbose_logger
            node.input_source = {}
            for edge in graph.in_edges(node):
                data = graph.get_edge_data(*edge)
                for sourceinfo, field in data["connect"]:
                    node.input_source[field] = (
                        os.path.join(edge[0].output_dir(),
                                     "result_%s.pklz" % edge[0].name),
                        sourceinfo,
                    )
                    if node and hasattr(node, '_mem_x'):
                        if isinstance(
                            node._mem_x,  # pylint: disable=protected-access
                            dict
                        ) and node._mem_x[  # pylint: disable=protected-access
                                          'file'] == field:
                            input_resultfile = node.input_source.get(field)
                            if input_resultfile:
                                # pylint: disable=protected-access
                                if isinstance(input_resultfile, tuple):
                                    input_resultfile = input_resultfile[0]
                                try:
                                    # memoize node._mem_gb if path
                                    # already exists
                                    node._apply_mem_x(_load_resultfile(
                                        input_resultfile
                                    ).inputs[field])
                                except (FileNotFoundError, KeyError,
                                        TypeError):
                                    self._handle_just_in_time_exception(node)

    def _get_dot(
        self, prefix=None, hierarchy=None, colored=False, simple_form=True,
        level=0
    ):
        """Create a dot file with connection info"""
        # pylint: disable=invalid-name,protected-access
        import networkx as nx

        if prefix is None:
            prefix = "  "
        if hierarchy is None:
            hierarchy = []
        colorset = [
            "#FFFFC8",  # Y
            "#0000FF",
            "#B4B4FF",
            "#E6E6FF",  # B
            "#FF0000",
            "#FFB4B4",
            "#FFE6E6",  # R
            "#00A300",
            "#B4FFB4",
            "#E6FFE6",  # G
            "#0000FF",
            "#B4B4FF",
        ]  # loop B
        if level > len(colorset) - 2:
            level = 3  # Loop back to blue
        quoted_prefix = f'"{prefix}"' if len(prefix.strip()) else prefix
        dotlist = [f'{quoted_prefix}label="{self.name}";']
        for node in nx.topological_sort(self._graph):
            fullname = ".".join(hierarchy + [node.fullname])
            nodename = fullname.replace(".", "_")
            if not isinstance(node, Workflow):
                node_class_name = get_print_name(node, simple_form=simple_form)
                if not simple_form:
                    node_class_name = ".".join(node_class_name.split(".")[1:])
                if hasattr(node, "iterables") and node.iterables:
                    dotlist.append(f'"{nodename}"[label="{node_class_name}", '
                                   "shape=box3d, style=filled, color=black, "
                                   "colorscheme=greys7 fillcolor=2];")
                else:
                    if colored:
                        dotlist.append(f'"{nodename}"[label="'
                                       f'{node_class_name}", style=filled,'
                                       f' fillcolor="{colorset[level]}"];')
                    else:
                        dotlist.append(f'"{nodename}"[label="'
                                       f'{node_class_name}"];')

        for node in nx.topological_sort(self._graph):
            if isinstance(node, Workflow):
                fullname = ".".join(hierarchy + [node.fullname])
                nodename = fullname.replace(".", "_")
                dotlist.append(f"subgraph \"cluster_{nodename}\" {{")
                if colored:
                    dotlist.append(f'{prefix}{prefix}edge [color="'
                                   f'{colorset[level + 1]}"];')
                    dotlist.append(f"{prefix}{prefix}style=filled;")
                    dotlist.append(f'{prefix}{prefix}fillcolor='
                                   f'"{colorset[level + 2]}";')
                dotlist.append(
                    node._get_dot(
                        prefix=prefix + prefix,
                        hierarchy=hierarchy + [self.name],
                        colored=colored,
                        simple_form=simple_form,
                        level=level + 3,
                    )
                )
                dotlist.append("}")
            else:
                for subnode in self._graph.successors(node):
                    if node._hierarchy != subnode._hierarchy:
                        continue
                    if not isinstance(subnode, Workflow):
                        nodefullname = ".".join(hierarchy + [node.fullname])
                        subnodefullname = ".".join(
                            hierarchy + [subnode.fullname])
                        nodename = nodefullname.replace(".", "_")
                        subnodename = subnodefullname.replace(".", "_")
                        for _ in self._graph.get_edge_data(
                            node, subnode
                        )["connect"]:
                            dotlist.append(f'"{nodename}" -> "{subnodename}";')
                        logger.debug("connection: %s", dotlist[-1])
        # add between workflow connections
        for u, v, d in self._graph.edges(data=True):
            uname = ".".join(hierarchy + [u.fullname])
            vname = ".".join(hierarchy + [v.fullname])
            for src, dest in d["connect"]:
                uname1 = uname
                vname1 = vname
                if isinstance(src, tuple):
                    srcname = src[0]
                else:
                    srcname = src
                if "." in srcname:
                    uname1 += "." + ".".join(srcname.split(".")[:-1])
                if "." in dest and "@" not in dest:
                    if not isinstance(v, Workflow):
                        if "datasink" not in str(
                            v._interface.__class__
                        ).lower():
                            vname1 += "." + ".".join(dest.split(".")[:-1])
                    else:
                        vname1 += "." + ".".join(dest.split(".")[:-1])
                if uname1.split(".")[:-1] != vname1.split(".")[:-1]:
                    dotlist.append(f'"{uname1.replace(".", "_")}" -> '
                                   f'"{vname1.replace(".", "_")}";')
                    logger.debug("cross connection: %s", dotlist[-1])
        return ("\n" + prefix).join(dotlist)

    def _handle_just_in_time_exception(self, node):
        # pylint: disable=protected-access
        if hasattr(self, '_local_func_scans'):
            node._apply_mem_x(
                self._local_func_scans)  # pylint: disable=no-member
        else:
            # TODO: handle S3 files
            node._apply_mem_x(UNDEFINED_SIZE)  # noqa: W0212

    def write_graph(
        self,
        dotfilename="graph.dot",
        graph2use="hierarchical",
        format="png",
        simple_form=True,
    ):
        graphtypes = ["orig", "flat", "hierarchical", "exec", "colored"]
        if graph2use not in graphtypes:
            raise ValueError(
                "Unknown graph2use keyword. Must be one of: " + str(graphtypes)
            )
        base_dir, dotfilename = os.path.split(dotfilename)
        if base_dir == "":
            if self.base_dir:
                base_dir = self.base_dir
                if self.name:
                    base_dir = os.path.join(base_dir, self.name)
            else:
                base_dir = os.getcwd()
        os.makedirs(base_dir, exist_ok=True)
        if graph2use in ["hierarchical", "colored"]:
            if self.name[:1].isdigit():  # these graphs break if int
                raise ValueError(f"{graph2use} graph failed, workflow name "
                                 "cannot begin with a number")
            dotfilename = os.path.join(base_dir, dotfilename)
            self.write_hierarchical_dotfile(
                dotfilename=dotfilename,
                colored=graph2use == "colored",
                simple_form=simple_form,
            )
            outfname = format_dot(dotfilename, format=format)
        else:
            graph = self._graph
            if graph2use in ["flat", "exec"]:
                graph = self._create_flat_graph()
            if graph2use == "exec":
                graph = generate_expanded_graph(deepcopy(graph))
            outfname = export_graph(
                graph,
                base_dir,
                dotfilename=dotfilename,
                format=format,
                simple_form=simple_form,
            )

        logger.info("Generated workflow graph: %s "
                    "(graph2use=%s, simple_form=%s).",
                    outfname, graph2use, simple_form)
        return outfname

    write_graph.__doc__ = pe.Workflow.write_graph.__doc__

    def write_hierarchical_dotfile(
        self, dotfilename=None, colored=False, simple_form=True
    ):
        # pylint: disable=invalid-name
        dotlist = [f"digraph \"{self.name}\"{{"]
        dotlist.append(self._get_dot(prefix="  ", colored=colored,
                                     simple_form=simple_form))
        dotlist.append("}")
        dotstr = "\n".join(dotlist)
        if dotfilename:
            with open(dotfilename, "wt", encoding="utf-8") as fp:
                fp.writelines(dotstr)
                fp.close()
        else:
            logger.info(dotstr)


def get_data_size(filepath, mode='xyzt'):
    """Function to return the size of a functional image (x * y * z * t)

    Parameters
    ----------
    filepath : str or path
        path to image file
        OR
        4-tuple
        stand-in dimensions (x, y, z, t)

    mode : str
        One of:
        * 'xyzt' (all dimensions multiplied) (DEFAULT)
        * 'xyz' (spatial dimensions multiplied)
        * 't' (number of TRs)

    Returns
    -------
    int or float
    """
    if isinstance(filepath, str):
        data_shape = load(filepath).shape
    elif isinstance(filepath, tuple) and len(filepath) == 4:
        data_shape = filepath
    if mode == 't':
        # if the data has muptiple TRs, return that number
        if len(data_shape) > 3:
            return data_shape[3]
        # otherwise return 1
        return 1
    if mode == 'xyz':
        return prod(data_shape[0:3]).item()
    return prod(data_shape).item()


def export_graph(
    graph_in,
    base_dir=None,
    show=False,
    use_execgraph=False,
    show_connectinfo=False,
    dotfilename="graph.dot",
    format="png",
    simple_form=True,
):
    """Displays the graph layout of the pipeline
    This function requires that pygraphviz and matplotlib are available on
    the system.
    Parameters
    ----------
    show : boolean
    Indicate whether to generate pygraphviz output fromn
    networkx. default [False]
    use_execgraph : boolean
    Indicates whether to use the specification graph or the
    execution graph. default [False]
    show_connectioninfo : boolean
    Indicates whether to show the edge data on the graph. This
    makes the graph rather cluttered. default [False]
    """
    import networkx as nx

    graph = deepcopy(graph_in)
    if use_execgraph:
        graph = generate_expanded_graph(graph)
        logger.debug("using execgraph")
    else:
        logger.debug("using input graph")
    if base_dir is None:
        base_dir = os.getcwd()

    os.makedirs(base_dir, exist_ok=True)
    out_dot = fname_presuffix(dotfilename, suffix="_detailed.dot",
                              use_ext=False, newpath=base_dir)
    _write_detailed_dot(graph, out_dot)

    # Convert .dot if format != 'dot'
    outfname, res = _run_dot(out_dot, format_ext=format)
    if res is not None and res.runtime.returncode:
        logger.warning("dot2png: %s", res.runtime.stderr)

    pklgraph = _create_dot_graph(graph, show_connectinfo, simple_form)
    simple_dot = fname_presuffix(dotfilename, suffix=".dot", use_ext=False,
                                 newpath=base_dir)
    nx.drawing.nx_pydot.write_dot(pklgraph, simple_dot)

    # Convert .dot if format != 'dot'
    simplefname, res = _run_dot(simple_dot, format_ext=format)
    if res is not None and res.runtime.returncode:
        logger.warning("dot2png: %s", res.runtime.stderr)

    if show:
        pos = nx.graphviz_layout(pklgraph, prog="dot")
        nx.draw(pklgraph, pos)
        if show_connectinfo:
            nx.draw_networkx_edge_labels(pklgraph, pos)

    return simplefname if simple_form else outfname


def _write_detailed_dot(graph, dotfilename):
    r"""
    Create a dot file with connection info ::
        digraph structs {
        node [shape=record];
        struct1 [label="<f0> left|<f1> middle|<f2> right"];
        struct2 [label="<f0> one|<f1> two"];
        struct3 [label="hello\nworld |{ b |{c|<here> d|e}| f}| g | h"];
        struct1:f1 -> struct2:f0;
        struct1:f0 -> struct2:f1;
        struct1:f2 -> struct3:here;
        }
    """
    # pylint: disable=invalid-name
    import networkx as nx

    text = ["digraph structs {", "node [shape=record];"]
    # write nodes
    edges = []
    for n in nx.topological_sort(graph):
        nodename = n.itername
        inports = []
        for u, v, d in graph.in_edges(nbunch=n, data=True):
            for cd in d["connect"]:
                if isinstance(cd[0], (str, bytes)):
                    outport = cd[0]
                else:
                    outport = cd[0][0]
                inport = cd[1]
                ipstrip = f"in{_replacefunk(inport)}"
                opstrip = f"out{_replacefunk(outport)}"
                edges.append(f'"{u.itername.replace(".", "")}":'
                             f'"{opstrip}":e -> '
                             f'"{v.itername.replace(".", "")}":'
                             f'"{ipstrip}":w;')
                if inport not in inports:
                    inports.append(inport)
        inputstr = (["{IN"]
                    + [f"|<in{_replacefunk(ip)}> {ip}" for
                       ip in sorted(inports)] + ["}"])
        outports = []
        for u, v, d in graph.out_edges(nbunch=n, data=True):
            for cd in d["connect"]:
                if isinstance(cd[0], (str, bytes)):
                    outport = cd[0]
                else:
                    outport = cd[0][0]
                if outport not in outports:
                    outports.append(outport)
        outputstr = (
            ["{OUT"]
            + [f"|<out{_replacefunk(oport)}> {oport}" for
               oport in sorted(outports)] + ["}"])
        srcpackage = ""
        if hasattr(n, "_interface"):
            pkglist = n.interface.__class__.__module__.split(".")
            if len(pkglist) > 2:
                srcpackage = pkglist[2]
        srchierarchy = ".".join(nodename.split(".")[1:-1])
        nodenamestr = (f"{{ {nodename.split('.')[-1]} | {srcpackage} | "
                       f"{srchierarchy} }}")
        text += [f'"{nodename.replace(".", "")}" [label='
                 f'"{"".join(inputstr)}|{nodenamestr}|{"".join(outputstr)}"];']
    # write edges
    for edge in sorted(edges):
        text.append(edge)
    text.append("}")
    with open(dotfilename, "wt", encoding="utf-8") as filep:
        filep.write("\n".join(text))
    return text
