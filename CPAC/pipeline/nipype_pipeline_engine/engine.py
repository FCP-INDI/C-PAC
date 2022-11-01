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

# Modifications Copyright (C) 2022 C-PAC Developers

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
from inspect import Parameter, Signature, signature
from logging import getLogger
from typing import Iterable, Tuple, Union
from nibabel import load
from nipype import logging
from nipype.interfaces.base.support import Bunch, InterfaceResult
from nipype.interfaces.utility import Function, Merge, Select
from nipype.pipeline import engine as pe
from nipype.pipeline.engine.utils import load_resultfile as _load_resultfile
from nipype.utils.functions import getsource
from numpy import prod
from traits.trait_base import Undefined
from traits.trait_handlers import TraitListObject
from CPAC.pipeline.random_state.seed import increment_seed, random_seed, \
                                            random_seed_flags
from CPAC.registration.guardrails import BestOf, registration_guardrail, \
                                         skip_if_first_try_succeeds

# set global default mem_gb
DEFAULT_MEM_GB = 2.0
UNDEFINED_SIZE = (42, 42, 42, 1200)

random_state_logger = getLogger('random')
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

    def __init__(self, *args, mem_gb=DEFAULT_MEM_GB, **kwargs):
        super().__init__(*args, mem_gb=mem_gb, **kwargs)
        self.logger = logging.getLogger("nipype.workflow")
        self.seed = random_seed()
        self.seed_applied = False
        self.input_data_shape = Undefined
        self._debug = False
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
            Parameter('mem_x', Parameter.KEYWORD_ONLY),
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
            * 't' (temporal)''']))  # noqa: E501  # pylint: disable=line-too-long

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
        if hasattr(self.interface, 'inputs'
                   ) and hasattr(self.interface.inputs, 'previous_failure'):
            if self.interface.inputs.previous_failure is False:
                return InterfaceResult(self.interface, Bunch(),
                                       self.inputs, self.outputs, None)
        if self.seed is not None:
            self._apply_random_seed()
            if self.seed_applied:
                random_state_logger.info('%s\t%s', '# (Atropos constant)' if
                                         'atropos' in self.name else
                                         str(self.seed), self.name)
        return super().run(updatehash)


class GuardrailedNode:
    '''A Node with QC guardrails.'''
    def __init__(self, wf, node, reference, registered):
        '''A Node with guardrails

        Parameters
        ----------
        wf : Workflow
            The parent workflow in which this node is guardrailed

        node : Node
            Node to guardrail

        reference : str
            key for reference image

        registered : str
            key for registered image
        '''
        self.guardrails = [registration_guardrail_node(
            f'{node.name}_guardrail')]
        self.node = node
        self.reference = reference
        self.registered = registered
        self.retries = []
        self.wf = wf
        self.wf.connect(self.node, registered,
                        self.guardrails[0], 'registered')
        if self.wf.num_tries > 1:
            if self.wf.retry_on_first_failure:
                self.guardrails.append(registration_guardrail_node(
                    f'{node.name}_guardrail'))
                self.retries.append(retry_clone(self.node))
            else:
                for i in range(self.wf.num_tries - 1):
                    self.retries.append(retry_clone(self.node, i + 2))
                    self.retries[0].interface = skip_if_first_try_succeeds(
                        self.retries[0].interface)
                    self.wf.connect(self.node, 'failed_qc',
                                    self.retries[0], 'previous_failure')
            for i, retry in enumerate(self.retries):
                self.wf.connect(retry, registered,
                                self.guardrails[i + 1], 'registered')

    def guardrail_selection(self, node, output_key):
        """Convenience method to :py:method:`Workflow.guardrail_selection`"""
        return self.wf.guardrail_selection(node, output_key)


class MapNode(Node, pe.MapNode):
    # pylint: disable=empty-docstring
    __doc__ = _doctest_skiplines(
        pe.MapNode.__doc__,
        {"    ...                           'functional3.nii']"}
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    __init__.__signature__ = Signature(parameters=[
        p[1] if p[0] != 'mem_gb' else (
            'mem_gb',
            Parameter('mem_gb', Parameter.POSITIONAL_OR_KEYWORD,
                      default=DEFAULT_MEM_GB)
        )[1] for p in signature(pe.Node).parameters.items()])


class Workflow(pe.Workflow):
    """Controls the setup and execution of a pipeline of processes."""

    def __init__(self, name, base_dir=None, debug=False, guardrail_config=None
                 ):
        """Create a workflow object.
        Parameters
        ----------
        name : alphanumeric string
            unique identifier for the workflow
        base_dir : string, optional
            path to workflow storage
        debug : boolean, optional
            enable verbose debug-level logging
        guardrail_config : dict, optional
        """
        import networkx as nx

        super().__init__(name, base_dir)
        self._debug = debug
        self.verbose_logger = getLogger('engine') if debug else None
        self._graph = nx.DiGraph()
        self._guardrail_config = {
        } if guardrail_config is None else guardrail_config
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

    def connect(self, *args, **kwargs):
        """Connects all the retry nodes and guardrails of guardrailed nodes,
        then connects the other nodes as usual

        .. seealso:: pe.Node.connect
        """
        if len(args) == 1:
            connection_list = args.pop(0)
        elif len(args) == 4:
            connection_list = [(args.pop(0), args.pop(2), [
                (args.pop(1), args.pop(3))])]
        new_connection_list = []
        for srcnode, destnode, connects in connection_list:
            if isinstance(srcnode, GuardrailedNode):
                for _from, _to in connects:
                    selected = srcnode.guardrail_selection(srcnode, _from)
                    self.connect(selected, 'out', destnode, _to)
            if isinstance(destnode, GuardrailedNode):
                for _from, _to in connects:
                    guardnodes = [destnode.node, *destnode.retries]
                    self.connect_many(guardnodes, [(srcnode, _from, _to)])
                    if _from == destnode.reference:
                        self.connect_many(guardnodes, [
                            (srcnode, _from, 'reference')])
            else:
                new_connection_list.extend([srcnode, destnode, connects])
        super().connect(new_connection_list, *args, **kwargs)

    def connect_many(self, nodes: Iterable['Node'],
                     connections: Iterable[Tuple['Node', Union[str, tuple],
                                                 str]]) -> None:
        """Method to generalize making the same connections to try and
        retry nodes.

        For each 3-tuple (``conn``) in ``connections``, will do
        ``wf.connect(conn[0], conn[1], node, conn[2])`` for each ``node``
        in ``nodes``

        Parameters
        ----------
        nodes : iterable of Nodes

        connections : iterable of 3-tuples of (Node, str or tuple, str)
        """
        wrong_conn_type_msg = (r'connect_many `connections` argument '
                               'must be an iterable of (Node, str or '
                               'tuple, str) tuples.')
        if not isinstance(connections, (list, tuple)):
            raise TypeError(f'{wrong_conn_type_msg}: Given {connections}')
        for node in nodes:
            if not isinstance(node, Node):
                raise TypeError('connect_many requires an iterable '
                                r'of nodes for the `nodes` parameter: '
                                f'Given {node}')
            for conn in connections:
                if not all((isinstance(conn, (list, tuple)), len(conn) == 3,
                            isinstance(conn[0], Node),
                            isinstance(conn[1], (tuple, str)),
                            isinstance(conn[2], str))):
                    raise TypeError(f'{wrong_conn_type_msg}: Given {conn}')
                self.connect(*conn[:2], node, conn[2])

    @property
    def guardrail(self):
        """Are guardrails on?

        Returns
        -------
        boolean
        """
        return any(self._guardrail_config['thresholds'].values())

    def guardrailed_node(self, node, reference, registered):
        """Method to return a GuardrailedNode in the given Workflow.

        .. seealso:: GuardrailedNode
        """
        return GuardrailedNode(self, node, reference, registered)

    def guardrail_selection(self, node: 'GuardrailedNode', output_key: str
                            ) -> Node:
        """Generate requisite Nodes for choosing a path through the graph
        with retries.

        Takes two nodes to choose an output from. These nodes are assumed
        to be guardrail nodes if `output_key` and `guardrail_node` are not
        specified.

        A ``nipype.interfaces.utility.Merge`` is generated, connecting
        ``output_key`` from ``node1`` and ``node2`` in that order.

        A ``nipype.interfaces.utility.Select`` node is generated taking the
        output from the generated ``Merge`` and using the ``failed_qc``
        output of ``guardrail_node`` (``node1`` if ``guardrail_node`` is
        unspecified).

        All relevant connections are made in the given Workflow.

        The ``Select`` node is returned; its output is keyed ``out`` and
        contains the value of the given ``output_key`` (``registered`` if
        unspecified).

        Parameters
        ----------
        node : GuardrailedNode

        output_key : key to select from Node

        Returns
        -------
        select : Node
        """
        name = node.node.name
        if output_key != 'registered':
            name = f'{name}_{output_key}'
        choices = Node(Merge(self.num_tries), run_without_submitting=True,
                       name=f'{name}_choices')
        select = Node(Select(), run_without_submitting=True,
                      name=f'choose_{name}')
        self.connect([(node.node, choices, [(output_key, 'in1')]),
                      (choices, select, [('out', 'inlist')])])
        if self._guardrail_config['best_of'] > 1:
            best_of = Node(BestOf(len(self.num_tries)))
            self.connect([(node.guardrail, best_of, [('error', 'error1')]),
                          (best_of, 'index', [(select, 'index')])])
            for i, retry in enumerate(node.retries):
                self.connect([(retry, choices, [(output_key, f'in{i+2}')]),
                              (retry.guardrail, best_of, [('error',
                                                           f'error{i+2}')])])
        elif self.retry_on_first_failure:
            self.connect([(node.retries[0], choices, [(output_key, 'in2')]),
                          (node.guardrail, select, [('failed_qc', 'index')])])
        return select

    def _handle_just_in_time_exception(self, node):
        # pylint: disable=protected-access
        if hasattr(self, '_local_func_scans'):
            node._apply_mem_x(
                self._local_func_scans)  # pylint: disable=no-member
        else:
            # TODO: handle S3 files
            node._apply_mem_x(UNDEFINED_SIZE)  # noqa: W0212

    @property
    def num_tries(self):
        """How many maximum tries?

        Returns
        -------
        int
        """
        if self.guardrail is False:
            return 1
        return 2 if (self.retry_on_first_failure
                     ) else self._guardrail_config['best_of']

    @property
    def retry_on_first_failure(self):
        """Retry iff first attempt fails?

        Returns
        -------
        bool
        """
        return (self._guardrail_config['best_of'] == 1 and
                self._guardrail_config['retry_on_first_failure'] is True)


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


def registration_guardrail_node(name=None, retry_num=0):
    """Convenience method to get a new registration_guardrail Node

    Parameters
    ----------
    name : str, optional

    retry_num : int, optional
        how many previous tries?

    Returns
    -------
    Node
    """
    if name is None:
        name = 'registration_guardrail'
    node = Node(Function(input_names=['registered', 'reference', 'retry_num'],
                         output_names=['registered', 'failed_qc', 'error'],
                         imports=['import logging',
                                  'from typing import Tuple',
                                  'from CPAC.qc import qc_masks, '
                                  'registration_guardrail_thresholds',
                                  'from CPAC.registration.guardrails '
                                  'import BadRegistrationError'],
                         function=registration_guardrail), name=name)
    if retry_num:
        node.inputs.retry_num = retry_num
    return node


def retry_clone(node: 'Node', index: int = 1) -> 'Node':
    """Function to clone a node, name the clone, and increment its
    random seed

    Parameters
    ----------
    node : Node

    index : int
        if multiple tries regardless of initial success, nth try
        (starting with 2)

    Returns
    -------
    Node
    """
    if index > 1:
        return increment_seed(node.clone(f'{node.name}_try{index}'))
    return increment_seed(node.clone(f'retry_{node.name}'))
