'''Module to import Nipype Pipeline engine and override some Classes.
See https://fcp-indi.github.io/docs/developer/nodes
for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.pipeline.engine.html
for Nipype's documentation.
'''  # noqa E501
import os
import re
import warnings
from logging import getLogger
from inspect import Parameter, Signature, signature
from nibabel import load
from nipype import logging
from nipype.interfaces.utility import Function
from nipype.pipeline import engine as pe
from nipype.pipeline.engine.utils import load_resultfile as _load_resultfile
from nipype.utils.functions import getsource
from numpy import prod
from traits.trait_base import Undefined
from traits.trait_handlers import TraitListObject

# set global default mem_gb
DEFAULT_MEM_GB = 2.0
UNDEFINED_SIZE = (42, 42, 42, 1200)

random_state_logger = getLogger('random')


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


class Node(pe.Node):
    __doc__ = _doctest_skiplines(
        pe.Node.__doc__,
        {"    >>> realign.inputs.in_files = 'functional.nii'"}
    )

    def __init__(self, *args, mem_gb=DEFAULT_MEM_GB, **kwargs):
        from CPAC.pipeline.random_state import random_seed
        super().__init__(*args, mem_gb=mem_gb, **kwargs)
        self.logger = logging.getLogger("nipype.workflow")
        self.seed = random_seed()
        self.seed_applied = False

        if 'mem_x' in kwargs and isinstance(
            kwargs['mem_x'], (tuple, list)
        ):
            mem_x = {}
            if len(kwargs['mem_x']) == 3:
                (
                    mem_x['multiplier'],
                    mem_x['file'],
                    mem_x['mode']
                ) = kwargs['mem_x']
            else:
                mem_x['mode'] = 'xyzt'
                if len(kwargs['mem_x']) == 2:
                    (
                        mem_x['multiplier'],
                        mem_x['file']
                    ) = kwargs['mem_x']
                else:
                    mem_x['multiplier'] = kwargs['mem_x']
                    mem_x['file'] = None
            setattr(self, '_mem_x', mem_x)
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
            * 't' (temporal)''']))  # noqa E501

    def _add_flags(self, flags):
        '''
        Parameters
        ----------
        flags : list
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
            except AttributeError as e:
                raise AttributeError(
                    f'{e.args[0]} in Node \'{self.name}\'') from e
            if self._check_mem_x_path(mem_x_path):
                # constant + mem_x[0] * t
                return self._apply_mem_x()
            raise FileNotFoundError(2, 'The memory estimate for Node '
                                    f"'{self.name}' depends on the input "
                                    f"'{self._mem_x['file']}' but "
                                    'no such file or directory', mem_x_path)
        return self._mem_gb

    def _check_mem_x_path(self, mem_x_path):
        '''Method to check if a supplied multiplier path exists.

        Parameters
        ----------
        mem_x_path : str, iterable, Undefined or None

        Returns
        -------
        bool
        '''
        mem_x_path = self._grab_first_path(mem_x_path)
        try:
            return mem_x_path is not Undefined and os.path.exists(
                mem_x_path)
        except (TypeError, ValueError):
            return False

    def _grab_first_path(self, mem_x_path):
        '''Method to grab the first path if multiple paths for given
        multiplier input

        Parameters
        ----------
        mem_x_path : str, iterable, Undefined or None

        Returns
        -------
        str, Undefined or None
        '''
        if (
            isinstance(mem_x_path, list) or
            isinstance(mem_x_path, TraitListObject) or
            isinstance(mem_x_path, tuple)
        ):
            mem_x_path = mem_x_path[0] if len(mem_x_path) else Undefined
        return mem_x_path

    def _mem_x_file(self):
        return getattr(self.inputs, getattr(self, '_mem_x', {}).get('file'))

    def _apply_mem_x(self, multiplicand=None):
        '''Method to calculate and memoize a Node's estimated memory
        footprint.

        Parameters
        ----------
        multiplier : str or int or float or list thereof or 4-tuple or None
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
            if isinstance(multiplicand, list):
                return sum([parse_multiplicand(part) for part in multiplicand])
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
            if self._check_mem_x_path(multiplicand):
                return get_data_size(
                    self._grab_first_path(multiplicand),
                    getattr(self, '_mem_x', {}).get('mode'))
            return 1

        if hasattr(self, '_mem_x'):
            if multiplicand is None:
                multiplicand = self._mem_x_file()
            self._mem_gb = (
                self._mem_gb +
                self._mem_x['multiplier'] *  # pylint: disable=no-member
                parse_multiplicand(multiplicand)
            )
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
        return self._mem_gb

    @property
    def mem_x(self):
        """Get dict of 'multiplier' (memory multiplier), 'file' (input file)
        and multiplier mode (spatial * temporal, spatial only or
        temporal only). Returns ``None`` if already consumed or not set."""
        return getattr(self, '_mem_x', None)

    def run(self, updatehash=False):
        if self.seed is not None:
            self._apply_random_seed()
            if self.seed_applied:
                random_state_logger.info('%s',
                                         '%s  # (Atropos constant)' %
                                         self.name if 'atropos' in
                                         self.name else self.name)
        return super().run(updatehash)


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
    def _configure_exec_nodes(self, graph):
        """Ensure that each node knows where to get inputs from"""
        for node in graph.nodes():
            node.input_source = {}
            for edge in graph.in_edges(node):
                data = graph.get_edge_data(*edge)
                for sourceinfo, field in data["connect"]:
                    node.input_source[field] = (
                        os.path.join(edge[0].output_dir(),
                                     "result_%s.pklz" % edge[0].name),
                        sourceinfo,
                    )
                    if node and hasattr(node, 'mem_x'):
                        if isinstance(
                            node.mem_x,
                            dict
                        ) and node.mem_x['file'] == field:
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
                                except FileNotFoundError:
                                    self._handle_just_in_time_exception(node)
                                except KeyError:
                                    warnings.warn(str(KeyError(
                                        f'Node {node.name} specifies memory '
                                        'allocation for input '
                                        f'{node.mem_x["file"]}, but no such '
                                        'input is specified for that Node.')))
                                    self._handle_just_in_time_exception(node)

    def _handle_just_in_time_exception(self, node):
        # pylint: disable=protected-access
        if hasattr(self, '_local_func_scans'):
            node._apply_mem_x(
                self._local_func_scans)  # pylint: disable=no-member
        else:
            # TODO: handle S3 files
            node._apply_mem_x(UNDEFINED_SIZE)  # noqa W0212


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
