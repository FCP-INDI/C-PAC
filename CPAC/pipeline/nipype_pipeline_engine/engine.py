'''Module to import Nipype Pipeline engine and override some Classes.
See https://fcp-indi.github.com/docs/developer/nodes
for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.pipeline.engine.html
for Nipype's documentation.
'''  # noqa E501
import os
import re
from inspect import Parameter, Signature, signature
from nibabel import load
from nipype.pipeline import engine as pe
from nipype.pipeline.engine.utils import load_resultfile as _load_resultfile
from numpy import prod
from traits.trait_base import Undefined
from traits.trait_handlers import TraitListObject

# set global default mem_gb
DEFAULT_MEM_GB = 2.0


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
        super().__init__(*args, mem_gb=mem_gb, **kwargs)
        if 'mem_x' in kwargs:
            setattr(self, '_mem_x', kwargs['mem_x'])

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

        mem_x : tuple or number

            tuple (len 2)
                (multiplier, path)
                (int or float, str)

                Multiplier for memory allocation such that `multiplier`
                times x * y * z * t of 4-D file at `path` plus
                `self._mem_gb` equals the total memory allocation for
                the node.

            tuple (len 1) or number
                (multiplier, ) or multiplier
                (number, ) or number
                Multiplier for memory allocation such that `multiplier`
                times `wf._largest_func` (x * y * z * t of the largest
                functional file in `wf` where `wf` is the Workflow
                containing this Node, if this Node is in a Workflow and
                `wf` has the attribute `_largest_func`)
            ''']))

    @property
    def mem_gb(self):
        """Get estimated memory (GB)"""
        if hasattr(self._interface, "estimated_memory_gb"):
            from nipype import logging
            logger = logging.getLogger("nipype.workflow")
            self._mem_gb = self._interface.estimated_memory_gb
            logger.warning(
                'Setting "estimated_memory_gb" on Interfaces has been '
                "deprecated as of nipype 1.0, please use Node.mem_gb."
            )
        if hasattr(self, '_mem_x'):
            if not isinstance(self._mem_x, tuple):
                self._mem_x = (self._mem_x, )
            if len(self._mem_x) == 1:
                return self._apply_mem_x()
            try:
                mem_x_path = getattr(self.inputs, self._mem_x[1])
            except AttributeError as e:
                raise AttributeError(
                    f'{e.args[0]} in Node \'{self.name}\'') from e
            if self._check_mem_x_path(mem_x_path):
                # constant + mem_x[0] * t
                return self._apply_mem_x(mem_x_path)
            raise FileNotFoundError(2, 'The memory estimate for Node '
                                    f"'{self.name}' depends on the input "
                                    f"'{self._mem_x[1]}' but no such file or "
                                    'directory', mem_x_path)
        return self._mem_gb

    def _check_mem_x_path(self, mem_x_path):
        '''Method to check if a supplied multiplicand path exists.

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
        multiplicand input

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

    def _apply_mem_x(self, multiplicand=None):
        '''Method to calculate and memoize a Node's estimated memory
        footprint.

        Parameters
        ----------
        multiplicand : str, int or None

        Returns
        -------
        number
            estimated memory usage (GB)
        '''
        if hasattr(self, '_mem_x'):
            if not isinstance(multiplicand, int):
                if self._check_mem_x_path(multiplicand):
                    multiplicand = get_data_size(
                        self._grab_first_path(multiplicand))
            self._mem_gb = self._mem_gb + self._mem_x[0] * multiplicand
            del self._mem_x
        return self._mem_gb

    @property
    def mem_x(self):
        """Get memory multiplier"""
        if hasattr(self, '_mem_x'):
            return self._mem_x
        return None


class MapNode(Node, pe.MapNode):
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
                            node.mem_x, tuple
                        ) and node.mem_x[1] == field:
                            input_resultfile = node.input_source.get(field)
                            if input_resultfile:
                                if isinstance(input_resultfile, tuple):
                                    input_resultfile = input_resultfile[0]
                                try:
                                    # memoize node._mem_gb if path
                                    # already exists
                                    multiplicand_path = _load_resultfile(
                                        input_resultfile
                                    ).inputs['in_file']
                                    node._apply_mem_x(multiplicand_path)
                                except FileNotFoundError:
                                    if hasattr(self, '_largest_func'):
                                        node._apply_mem_x(self._largest_func)
                                    else:
                                        # memoize the path otherwise
                                        node._mem_x = (node._mem_x[0],
                                                       multiplicand_path)


def get_data_size(filepath):
    """Function to return the size of a functional image (x * y * z * t)

    # !!!
    # Temporarily returns just the time dimension, defaulting to 1200.
    # !!!

    Parameters
    ----------
    filepath : str or path

    Returns
    -------
    int
    """
    return prod(load(filepath).shape)
