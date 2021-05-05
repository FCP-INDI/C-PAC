'''Module to import Nipype Pipeline engine and override some Classes.
See https://fcp-indi.github.com/docs/developer/nodes
for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.pipeline.engine.html
for Nipype's documentation.'''  # noqa E501
import os
import re
from inspect import Parameter, Signature, signature
from nipype.pipeline import engine as pe
from traits.trait_base import Undefined

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

    __init__.__signature__ = Signature(parameters=[
        p[1] if p[0] != 'mem_gb' else (
            'mem_gb',
            Parameter('mem_gb', Parameter.POSITIONAL_OR_KEYWORD,
                      default=DEFAULT_MEM_GB)
        )[1] for p in signature(pe.Node).parameters.items()])

    __init__.__doc__ = re.sub(r'(?<!\s):', ' :', '\n'.join([
        pe.Node.__init__.__doc__.rstrip(),
        '''
        mem_gb : int or float
            Estimate (in GB) of constant memory to allocate for this node.

        mem_x : tuple
            (int or float, str)
            Multiplier for memory allocation such that
            `mem_x[0]` times
            the number of timepoints in file at `mem+x[1]` plus
            `mem_gb` equals
            the total memory allocation for the node.
            (TEMPORARY: will replace number of timepoints with
            spatial dimensions times timepoints)''']))

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
            try:
                mem_x_path = getattr(self.inputs, self._mem_x[1])
            except AttributeError as e:
                raise AttributeError(
                    f'{e.args[0]} in Node \'{self.name}\'') from e
            if self._check_mem_x_path(mem_x_path):
                # constant + mem_x[0] * t
                return self._apply_mem_x(mem_x_path)
            else:
                mem_x_path = self.input_source[self._mem_x[1]]
                if self._check_mem_x_path(mem_x_path):
                    # constant + mem_x[0] * t
                    return self._apply_mem_x(mem_x_path)
                else:
                    # constant + mem_x[0] * 300
                    return self._mem_gb + self._mem_x[0] * 300

        return self._mem_gb

    def _check_mem_x_path(self, mem_x_path):
        if isinstance(mem_x_path, list):
            mem_x_path = mem_x_path[0] if len(mem_x_path) else Undefined
        try:
            return mem_x_path is not Undefined and os.path.exists(
                mem_x_path)
        except TypeError:
            return False

    def _apply_mem_x(self, mem_x_path):
        from CPAC.vmhc.utils import get_img_nvols
        print(f'mem_x_path: {mem_x_path}')
        if os.path.exists(mem_x_path):
            self._mem_gb = self._mem_gb + self._mem_x[0] * get_img_nvols(
                mem_x_path)
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
                    if node and hasattr(
                        node, 'mem_x'
                    ) and isinstance(
                        node.mem_x, tuple
                    ) and node.mem_x[1] == field:
                        node._apply_mem_x(node.input_source[field][0])
