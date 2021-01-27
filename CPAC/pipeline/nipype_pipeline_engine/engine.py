'''Module to import Nipype Pipeline engine and override some Classes.
See https://fcp-indi.github.com/docs/developer/nodes
for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.pipeline.engine.html
for Nipype's documentation.'''  # noqa E501
from nipype.pipeline import engine as pe
from functools import partialmethod

# set global default mem_gb
DEFAULT_MEM_GB = 2.0


def _doctest_skiplines(docstring, lines_to_skip):
    '''
    Function to add '  # doctest: +SKIP' to the end of docstring lines
    to skip in imported docstrings.

    Parameters
    ----------
    docstring: str

    lines_to_skip: set or list

    Returns
    -------
    docstring: str

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

    __init__ = partialmethod(pe.Node.__init__, mem_gb=DEFAULT_MEM_GB)


class MapNode(pe.MapNode):
    __doc__ = _doctest_skiplines(
        f'mem_gb={DEFAULT_MEM_GB}\n\n{pe.MapNode.__doc__}',
        {"    ...                           'functional3.nii']"}
    )

    __init__ = partialmethod(pe.MapNode.__init__, mem_gb=DEFAULT_MEM_GB)
