from typing import Callable

from nipype.interfaces.base.support import Bunch
from traits.has_traits import HasTraits

from . import Configuration


def _truncate_large_cpac_internals(obj: object) -> object:
    """Recursively replace json_data values and
    Config objects.
    """
    if isinstance(obj, Configuration):
        return '[C-PAC config]'
    if isinstance(obj, list):
        return [_truncate_large_cpac_internals(i) for i in obj]
    if isinstance(obj, tuple):
        return (_truncate_large_cpac_internals(i) for i in obj)
    if isinstance(obj, HasTraits):
        return _truncate_large_cpac_internals(obj.trait_get())
    if isinstance(obj, (Bunch, dict)):
        return {
            str(k): _truncate_large_cpac_internals(v)
            if k != "json_data" else
            "[Truncated]"
            for k, v in obj.items()
        }
    return obj


def cpac_flowdump_serializer(
        flowdump_serializer: Callable[[object], object],
        obj: object
) -> object:
    """
    Custom flowdump serializer that removes `json_data` fields
    and CPAC configs from workflow json files as these are repeated
    for every node (and increase file size dramatically).
    """
    return flowdump_serializer(_truncate_large_cpac_internals(obj))
