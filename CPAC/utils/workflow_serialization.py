from typing import Callable

from . import Configuration


def cpac_flowdump_serializer(
        flowdump_serializer: Callable[[object], object],
        obj: object
) -> object:
    """
    Custom flowdump serializer that removes `json_data` fields
    and CPAC configs from workflow json files as these are repeated
    for every node (and increase file size dramatically).
    """
    if isinstance(obj, dict):
        if 'json_data' in obj:
            obj_clone = obj.copy()
            obj_clone['json_data'] = '[truncated]'
            obj = obj_clone
        return flowdump_serializer(obj)
    if isinstance(obj, Configuration):
        return '[C-PAC config]'
    return flowdump_serializer(obj)
