from typing import Callable, List, Union, Optional, Dict, Any


class NodeBlockFunction(object):
    """
    Stores a reference to the nodeblock function and all of its meta-data.
    """

    def __init__(
            self,
            func: Callable,
            name: Optional[str] = None,
            config: Optional[List[str]] = None,
            switch: Optional[Union[List[str], List[List[str]]]] = None,
            option_key: Optional[Union[str, List[str]]] = None,
            option_val: Optional[Union[str, List[str]]] = None,
            inputs: Optional[List[Union[str, list, tuple]]] = None,
            outputs: Optional[Union[List[str], Dict[str, Any]]] = None
    ) -> None:
        self.func = func
        """Nodeblock function reference."""
        self.name: Optional[str] = name
        """Used in the graph and logging to identify the NodeBlock and its component nodes."""
        self.config: Optional[List[str]] = config
        """
        Indicates the nested keys in a C-PAC pipeline configuration should configure a NodeBlock built from this
        function. If config is sct to "None", then all other configuration-related entities must be specified from the
        root of the configuration.
        """
        self.switch: Optional[Union[List[str], List[List[str]]]] = switch
        """
        Indicates any keys that should evaluate to True for this NodeBlock to be active. A list of lists of strings
        indicates multiple switches that must all be True to run, and is currently only an option if config is set to
        "None".
        """
        self.option_key: Optional[Union[str, List[str]]] = option_key
        """
        Indicates the nested keys (starting at the nested key indicated by config) that should configure this NodeBlock.
        """
        self.option_val: Optional[Union[str, List[str]]] = option_val
        """Indicates values for which this NodeBlock should be active."""
        self.inputs: Optional[List[Union[str, list, tuple]]] = inputs
        """ResourcePool keys indicating files needed for the NodeBlock's functionality."""
        self.outputs: Optional[Union[List[str], Dict[str, Any]]] = outputs
        """
        ResourcePool keys indicating files generated or updated by the NodeBlock, optionally including metadata
        for the outputs' respective sidecars.
        """

    # all node block functions have this signature
    def __call__(self, wf, cfg, strat_pool, pipe_num, opt=None):
        return self.func(wf, cfg, strat_pool, pipe_num, opt)

    def legacy_nodeblock_dict(self):
        """
        Returns nodeblock metadata as a dictionary. Helper for compatibility reasons.
        """
        return {
            'name': self.name,
            'config': self.config,
            'switch': self.switch,
            'option_key': self.option_key,
            'option_val': self.option_val,
            'inputs': self.inputs,
            'outputs': self.outputs
        }


def nodeblock(
        name: Optional[str] = None,
        config: Optional[List[str]] = None,
        switch: Optional[Union[List[str], List[List[str]]]] = None,
        option_key: Optional[Union[str, List[str]]] = None,
        option_val: Optional[Union[str, List[str]]] = None,
        inputs: Optional[List[Union[str, list, tuple]]] = None,
        outputs: Optional[Union[List[str], Dict[str, Any]]] = None
):
    """
    Define a node block: Connections to the pipeline configuration and to other node blocks.

    Parameters
    ----------
    name
        Used in the graph and logging to identify the NodeBlock and its component nodes.
    config
        Indicates the nested keys in a C-PAC pipeline configuration should configure a NodeBlock built from this
        function. If config is sct to "None", then all other configuration-related entities must be specified from the
        root of the configuration.
    switch
        Indicates any keys that should evaluate to True for this NodeBlock to be active. A list of lists of strings
        indicates multiple switches that must all be True to run, and is currently only an option if config is set to
        "None".
    option_key
        Indicates the nested keys (starting at the nested key indicated by config) that should configure this NodeBlock.
    option_val
        Indicates values for which this NodeBlock should be active.
    inputs
        ResourcePool keys indicating files needed for the NodeBlock's functionality.
    outputs
        ResourcePool keys indicating files generated or updated by the NodeBlock, optionally including metadata
        for the outputs' respective sidecars.

    Returns
    -------

    """
    return lambda func: NodeBlockFunction(
        func,
        name if name is not None else func.__name__,
        config,
        switch,
        option_key,
        option_val,
        inputs,
        outputs
    )
