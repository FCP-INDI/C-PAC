from typing import Callable, List, Union, Optional, Dict


class NodeBlockData(object):
    def __init__(
            self,
            func: Callable,
            name: Optional[str] = None,
            config: Optional[List[str]] = None,
            switch: Optional[Union[List[str], List[List[str]]]] = None,
            option_key: Optional[Union[str, List[str]]] = None,
            option_val: Optional[Union[str, List[str]]] = None,
            inputs: Optional[List[Union[str, list, tuple]]] = None,
            outputs: Optional[Union[List[str], Dict[str]]] = None
    ) -> None:
        self.func = func
        self.name: Optional[str] = name
        self.config: Optional[List[str]] = config
        self.switch: Optional[Union[List[str], List[List[str]]]] = switch
        self.option_key: Optional[Union[str, List[str]]] = option_key
        self.option_val: Optional[Union[str, List[str]]] = option_val
        self.inputs: Optional[List[Union[str, list, tuple]]] = inputs
        self.outputs: Optional[Union[List[str], Dict[str]]] = outputs

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def legacy_nodeblock_dict(self):
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
        outputs: Optional[Union[List[str], Dict[str]]] = None
):
    return lambda func: NodeBlockData(func, name, config, switch, option_key, option_val, inputs, outputs)
