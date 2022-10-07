# Copyright (C) 2022  C-PAC Developers

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
"""Guardrails to protect against bad registrations"""
from copy import deepcopy
from nipype.interfaces.ants import Registration
from nipype.interfaces.fsl import FLIRT
from nipype.interfaces.utility import Function
from CPAC.pipeline.nipype_pipeline_engine import Node, Workflow
from CPAC.pipeline.nipype_pipeline_engine.utils import connect_from_spec
from CPAC.qc import qc_masks, REGISTRATION_GUARDRAIL_THRESHOLDS


_SPEC_KEYS = {
    FLIRT: {'reference': 'reference', 'registered': 'out_file'},
    Registration: {'reference': 'reference', 'registered': 'out_file'}}


class BadRegistrationError(ValueError):
    """Exception for when a QC measure for a registration falls below a
    specified threshold"""
    def __init__(self, *args, metric=None, value=None, threshold=None,
                 **kwargs):
        """
        Parameters
        ----------
        metric : str
            QC metric

        value : float
            calculated QC value

        threshold : float
            specified threshold
        """
        msg = "Registration failed quality control"
        if all(arg is not None for arg in (metric, value, threshold)):
            msg += f" ({metric}: {value} < {threshold})"
        msg += "."
        super().__init__(msg, *args, **kwargs)


def registration_guardrail(registered: str, reference: str, retry: bool = False
                           ) -> str:
    """Check QC metrics post-registration and throw an exception if
    metrics are below given thresholds.

    If inputs point to images that are not masks, images will be
    binarized before being compared.

    .. seealso::

       :py:mod:`CPAC.qc.qcmetrics`
          Documentation of the :py:mod:`CPAC.qc.qcmetrics` module.

    Parameters
    ----------
    registered, reference : str
        path to mask

    retry : bool
        can retry?

    Returns
    -------
    registered_mask : str
        path to mask
    """
    qc_metrics = qc_masks(registered, reference)
    for metric, threshold in REGISTRATION_GUARDRAIL_THRESHOLDS.items():
        if threshold is not None:
            value = qc_metrics.get(metric)
            if isinstance(value, list):
                value = value[0]
            if value < threshold:
                with open(f'{registered}.failed_qc', 'w',
                          encoding='utf-8') as _f:
                    _f.write(True)
                if retry:
                    registered = f'{registered}-failed'
                else:
                    raise BadRegistrationError(metric=metric, value=value,
                                               threshold=threshold)
    return registered


def registration_guardrail_node(name=None):
    """Convenience method to get a new registration_guardrail Node

    Parameters
    ----------
    name : str, optional

    Returns
    -------
    Node
    """
    if name is None:
        name = 'registration_guardrail'
    return Node(Function(input_names=['registered',
                                      'reference'],
                         output_names=['registered'],
                         imports=['from CPAC.qc import qc_masks, '
                                  'REGISTRATION_GUARDRAIL_THRESHOLDS',
                                  'from CPAC.registration.guardrails '
                                  'import BadRegistrationError'],
                         function=registration_guardrail), name=name)


def registration_guardrail_workflow(registration_node, retry=True):
    """A workflow to handle hitting a registration guardrail

    Parameters
    ----------
    name : str

    registration_node : Node

    retry : bool, optional

    Returns
    -------
    Workflow
    """
    name = f'{registration_node.name}_guardrail'
    wf = Workflow(name=f'{name}_wf')
    outputspec = deepcopy(registration_node.outputs)
    guardrail = registration_guardrail_node(name)
    outkey = spec_key(registration_node, 'registered')
    wf.connect([
        (registration_node, guardrail, [
            (spec_key(registration_node, 'reference'), 'reference')]),
        (registration_node, guardrail, [(outkey, 'registered')])])
    if retry:
        wf = retry_registration(wf, registration_node,
                                guardrail.outputs.registered)
    else:
        wf.connect(guardrail, 'registered', outputspec, outkey)
        wf = connect_from_spec(wf, outputspec, registration_node, outkey)
    return wf


def retry_registration(wf, registration_node, registered):
    """Function conditionally retry registration if previous attempt failed

    Parameters
    ----------
    wf : Workflow

    registration_node : Node

    registered : str

    Returns
    -------
    Workflow
    """
    name = f'retry_{registration_node.name}'
    retry_node = Node(Function(function=retry_registration_node,
                               inputs=['registered', 'registration_node'],
                               outputs=['registered']), name=name)
    retry_node.inputs.registration_node = registration_node
    inputspec = registration_node.inputs
    outputspec = registration_node.outputs
    outkey = spec_key(registration_node, 'registered')
    guardrail = registration_guardrail_node(f'{name}_guardrail')
    wf = connect_from_spec(wf, inputspec, retry_node)
    wf.connect([
        (inputspec, guardrail, [
            (spec_key(retry_node, 'reference'), 'reference')]),
        (retry_node, guardrail, [(outkey, 'registered')]),
        (guardrail, outputspec, [('registered', outkey)])])
    wf = connect_from_spec(wf, retry_node, outputspec, registered)
    return wf


def retry_registration_node(registered, registration_node):
    """Retry registration if previous attempt failed

    Parameters
    ----------
    registered : str

    registration_node : Node

    Returns
    -------
    Node
    """
    from CPAC.pipeline.random_state.seed import MAX_SEED, random_seed
    seed = random_seed()
    if registered.endswith('-failed'):
        retry_node = registration_node.clone(
            name=f'{registration_node.name}-retry')
        if isinstance(seed, int):
            if seed < MAX_SEED:  # increment random seed
                retry_node.seed = seed + 1
            else:  # loop back to minumum seed
                retry_node.seed = 1
        return retry_node
    return registration_node


def spec_key(interface, guardrail_key):
    """Function to get the canonical key to connect to a guardrail

    Parameters
    ----------
    interface : Interface or Node

    guardrail_key : str

    Returns
    -------
    str
    """
    if isinstance(interface, Node):
        interface = interface.interface
    return _SPEC_KEYS.get(interface, {}).get(guardrail_key, guardrail_key)
