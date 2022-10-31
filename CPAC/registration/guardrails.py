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
import logging
from typing import Tuple
from nipype.interfaces.utility import Function, Merge, Select
# pylint: disable=unused-import
from CPAC.pipeline.nipype_pipeline_engine import Node, Workflow
from CPAC.pipeline.random_state.seed import increment_seed
from CPAC.qc import qc_masks, registration_guardrail_thresholds
from CPAC.registration.exceptions import BadRegistrationError
from CPAC.registration.utils import hardcoded_reg
from CPAC.utils.docs import retry_docstring


# noqa: F401
def guardrail_selection(wf: 'Workflow', node1: 'Node', node2: 'Node',
                        output_key: str = 'registered',
                        guardrail_node: 'Node' = None) -> Node:
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
    wf : Workflow

    node1, node2 : Node
        first try, retry

    output_key : str
        field to choose

    guardrail_node : Node
        guardrail to collect 'failed_qc' from if not node1

    Returns
    -------
    select : Node
    """
    # pylint: disable=redefined-outer-name,reimported,unused-import
    from CPAC.pipeline.nipype_pipeline_engine import Node, Workflow
    if guardrail_node is None:
        guardrail_node = node1
    name = node1.name
    if output_key != 'registered':
        name = f'{name}_{output_key}'
    choices = Node(Merge(2), run_without_submitting=True,
                   name=f'{name}_choices')
    select = Node(Select(), run_without_submitting=True,
                  name=f'choose_{name}')
    wf.connect([(node1, choices, [(output_key, 'in1')]),
                (node2, choices, [(output_key, 'in2')]),
                (choices, select, [('out', 'inlist')]),
                (guardrail_node, select, [('failed_qc', 'index')])])
    return select


def registration_guardrail(registered: str, reference: str,
                           retry: bool = False, retry_num: int = 0
                           ) -> Tuple[str, int]:
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

    retry : bool, optional
        can retry?

    retry_num : int, optional
        how many previous tries?

    Returns
    -------
    registered_mask : str
        path to mask

    failed_qc : int
        metrics met specified thresholds?, used as index for selecting
        outputs
        .. seealso::

           :py:mod:`guardrail_selection`
    """
    logger = logging.getLogger('nipype.workflow')
    qc_metrics = qc_masks(registered, reference)
    failed_qc = 0
    for metric, threshold in registration_guardrail_thresholds().items():
        if threshold is not None:
            value = qc_metrics.get(metric)
            if isinstance(value, list):
                value = value[0]
            if value < threshold:
                failed_qc = 1
                with open(f'{registered}.failed_qc', 'w',
                          encoding='utf-8') as _f:
                    _f.write(f'{metric}: {value} < {threshold}')
                if retry:
                    registered = f'{registered}-failed'
                else:
                    bad_registration = BadRegistrationError(
                        metric=metric, value=value, threshold=threshold)
                    logger.error(str(bad_registration))
                    if retry_num:
                        # if we've already retried, raise the error
                        raise bad_registration
    return registered, failed_qc


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
                         output_names=['registered', 'failed_qc'],
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


def retry_clone(node: 'Node') -> 'Node':
    """Function to clone a node, name the clone, and increment its
    random seed

    Parameters
    ----------
    node : Node

    Returns
    -------
    Node
    """
    return increment_seed(node.clone(f'retry_{node.name}'))


# pylint: disable=missing-function-docstring,too-many-arguments
@retry_docstring(hardcoded_reg)
def retry_hardcoded_reg(moving_brain, reference_brain, moving_skull,
                        reference_skull, ants_para, moving_mask=None,
                        reference_mask=None, fixed_image_mask=None,
                        interp=None, reg_with_skull=0, previous_failure=False):
    if not previous_failure:
        return [], None
    return hardcoded_reg(moving_brain, reference_brain, moving_skull,
                         reference_skull, ants_para, moving_mask,
                         reference_mask, fixed_image_mask, interp,
                         reg_with_skull)
