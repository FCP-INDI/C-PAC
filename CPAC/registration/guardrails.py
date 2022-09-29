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
from nipype.interfaces.utility import Function
from CPAC.pipeline.nipype_pipeline_engine import Node
from CPAC.qc import qc_masks, REGISTRATION_GUARDRAIL_THRESHOLDS


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


def registration_guardrail_node(name=None, ):
    """Convenience method to get a new registration_guardrail Node"""
    if name is None:
        name = 'registration_guardrail'
    return Node(Function(input_names=['registered_mask',
                                      'reference_mask'],
                         output_names=['registered_mask'],
                         function=registration_guardrail), name=name)


def registration_guardrail(registered: str, reference: str) -> str:
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

    Returns
    -------
    registered_mask : str
        path to mask
    """
    qc_metrics = qc_masks(registered, reference)
    for metric, threshold in REGISTRATION_GUARDRAIL_THRESHOLDS.items():
        if threshold is not None:
            value = qc_metrics.get(metric)
            if value < threshold:
                raise BadRegistrationError(metric=metric, value=value,
                                           threshold=threshold)
    return registered
