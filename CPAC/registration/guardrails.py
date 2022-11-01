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
from nipype.interfaces.base import TraitedSpec, traits
from nipype.interfaces.io import add_traits, IOBase
from nipype.interfaces.utility.base import MergeInputSpec
from CPAC.qc import qc_masks, registration_guardrail_thresholds
from CPAC.registration.exceptions import BadRegistrationError
from CPAC.registration.utils import hardcoded_reg
from CPAC.utils.docs import retry_docstring


class BestOfOutputSpec(TraitedSpec):
    """Outputspec for :py:class:`BestOf`"""
    index = traits.Int(desc="0-indexed index of minimum error")


class BestOf(IOBase):
    """Returns the index of the smallest 'error'. Inputs are 1-indexed
    and output is 0-indexed to mirror Merge and Select.

    .. seealso::

       nipype.interfaces.utility.base.Merge

       nipype.interfaces.utility.base.Select

    Example
    -------
    >>> best_of = BestOf(3)
    >>> best_of.inputs.error1 = 0.5
    >>> best_of.inputs.error2 = 0.1
    >>> best_of.inputs.error3 = 0.2
    >>> res = best_of.run()
    >>> res.outputs.index
    1
    """
    input_spec = MergeInputSpec
    output_spec = BestOfOutputSpec

    def __init__(self, numinputs=0, **inputs):
        super().__init__(**inputs)
        self._numinputs = numinputs
        if numinputs >= 1:
            input_names = [f"error{(i + 1)}" for i in range(numinputs)]
        else:
            input_names = []
        add_traits(self.inputs, input_names)

    def _getval(self, idx):
        return getattr(self.inputs, f"error{idx + 1}", 1)

    def _list_outputs(self):
        outputs = self._outputs().get()
        if self._numinputs >= 1:
            values = [self._getval(idx) for idx in range(self._numinputs)]
            outputs["index"] = values.index(min(values))
        return outputs


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
        outputs with ``retry_on_first_failure``

    error : float:
        sum of distance from thresholded QC metrics (min=0, max=count(metrics))
    """
    logger = logging.getLogger('nipype.workflow')
    qc_metrics = qc_masks(registered, reference)
    failed_qc = 0
    error = 0
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
            error += (1 - value)
    return registered, failed_qc, error


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


def skip_if_first_try_succeeds(interface):
    """Set an interface up to skip if a previous attempt succeeded"""
    if hasattr(interface, 'input_spec'):
        interface.inputs.add_trait('previous_failure', traits.Bool())
    return interface
