# Copyright (C) 2018-2025  C-PAC Developers

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
"""Specify the resources that C-PAC writes to the output direcotry."""

from importlib.resources import files
from typing import ClassVar

import pandas as pd


class Outputs:
    """Settle some things about the resource pool reference and the output directory."""

    reference_csv = str(files("CPAC").joinpath("resources/cpac_outputs.tsv"))

    try:
        reference = pd.read_csv(reference_csv, delimiter="\t", keep_default_na=False)
    except Exception as e:
        err = (
            "\n[!] Could not access or read the cpac_outputs.tsv "
            f"resource file:\n{reference_csv}\n\nError details {e}\n"
        )
        raise Exception(err)

    # all outputs
    any = list(reference.Resource)

    # extra outputs that we don't write to the output directory, unless the
    # user selects to do so
    debugging = list(reference[reference["Optional: Debugging"] == "Yes"]["Resource"])

    # functional data that are 4D time series, instead of derivatives
    functional_timeseries = list(
        reference[reference["4D Time Series"] == "Yes"]["Resource"]
    )

    anat: ClassVar[list[str]] = list(
        reference[reference["Sub-Directory"] == "anat"]["Resource"]
    )
    func: ClassVar[list[str]] = list(
        reference[reference["Sub-Directory"] == "func"]["Resource"]
    )

    # outputs to send into smoothing, if smoothing is enabled, and
    # outputs to write out if the user selects to write non-smoothed outputs
    _template_filter = reference["Space"] == "template"
    _epitemplate_filter = reference["Space"] == "EPI template"
    _symtemplate_filter = reference["Space"] == "symmetric template"
    _T1w_native_filter = reference["Space"] == "T1w"
    _bold_native_filter = reference["Space"] == "functional"
    _long_native_filter = reference["Space"] == "longitudinal T1w"
    _nonsmoothed_filter = reference["To Smooth"] == "Yes"
    _zstd_filter = reference["To z-std"] == "Yes"
    _corr_filter = reference["Type"] == "correlation"

    all_template_filter = _template_filter | _epitemplate_filter | _symtemplate_filter
    all_native_filter = _T1w_native_filter | _bold_native_filter | _long_native_filter

    bold_native: ClassVar[list[str]] = list(reference[_bold_native_filter]["Resource"])

    native_nonsmooth = list(
        reference[all_native_filter & _nonsmoothed_filter]["Resource"]
    )
    template_nonsmooth = list(
        reference[all_template_filter & _nonsmoothed_filter]["Resource"]
    )

    to_smooth = list(reference[_nonsmoothed_filter]["Resource"])
    to_zstd = list(reference[_zstd_filter & ~_corr_filter]["Resource"])
    to_fisherz = list(reference[_zstd_filter & _corr_filter]["Resource"])

    # don't write these, unless the user selects to write native-space outputs
    native_smooth = list(
        reference[~all_template_filter & ~_nonsmoothed_filter]["Resource"]
    )

    # ever used??? contains template-space, smoothed, both raw and z-scored
    template_smooth = list(
        reference[all_template_filter & ~_nonsmoothed_filter]["Resource"]
    )

    _bold_filter = reference["Type"] == "bold"
    _ts_filter = reference["4D Time Series"] == "Yes"
    bold_ts = list(
        reference[_bold_filter & _bold_native_filter & _ts_filter]["Resource"]
    )

    # outputs to send into z-scoring, if z-scoring is enabled, and
    # outputs to write out if user selects to write non-z-scored outputs
    native_raw = list(
        reference[all_native_filter & (reference["To z-std"] == "Yes")]["Resource"]
    )

    template_raw = list(
        reference[all_template_filter & (reference["To z-std"] == "Yes")]["Resource"]
    )

    def _is_cifti(_file_key):
        return _file_key.upper().startswith("CIFTI ")

    ciftis = reference[reference.File.map(_is_cifti)][["Resource", "File"]]
    ciftis = {
        cifti.Resource: cifti.File.split(" ")[-1]
        for cifti in ciftis.itertuples()
        if " " in cifti.File
    }

    def _is_gifti(_file_key):
        return _file_key.upper().startswith("GIFTI ")

    giftis = reference[reference.File.map(_is_gifti)][["Resource", "File"]]
    giftis = {
        gifti.Resource: gifti.File.split(" ")[-1]
        for gifti in giftis.itertuples()
        if " " in gifti.File
    }


def group_derivatives(pull_func: bool = False) -> list[str]:
    """Gather keys for anatomical and functional derivatives for group analysis."""
    derivatives: list[str] = Outputs.func + Outputs.anat
    if pull_func:
        derivatives = derivatives + Outputs.bold_native
    return derivatives
