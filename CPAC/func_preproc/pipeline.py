# Copyright (C) 2025  C-PAC Developers

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
"""Build a functional preprocessing pipeline."""

from typing import cast

from CPAC.distortion_correction.distortion_correction import (
    distcor_blip_afni_qwarp,
    distcor_blip_fsl_topup,
    distcor_phasediff_fsl_fugue,
)
from CPAC.func_preproc.func_motion import (
    stack_motion_blocks,
)
from CPAC.func_preproc.func_preproc import (
    bold_mask_afni,
    bold_mask_anatomical_based,
    bold_mask_anatomical_refined,
    bold_mask_ccs,
    bold_mask_fsl,
    bold_mask_fsl_afni,
    bold_masking,
    func_despike,
    func_mean,
    func_normalize,
    func_reorient,
    func_scaling,
    func_slice_time,
    func_truncate,
)
from CPAC.pipeline.engine import ResourcePool
from CPAC.pipeline.nodeblock import NODEBLOCK_STACK, NodeBlockFunction
from CPAC.registration.registration import (
    coregistration_prep_fmriprep,
    coregistration_prep_mean,
    coregistration_prep_vol,
    mask_sbref,
)
from CPAC.utils.configuration import Configuration


def stack_func_preproc_blocks(
    sub_dict, cfg: Configuration, rpool: ResourcePool
) -> NODEBLOCK_STACK:
    """Stack functional preprocessing nodeblocks."""
    func_blocks: dict[str, NODEBLOCK_STACK] = {}
    func_blocks["init"] = [func_reorient, func_scaling, func_truncate]
    func_blocks["preproc"] = [func_despike, func_slice_time]
    if not rpool.check_rpool("desc-mean_bold"):
        func_blocks["preproc"].append(func_mean)
    func_blocks["mask"] = (
        []
        if rpool.check_rpool("space-bold_desc-brain_mask")
        else [
            [
                bold_mask_afni,
                bold_mask_fsl,
                bold_mask_fsl_afni,
                bold_mask_anatomical_refined,
                bold_mask_anatomical_based,
                bold_mask_ccs,
            ],
            bold_masking,
        ]
    )
    func_blocks["prep"] = [
        func_normalize,
        [
            coregistration_prep_vol,
            coregistration_prep_mean,
            coregistration_prep_fmriprep,
        ],
        mask_sbref,
    ]
    distcor_blocks: NODEBLOCK_STACK = []
    if "fmap" in sub_dict:
        fmap_keys = sub_dict["fmap"]
        if "phasediff" in fmap_keys or "phase1" in fmap_keys:
            if "magnitude" in fmap_keys or "magnitude1" in fmap_keys:
                distcor_blocks.append(distcor_phasediff_fsl_fugue)
        if len(fmap_keys) == 2:  # noqa: PLR2004
            for key in fmap_keys:
                if "epi_" not in key:
                    break
            else:
                distcor_blocks.append(distcor_blip_afni_qwarp)
                distcor_blocks.append(distcor_blip_fsl_topup)

    if distcor_blocks:
        if len(distcor_blocks) > 1:
            distcor_blocks = [cast(NodeBlockFunction, distcor_blocks)]
        func_blocks["prep"] += distcor_blocks

    return stack_motion_blocks(func_blocks, cfg)
