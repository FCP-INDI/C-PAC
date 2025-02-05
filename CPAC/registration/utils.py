# Copyright (C) 2014-2025  C-PAC Developers

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
"""Utilities for registration."""

import os
import subprocess
from typing import overload

import numpy as np
from voluptuous import RequiredFieldInvalid
from nipype.interfaces.ants import ApplyTransforms as AntsApplyTransforms
from nipype.interfaces.utility import Merge
from nipype.pipeline.engine import Node as NipypeNode, Workflow

from CPAC.pipeline.nipype_pipeline_engine import Node
from CPAC.utils.interfaces import Function


def single_ants_xfm_to_list(transform):
    """Convert a single ANTs transform to a list."""
    return [transform]


def interpolation_string(interpolation, reg_tool):
    """Translate interpolation string to reg_tool-specific string."""
    if reg_tool == "ants":
        pass
    elif reg_tool == "fsl":
        # translate to FSL
        # warning: flirt requires 'nearestneighbour', but FSL applywarp uses
        #          'nn', so this is designed for applywarp, as all FSL xfm's
        #          in C-PAC are now converted to .nii.gz
        interpolation = interpolation.replace("NearestNeighbor", "nn")
    return interpolation


def combine_inputs_into_list(input1, input2, input3):
    """Combine inputs into a list."""
    return [input1, input2, input3]


def seperate_warps_list(warp_list, selection):
    """Select the warp from the warp list."""
    selected_warp = None
    for warp in warp_list:
        if selection == "Warp":
            if "3Warp" in warp or "2Warp" in warp or "1Warp" in warp:
                selected_warp = warp
        elif selection in warp:
            selected_warp = warp
    return selected_warp


def check_transforms(transform_list):
    """Check if the transform list is empty."""
    transform_number = list(filter(None, transform_list))
    return [(transform_number[index]) for index in range(len(transform_number))], len(
        transform_number
    )


def generate_inverse_transform_flags(transform_list):
    """List whether each transform has an inverse."""
    inverse_transform_flags = []
    for transform in transform_list:
        # check `blip_warp_inverse` file name and rename it
        if "WARPINV" in transform:
            inverse_transform_flags.append(False)
        if "updated_affine" in transform:
            inverse_transform_flags.append(True)
        if "Initial" in transform:
            inverse_transform_flags.append(True)
        if "Rigid" in transform:
            inverse_transform_flags.append(True)
        if "Affine" in transform:
            inverse_transform_flags.append(True)
        if "InverseWarp" in transform:
            inverse_transform_flags.append(False)
    return inverse_transform_flags


def hardcoded_reg(
    moving_brain,
    reference_brain,
    moving_skull,
    reference_skull,
    ants_para,
    moving_mask=None,
    reference_mask=None,
    fixed_image_mask=None,
    interp=None,
    reg_with_skull=0,
):
    """Run ANTs registration."""
    # TODO: expand transforms to cover all in ANTs para

    regcmd = ["antsRegistration"]
    for para_index in range(len(ants_para)):
        for para_type in ants_para[para_index]:
            if para_type == "dimensionality":
                if ants_para[para_index][para_type] not in [2, 3, 4]:
                    err_msg = (
                        "Dimensionality specified in ANTs parameters:"
                        f" {ants_para[para_index][para_type]}, is not supported."
                        " Change to 2, 3, or 4 and try again"
                    )
                    raise ValueError(err_msg)
                regcmd.append("--dimensionality")
                regcmd.append(str(ants_para[para_index][para_type]))

            elif para_type == "verbose":
                if ants_para[para_index][para_type] not in [0, 1]:
                    err_msg = (
                        "Verbose output option in ANTs parameters:"
                        f" {ants_para[para_index][para_type]}, is not supported."
                        " Change to 0 or 1 and try again"
                    )
                    raise ValueError(err_msg)
                regcmd.append("--verbose")
                regcmd.append(str(ants_para[para_index][para_type]))

            elif para_type == "float":
                if ants_para[para_index][para_type] not in [0, 1]:
                    err_msg = (
                        "Float option in ANTs parameters:"
                        f" {ants_para[para_index][para_type]}, is not supported."
                        " Change to 0 or 1 and try again"
                    )
                    raise ValueError(err_msg)
                regcmd.append("--float")
                regcmd.append(str(ants_para[para_index][para_type]))

            elif para_type == "collapse-output-transforms":
                if ants_para[para_index][para_type] not in [0, 1]:
                    err_msg = (
                        "collapse-output-transforms specified in ANTs parameters:"
                        f" {ants_para[para_index][para_type]}, is not supported."
                        " Change to 0 or 1 and try again"
                    )
                    raise ValueError(err_msg)
                regcmd.append("--collapse-output-transforms")
                regcmd.append(str(ants_para[para_index][para_type]))

            elif para_type == "winsorize-image-intensities":
                if (
                    ants_para[para_index][para_type]["lowerQuantile"] is None
                    or ants_para[para_index][para_type]["upperQuantile"] is None
                ):
                    err_msg = (
                        "Please specifiy lowerQuantile and upperQuantile of ANTs"
                        " parameters --winsorize-image-intensities in pipeline config."
                    )
                    raise RequiredFieldInvalid(err_msg)
                regcmd.append("--winsorize-image-intensities")
                _quantile = ants_para[para_index][para_type]
                regcmd.append(
                    f"[{_quantile['lowerQuantile']},{_quantile['upperQuantile']}]"
                )

            elif para_type == "initial-moving-transform":
                if ants_para[para_index][para_type]["initializationFeature"] is None:
                    err_msg = (
                        "Please specifiy initializationFeature of ANTs parameters in"
                        " pipeline config."
                    )
                    raise RequiredFieldInvalid(err_msg)
                regcmd.append("--initial-moving-transform")
                initialization_feature = ants_para[para_index][para_type][
                    "initializationFeature"
                ]
                if reg_with_skull == 1:
                    regcmd.append(
                        f"[{reference_skull},{moving_skull},{initialization_feature}]"
                    )
                else:
                    regcmd.append(
                        f"[{reference_brain},{moving_brain},{initialization_feature}]"
                    )

            elif para_type == "transforms":
                for trans_index in range(len(ants_para[para_index][para_type])):
                    for trans_type in ants_para[para_index][para_type][trans_index]:
                        regcmd.append("--transform")
                        if trans_type in ("Rigid", "Affine"):
                            if (
                                ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["gradientStep"]
                                is None
                            ):
                                err_msg = (
                                    f"Please specifiy {trans_type} Gradient Step of"
                                    " ANTs parameters in pipeline config."
                                )
                                raise RequiredFieldInvalid(err_msg)
                            gradient_step = ants_para[para_index][para_type][
                                trans_index
                            ][trans_type]["gradientStep"]
                            regcmd.append(f"{trans_type}[{gradient_step}]")

                        if trans_type == "SyN":
                            if (
                                ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["gradientStep"]
                                is None
                            ):
                                err_msg = f"Please specifiy {trans_type} Gradient Step of ANTs parameters in pipeline config."
                                raise RequiredFieldInvalid(err_msg)
                            SyN_para = []
                            SyN_para.append(
                                str(
                                    ants_para[para_index][para_type][trans_index][
                                        trans_type
                                    ]["gradientStep"]
                                )
                            )
                            if (
                                ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["updateFieldVarianceInVoxelSpace"]
                                is not None
                            ):
                                SyN_para.append(
                                    str(
                                        ants_para[para_index][para_type][trans_index][
                                            trans_type
                                        ]["updateFieldVarianceInVoxelSpace"]
                                    )
                                )
                            if (
                                ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["totalFieldVarianceInVoxelSpace"]
                                is not None
                            ):
                                SyN_para.append(
                                    str(
                                        ants_para[para_index][para_type][trans_index][
                                            trans_type
                                        ]["totalFieldVarianceInVoxelSpace"]
                                    )
                                )
                            SyN_para = ",".join([str(elem) for elem in SyN_para])
                            regcmd.append(f"{trans_type}[{SyN_para}]")

                        if (
                            ants_para[para_index][para_type][trans_index][trans_type][
                                "metric"
                            ]["type"]
                            == "MI"
                        ):
                            if (
                                ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["metric"]["metricWeight"]
                                is None
                                or ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["metric"]["numberOfBins"]
                                is None
                            ):
                                err_msg = (
                                    "Please specifiy metricWeight and numberOfBins for"
                                    " metric MI of ANTs parameters in pipeline config."
                                )
                                raise RequiredFieldInvalid(err_msg)
                            MI_para = []
                            _metric = ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["metric"]
                            MI_para.append(
                                f"{_metric['metricWeight']},{_metric['numberOfBins']}"
                            )
                            if "samplingStrategy" in ants_para[para_index][para_type][
                                trans_index
                            ][trans_type]["metric"] and ants_para[para_index][
                                para_type
                            ][trans_index][trans_type]["metric"][
                                "samplingStrategy"
                            ] in ["None", "Regular", "Random"]:
                                MI_para.append(
                                    str(
                                        ants_para[para_index][para_type][trans_index][
                                            trans_type
                                        ]["metric"]["samplingStrategy"]
                                    )
                                )
                            if (
                                "samplingPercentage"
                                in ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["metric"]
                                and ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["metric"]["samplingPercentage"]
                                is not None
                            ):
                                MI_para.append(
                                    str(
                                        ants_para[para_index][para_type][trans_index][
                                            trans_type
                                        ]["metric"]["samplingPercentage"]
                                    )
                                )
                            MI_para = ",".join([str(elem) for elem in MI_para])
                            regcmd.append("--metric")
                            if reg_with_skull == 1:
                                regcmd.append(
                                    f"MI[{reference_skull},{moving_skull},{MI_para}]"
                                )
                            else:
                                regcmd.append(
                                    f"MI[{reference_brain},{moving_brain},{MI_para}]"
                                )

                        if (
                            ants_para[para_index][para_type][trans_index][trans_type][
                                "metric"
                            ]["type"]
                            == "CC"
                        ):
                            if (
                                ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["metric"]["metricWeight"]
                                is None
                                or ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["metric"]["radius"]
                                is None
                            ):
                                err_msg = (
                                    "Please specifiy metricWeight and radius for metric"
                                    " CC of ANTs parameters in pipeline config."
                                )
                                raise RequiredFieldInvalid(err_msg)
                            CC_para = []
                            _metric = ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["metric"]
                            CC_para.append(
                                f"{_metric['metricWeight']},{_metric['radius']}"
                            )
                            if "samplingStrategy" in ants_para[para_index][para_type][
                                trans_index
                            ][trans_type]["metric"] and ants_para[para_index][
                                para_type
                            ][trans_index][trans_type]["metric"][
                                "samplingStrategy"
                            ] in ["None", "Regular", "Random"]:
                                CC_para.append(
                                    str(
                                        ants_para[para_index][para_type][trans_index][
                                            trans_type
                                        ]["metric"]["samplingStrategy"]
                                    )
                                )
                            if (
                                "samplingPercentage"
                                in ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["metric"]
                                and ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["metric"]["samplingPercentage"]
                                is not None
                            ):
                                CC_para.append(
                                    str(
                                        ants_para[para_index][para_type][trans_index][
                                            trans_type
                                        ]["metric"]["samplingPercentage"]
                                    )
                                )
                            CC_para = ",".join([str(elem) for elem in CC_para])
                            regcmd.append("--metric")
                            regcmd.append(
                                f"CC[{reference_skull},{moving_skull},{CC_para}]"
                            )

                        if (
                            "convergence"
                            in ants_para[para_index][para_type][trans_index][trans_type]
                        ):
                            convergence_para = []
                            if (
                                ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["convergence"]["iteration"]
                                is None
                            ):
                                err_msg = (
                                    "Please specifiy convergence iteration of ANTs"
                                    " parameters in pipeline config."
                                )
                                raise RequiredFieldInvalid(err_msg)
                            convergence_para.append(
                                str(
                                    ants_para[para_index][para_type][trans_index][
                                        trans_type
                                    ]["convergence"]["iteration"]
                                )
                            )
                            if (
                                "convergenceThreshold"
                                in ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["convergence"]
                                and ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["convergence"]["convergenceThreshold"]
                                is not None
                            ):
                                convergence_para.append(
                                    str(
                                        ants_para[para_index][para_type][trans_index][
                                            trans_type
                                        ]["convergence"]["convergenceThreshold"]
                                    )
                                )
                            if (
                                "convergenceWindowSize"
                                in ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["convergence"]
                                and ants_para[para_index][para_type][trans_index][
                                    trans_type
                                ]["convergence"]["convergenceWindowSize"]
                                is not None
                            ):
                                convergence_para.append(
                                    str(
                                        ants_para[para_index][para_type][trans_index][
                                            trans_type
                                        ]["convergence"]["convergenceWindowSize"]
                                    )
                                )
                            convergence_para = ",".join(
                                [str(elem) for elem in convergence_para]
                            )
                            regcmd.append("--convergence")
                            regcmd.append(f"[{convergence_para}]")

                        if (
                            "smoothing-sigmas"
                            in ants_para[para_index][para_type][trans_index][trans_type]
                            and ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["smoothing-sigmas"]
                            is not None
                        ):
                            regcmd.append("--smoothing-sigmas")
                            regcmd.append(
                                str(
                                    ants_para[para_index][para_type][trans_index][
                                        trans_type
                                    ]["smoothing-sigmas"]
                                )
                            )

                        if (
                            "shrink-factors"
                            in ants_para[para_index][para_type][trans_index][trans_type]
                            and ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["shrink-factors"]
                            is not None
                        ):
                            regcmd.append("--shrink-factors")
                            regcmd.append(
                                str(
                                    ants_para[para_index][para_type][trans_index][
                                        trans_type
                                    ]["shrink-factors"]
                                )
                            )

                        if (
                            "use-histogram-matching"
                            in ants_para[para_index][para_type][trans_index][trans_type]
                        ):
                            if ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["use-histogram-matching"]:
                                regcmd.append("--use-histogram-matching")
                                regcmd.append("1")

                        if (
                            "winsorize-image-intensities"
                            in ants_para[para_index][para_type][trans_index][trans_type]
                            and ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["winsorize-image-intensities"]["lowerQuantile"]
                            is not None
                            and ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["winsorize-image-intensities"]["upperQuantile"]
                            is not None
                        ):
                            regcmd.append("--winsorize-image-intensities")
                            _quantile = ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["winsorize-image-intensities"]
                            regcmd.append(
                                f"[{_quantile['lowerQuantile']},{_quantile['upperQuantile']}]"
                            )

                        if (
                            "masks"
                            in ants_para[para_index][para_type][trans_index][trans_type]
                            and ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["masks"]
                            is not None
                        ):
                            if ants_para[para_index][para_type][trans_index][
                                trans_type
                            ]["masks"]:
                                regcmd.append("--masks")
                                regcmd.append(f"[{reference_mask},{moving_mask}]")
                            else:
                                regcmd.append("--masks")
                                regcmd.append("[NULL,NULL]")

            elif para_type == "masks":
                # lesion preproc has
                if fixed_image_mask is not None:
                    regcmd.append("--masks")
                    regcmd.append(str(fixed_image_mask))
                else:
                    if (
                        not ants_para[para_index][para_type]["fixed_image_mask"]
                        and ants_para[para_index][para_type]["moving_image_mask"]
                    ):
                        err_msg = (
                            "Masks option in ANTs parameters:"
                            f" {ants_para[para_index][para_type]} is not supported."
                            " Please set `fixed_image_mask` as True. Or set both"
                            " `fixed_image_mask` and `moving_image_mask` as False"
                        )
                        raise NotImplementedError(err_msg)
                    if (
                        ants_para[para_index][para_type]["fixed_image_mask"]
                        and ants_para[para_index][para_type]["moving_image_mask"]
                    ):
                        regcmd.append("--masks")
                        regcmd.append(
                            "[" + str(reference_mask) + "," + str(moving_mask) + "]"
                        )
                    elif (
                        ants_para[para_index][para_type]["fixed_image_mask"]
                        and ants_para[para_index][para_type]["moving_image_mask"]
                    ):
                        regcmd.append("--masks")
                        regcmd.append("[" + str(reference_mask) + "]")
                    else:
                        continue

    if interp is not None:
        regcmd.append("--interpolation")
        regcmd.append(f"{interp}")

    regcmd.append("--output")
    regcmd.append("[transform,transform_Warped.nii.gz]")

    # write out the actual command-line entry for testing/validation later
    command_file = os.path.join(os.getcwd(), "command.txt")
    with open(command_file, "wt") as f:
        f.write(" ".join(regcmd))

    try:
        subprocess.check_output(regcmd)
    except Exception as e:
        msg = (
            "[!] ANTS registration did not complete successfully."
            f"\n\nError details:\n{e}\n{e.output}\n"
        )
        raise RuntimeError(msg)

    warp_list = []
    warped_image = None

    files = [f for f in os.listdir(".") if os.path.isfile(f)]

    for f in files:
        if ("transform" in f) and ("Warped" not in f):
            warp_list.append(os.getcwd() + "/" + f)
        if "Warped" in f:
            warped_image = os.getcwd() + "/" + f

    if not warped_image:
        msg = (
            "\n\n[!] No registration output file found. ANTS registration may not have"
            " completed successfully.\n\n"
        )
        raise RuntimeError(msg)

    return warp_list, warped_image


def change_itk_transform_type(input_affine_file):
    """Produce an updated affine file for ANTs compatibility.

    This function takes in the affine.txt produced by the c3d_affine_tool
    (which converted an FSL FLIRT affine.mat into the affine.txt).

    It then modifies the 'Transform Type' of this affine.txt so that it is
    compatible with the antsApplyTransforms tool and produces a new affine
    file titled 'updated_affine.txt'
    """
    new_file_lines = []

    with open(input_affine_file) as f:
        for line in f:
            if "Transform:" in line:
                if "MatrixOffsetTransformBase_double_3_3" in line:
                    transform_line = "Transform: AffineTransform_double_3_3\n"
                    new_file_lines.append(transform_line)
            else:
                new_file_lines.append(line)

    updated_affine_file = os.path.join(os.getcwd(), "updated_affine.txt")

    with open(updated_affine_file, "wt") as f:
        for line in new_file_lines:
            f.write(line)

    return updated_affine_file


def one_d_to_mat(one_d_filename):
    """Convert a .1D file to a .mat directory.

    Parameters
    ----------
    one_d_filename : str
        The filename of the .1D file to convert

    Returns
    -------
    mat_filenames : list of str
        The of paths in the .mat directory created
    """
    mat_dirname = one_d_filename.replace(".1D", ".mat")
    with open(one_d_filename, "r") as one_d_file:
        rows = [
            np.reshape(row, (4, 4)).astype("float")
            for row in [
                [term.strip() for term in row.split(" ") if term.strip()] + [0, 0, 0, 1]
                for row in [
                    line.strip()
                    for line in one_d_file.readlines()
                    if not line.startswith("#")
                ]
            ]
        ]
    try:
        os.mkdir(mat_dirname)
    except FileExistsError:
        pass
    for i, row in enumerate(rows):
        np.savetxt(
            os.path.join(mat_dirname, f"MAT_{i:04}"), row, fmt="%.5f", delimiter=" "
        )
        mat_filenames = [
            os.path.join(mat_dirname, filename) for filename in os.listdir(mat_dirname)
        ]
        mat_filenames.sort()
    return mat_filenames


def run_ants_apply_warp(  # noqa: PLR0913
    moving_image,
    reference,
    initial=None,
    rigid=None,
    affine=None,
    nonlinear=None,
    func_to_anat=None,
    anatomical_brain=None,
    dim=3,
    interp="Linear",
    inverse=False,
):
    """Apply a transform using ANTs transforms."""
    import os
    import subprocess

    if func_to_anat:
        # this assumes the func->anat affine transform is FSL-based and needs
        # to be converted to ITK format via c3d_affine_tool
        cmd = [
            "c3d_affine_tool",
            "-ref",
            anatomical_brain,
            "-src",
            moving_image,
            func_to_anat,
            "-fsl2ras",
            "-oitk",
            "affine.txt",
        ]
        subprocess.check_output(cmd)
        func_to_anat = change_itk_transform_type(
            os.path.join(os.getcwd(), "affine.txt")
        )

    out_image = os.path.join(
        os.getcwd(),
        moving_image[moving_image.rindex("/") + 1 : moving_image.rindex(".nii.gz")]
        + "_warp.nii.gz",
    )

    cmd = [
        "antsApplyTransforms",
        "-d",
        str(dim),
        "-i",
        moving_image,
        "-r",
        reference,
        "-o",
        out_image,
        "-n",
        interp,
    ]

    if nonlinear:
        cmd.append("-t")
        if inverse:
            cmd.append(f"[{os.path.abspath(nonlinear)}, 1]")
        else:
            cmd.append(os.path.abspath(nonlinear))

    if affine:
        cmd.append("-t")
        if inverse:
            cmd.append(f"[{os.path.abspath(affine)}, 1]")
        else:
            cmd.append(os.path.abspath(affine))

    if rigid:
        cmd.append("-t")
        if inverse:
            cmd.append(f"[{os.path.abspath(rigid)}, 1]")
        else:
            cmd.append(os.path.abspath(rigid))

    if initial:
        cmd.append("-t")
        if inverse:
            cmd.append(f"[{os.path.abspath(initial)}, 1]")
        else:
            cmd.append(os.path.abspath(initial))

    if func_to_anat:
        cmd.append("-t")
        if inverse:
            cmd.append(f"[{os.path.abspath(func_to_anat)}, 1]")
        else:
            cmd.append(os.path.abspath(func_to_anat))

    subprocess.check_output(cmd)

    return out_image


def run_c3d(reference_file, source_file, transform_file):
    """Run c3d_affine_tool to convert an FSL FLIRT affine transform to ITK."""
    import os
    import subprocess

    itk_transform = os.path.join(os.getcwd(), "affine.txt")

    cmd = [
        "c3d_affine_tool",
        "-ref",
        reference_file,
        "-src",
        source_file,
        transform_file,
        "-fsl2ras",
        "-oitk",
        itk_transform,
    ]
    subprocess.check_output(cmd)

    return itk_transform


def run_c4d(input_name, output_name):
    """Run c4d to split a 4D image into 3D images."""
    import os

    output1 = os.path.join(os.getcwd(), output_name + "1.nii.gz")
    output2 = os.path.join(os.getcwd(), output_name + "2.nii.gz")
    output3 = os.path.join(os.getcwd(), output_name + "3.nii.gz")

    cmd = f"c4d -mcs {input_name} -oo {output1} {output2} {output3}"
    os.system(cmd)

    return output1, output2, output3


@overload
def prepend_space(resource: list[str], space: str) -> list[str]: ...
@overload
def prepend_space(resource: str, space: str) -> str: ...
def prepend_space(resource: str | list[str], space: str) -> str | list[str]:
    """Given a resource or list of resources, return same but with updated space."""
    if isinstance(resource, list):
        return [prepend_space(_, space) for _ in resource]
    prefix = "longitudinal-template_" if space == "longitudinal" else ""
    if "space" not in resource:
        return f"{prefix}space-{space}_{resource}"
    pre, post = resource.split("space-")
    _old_space, post = post.split("_", 1)
    return f"{prefix}space-{space}_".join([pre, post])


def collect_xfms(
    wf: Workflow,
    name: str,
    node_inputs: list[tuple[Node | Workflow, list[str]]],
    *,
    mem_gb: float = 0.8,
    mem_x: tuple[float, str] = (263474863123069 / 37778931862957161709568, "in1"),
) -> Node:
    """Collect transforms to compose into one."""
    if len(node_inputs) == 1:
        input_node, inputs = node_inputs
    else:
        msg = "Combining transforms from multiple nodes not yet implemented."
        raise NotImplementedError(msg)
    numinputs = len(inputs)
    node = Node(interface=Merge(numinputs), name=name, mem_gb=mem_gb, mem_x=mem_x)
    # check transform list to exclude Nonetype (missing) init/rig/affine
    check_transform = Node(
        Function(
            input_names=["transform_list"],
            output_names=["checked_transform_list", "list_length"],
            function=check_transforms,
        ),
        name="check_transforms",
        mem_gb=6,
    )
    wf.connect(
        [
            (
                input_node,
                node,
                [(node_input, f"in{i + 1}") for i, node_input in enumerate(inputs)],
            ),
            (node, check_transform, [("out", "transform_list")]),
        ]
    )
    return node


def compose_ants_warp(  # noqa: PLR0913
    wf: Workflow,
    name: str,
    input_node: NipypeNode | Workflow,
    warp_from: str,
    warp_to: str,
    inputs: list[tuple[Node | Workflow, list[str]]],
    inv: bool = False,
    *,
    dimension: int = 3,
    input_image_type: int = 0,
    mem_gb: float = 1.155,
    mem_x: tuple[float, str] = (
        1708448960473801 / 1208925819614629174706176,
        "input_image",
    ),
) -> Node:
    """Create a Node to combine xfms."""
    node = Node(
        interface=AntsApplyTransforms(
            dimension=dimension,
            input_image_type=input_image_type,
            print_out_composite_warp_file=True,
        ),
        name=f"write_composite_{name}",
        mem_gb=mem_gb,
        mem_x=mem_x,
    )
    wf.connect(
        [
            (
                input_node,
                node,
                [
                    (warp_from, "input_image"),
                    (warp_to, "reference_image"),
                    ("interpolation", "interpolation"),
                ],
            )
        ]
    )
    collect_transforms = collect_xfms(wf, f"collect_{name}", inputs)
    wf.connect(collect_transforms, "checked_transform_list", node, "transforms")
    if inv:
        # generate inverse transform flags, which depends on the
        # number of transforms
        inverse_transform_flags = Node(
            Function(
                input_names=["transform_list"],
                output_names=["inverse_transform_flags"],
                function=generate_inverse_transform_flags,
            ),
            name="inverse_transform_flags",
        )
        wf.connect(
            collect_transforms,
            "checked_transform_list",
            inverse_transform_flags,
            "transform_list",
        )
        wf.connect(
            inverse_transform_flags,
            "inverse_transform_flags",
            node,
            "invert_transform_flags",
        )
    return node
