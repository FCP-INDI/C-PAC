# Copyright (C) 2019-2024  C-PAC Developers

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
"""General utilities for nuisance regression."""
from collections import OrderedDict
import os
import re
from typing import Optional

from nipype.interfaces import afni, ants, fsl
import nipype.interfaces.utility as util
from nipype.pipeline.engine import Workflow

from CPAC.nuisance.utils.crc import encode as crc_encode
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.registration.utils import generate_inverse_transform_flags
from CPAC.utils.interfaces.fsl import Merge as fslMerge
from CPAC.utils.interfaces.function import Function
from CPAC.utils.monitoring import IFLOGGER


def find_offending_time_points(
    fd_j_file_path=None,
    fd_p_file_path=None,
    dvars_file_path=None,
    fd_j_threshold=None,
    fd_p_threshold=None,
    dvars_threshold=None,
    number_of_previous_trs_to_censor=0,
    number_of_subsequent_trs_to_censor=0,
):
    """
    Find time points whose FD and/or DVARS are > threshold.

    :param fd_j_file_path: path to TSV containing framewise displacement as a
        single column. If not specified, it will not be used.
    :param fd_p_file_path: path to TSV containing framewise displacement as a
        single column. If not specified, it will not be used.
    :param dvars_file_path: path to TSV containing DVARS as a single column.
        If not specified, it will not be used.
    :param fd_j_threshold: threshold to apply to framewise displacement (Jenkinson),
        it can be a value such as 0.2 or a floating point multiple of the
        standard deviation specified as, e.g. '1.5SD'.
    :param fd_p_threshold: threshold to apply to framewise displacement (Power),
        it can be a value such as 0.2 or a floating point multiple of the
        standard deviation specified as, e.g. '1.5SD'.
    :param dvars_threshold: threshold to apply to DVARS, can be a value such
        as 0.5 or a floating point multiple of the standard deviation specified
        as, e.g. '1.5SD'.
    :param number_of_previous_trs_to_censor: extent of censorship window before
        the censor.
    :param number_of_subsequent_trs_to_censor: extent of censorship window after
        the censor.

    :return: File path to TSV file containing the volumes to be censored.
    """
    import os
    import re

    import numpy as np

    offending_time_points = set()
    time_course_len = 0

    motion_measures = ["FDJ", "FDP", "DVARS"]
    file_paths = [fd_j_file_path, fd_p_file_path, dvars_file_path]
    thresholds = [fd_j_threshold, fd_p_threshold, dvars_threshold]

    for motion_measure, file_path, _threshold in zip(
        motion_measures, file_paths, thresholds
    ):
        threshold = _threshold
        if not file_path:
            continue

        if not os.path.isfile(file_path):
            msg = f"File {file_path} could not be found."
            raise ValueError(msg)

        if not threshold:
            msg = "Method requires the specification of a threshold, none received"
            raise ValueError(msg)

        metric = np.loadtxt(file_path)
        if motion_measure == "DVARS":
            metric = np.array([0.0, *metric.tolist()])

        if not time_course_len:
            time_course_len = metric.shape[0]
        else:
            assert (
                time_course_len == metric.shape[0]
            ), "Threshold metric files does not have same size."

        try:
            threshold_sd = re.match(r"([0-9]*\.*[0-9]*)\s*SD", str(threshold))

            if threshold_sd:
                threshold_sd = float(threshold_sd.groups()[0])
                threshold = metric.mean() + threshold_sd * metric.std()
            else:
                threshold = float(threshold)
        except (AttributeError, re.error, IndexError, TypeError, ValueError):
            msg = f"Could not translate threshold {threshold} into a meaningful value"
            raise ValueError(msg)

        offending_time_points |= set(np.where(metric > threshold)[0].tolist())

    extended_censors = []
    for censor in offending_time_points:
        extended_censors += list(
            range(
                (censor - number_of_previous_trs_to_censor),
                (censor + number_of_subsequent_trs_to_censor + 1),
            )
        )

    extended_censors = [
        censor
        for censor in np.unique(extended_censors)
        if 0 <= censor < time_course_len
    ]

    censor_vector = np.ones((time_course_len, 1))
    censor_vector[extended_censors] = 0

    out_file_path = os.path.join(os.getcwd(), "censors.tsv")
    np.savetxt(out_file_path, censor_vector, fmt="%d", header="censor", comments="")

    return out_file_path


def compute_threshold(in_file, mask, threshold):
    """Return a given threshold."""
    return threshold


def compute_pct_threshold(in_file, mask, threshold_pct):
    """Compute the threshold based on the percentile of the data."""
    import numpy as np
    import nibabel as nib

    m = nib.load(mask).get_fdata().astype(bool)
    if not np.any(m):
        return 0.0
    d = nib.load(in_file).get_fdata()[m]
    return np.percentile(d, 100.0 - threshold_pct)


def compute_sd_threshold(in_file, mask, threshold_sd):
    """Compute the threshold based on the mean and standard deviation of the data."""
    import numpy as np
    import nibabel as nib

    m = nib.load(mask).get_fdata().astype(bool)
    if not np.any(m):
        return 0.0
    d = nib.load(in_file).get_fdata()[m]
    return d.mean() + threshold_sd * d.std()


def temporal_variance_mask(
    threshold, by_slice=False, erosion=False, degree=1
) -> Workflow:
    """Create a mask based on the temporal variance of the data."""
    threshold_method = "VAR"

    if isinstance(threshold, str):
        regex_match = {
            "SD": r"([0-9]+(\.[0-9]+)?)\s*SD",
            "PCT": r"([0-9]+(\.[0-9]+)?)\s*PCT",
        }

        for method, regex in regex_match.items():
            matched = re.match(regex, threshold)
            if matched:
                threshold_method = method
                threshold_value = matched.groups()[0]

    try:
        threshold_value = float(threshold_value)
    except (TypeError, ValueError):
        msg = (
            f"Error converting threshold value {threshold_value} from {threshold} to a "
            "floating point number. The threshold value can "
            "contain SD or PCT for selecting a threshold based on "
            "the variance distribution, otherwise it should be a "
            "floating point number."
        )
        raise ValueError(msg)

    if threshold_value < 0:
        msg = f"Threshold value should be positive, instead of {threshold_value}."
        raise ValueError(msg)

    if threshold_method == "PCT" and threshold_value >= 100.0:  # noqa: PLR2004
        msg = f"Percentile should be less than 100, received {threshold_value}."
        raise ValueError(msg)

    threshold = threshold_value

    wf = pe.Workflow(name="tcompcor")

    input_node = pe.Node(
        util.IdentityInterface(fields=["functional_file_path", "mask_file_path"]),
        name="inputspec",
    )
    output_node = pe.Node(util.IdentityInterface(fields=["mask"]), name="outputspec")

    # C-PAC default performs linear regression while nipype performs quadratic regression
    detrend = pe.Node(
        afni.Detrend(args=f"-polort {degree}", outputtype="NIFTI"),
        name="detrend",
    )
    wf.connect(input_node, "functional_file_path", detrend, "in_file")

    std = pe.Node(afni.TStat(args="-nzstdev", outputtype="NIFTI"), name="std")
    wf.connect(input_node, "mask_file_path", std, "mask")
    wf.connect(detrend, "out_file", std, "in_file")

    var = pe.Node(afni.Calc(expr="a*a", outputtype="NIFTI"), name="var")
    wf.connect(std, "out_file", var, "in_file_a")

    if by_slice:
        slices = pe.Node(fsl.Slice(), name="slicer")
        wf.connect(var, "out_file", slices, "in_file")

        mask_slices = pe.Node(fsl.Slice(), name="mask_slicer")
        wf.connect(input_node, "mask_file_path", mask_slices, "in_file")

        mapper = pe.MapNode(
            util.IdentityInterface(fields=["out_file", "mask_file"]),
            name="slice_mapper",
            iterfield=["out_file", "mask_file"],
        )
        wf.connect(slices, "out_files", mapper, "out_file")
        wf.connect(mask_slices, "out_files", mapper, "mask_file")

    else:
        mapper_list = pe.Node(util.Merge(1), name="slice_mapper_list")
        wf.connect(var, "out_file", mapper_list, "in1")

        mask_mapper_list = pe.Node(util.Merge(1), name="slice_mask_mapper_list")
        wf.connect(input_node, "mask_file_path", mask_mapper_list, "in1")

        mapper = pe.Node(
            util.IdentityInterface(fields=["out_file", "mask_file"]),
            name="slice_mapper",
        )
        wf.connect(mapper_list, "out", mapper, "out_file")
        wf.connect(mask_mapper_list, "out", mapper, "mask_file")

    if threshold_method == "PCT":
        threshold_node = pe.MapNode(
            Function(
                input_names=["in_file", "mask", "threshold_pct"],
                output_names=["threshold"],
                function=compute_pct_threshold,
                as_module=True,
            ),
            name="threshold_value",
            iterfield=["in_file", "mask"],
        )
        threshold_node.inputs.threshold_pct = threshold_value
        wf.connect(mapper, "out_file", threshold_node, "in_file")
        wf.connect(mapper, "mask_file", threshold_node, "mask")

    elif threshold_method == "SD":
        threshold_node = pe.MapNode(
            Function(
                input_names=["in_file", "mask", "threshold_sd"],
                output_names=["threshold"],
                function=compute_sd_threshold,
                as_module=True,
            ),
            name="threshold_value",
            iterfield=["in_file", "mask"],
        )
        threshold_node.inputs.threshold_sd = threshold_value
        wf.connect(mapper, "out_file", threshold_node, "in_file")
        wf.connect(mapper, "mask_file", threshold_node, "mask")

    else:
        threshold_node = pe.MapNode(
            Function(
                input_names=["in_file", "mask", "threshold"],
                output_names=["threshold"],
                function=compute_threshold,
                as_module=True,
            ),
            name="threshold_value",
            iterfield=["in_file", "mask"],
        )
        threshold_node.inputs.threshold = threshold_value
        wf.connect(mapper, "out_file", threshold_node, "in_file")
        wf.connect(mapper, "mask_file", threshold_node, "mask")

    threshold_mask = pe.MapNode(
        interface=fsl.maths.Threshold(),
        name="threshold",
        iterfield=["in_file", "thresh"],
    )
    threshold_mask.inputs.args = "-bin"
    wf.connect(mapper, "out_file", threshold_mask, "in_file")
    wf.connect(threshold_node, "threshold", threshold_mask, "thresh")

    merge_slice_masks = pe.Node(interface=fslMerge(), name="merge_slice_masks")
    merge_slice_masks.inputs.dimension = "z"
    wf.connect(threshold_mask, "out_file", merge_slice_masks, "in_files")

    wf.connect(merge_slice_masks, "merged_file", output_node, "mask")

    return wf


def generate_summarize_tissue_mask(
    nuisance_wf,
    pipeline_resource_pool,
    regressor_descriptor,
    regressor_selector,
    csf_mask_exist,
    use_ants=True,
    ventricle_mask_exist=True,
    all_bold=False,
):
    """
    Add tissue mask generation into pipeline according to the selector.

    :param nuisance_wf: Nuisance regressor workflow.
    :param pipeline_resource_pool: dictionary of available resources.
    :param regressor_descriptor: dictionary of steps to build, including keys:
        'tissue', 'resolution', 'erosion'
    :param regressor_selector: dictionary with the original selector

    :return: the full path of the 3D nifti file containing the mask created by
        this operation.
    """
    steps = [
        key
        for key in ["tissue", "resolution", "erosion"]
        if key in regressor_descriptor
    ]

    full_mask_key = "_".join(regressor_descriptor[s] for s in steps)

    for step_i, step in enumerate(steps):
        mask_key = "_".join(regressor_descriptor[s] for s in steps[: step_i + 1])

        if mask_key in pipeline_resource_pool:
            continue

        node_mask_key = re.sub(r"[^\w]", "_", mask_key)

        prev_mask_key = "_".join(regressor_descriptor[s] for s in steps[:step_i])

        if step == "tissue":
            pass

        elif step == "resolution":
            if all_bold:
                pass

            if csf_mask_exist:
                mask_to_epi = pe.Node(
                    interface=fsl.FLIRT(),
                    name=f"{node_mask_key}_flirt",
                    mem_gb=3.63,
                    mem_x=(3767129957844731 / 1208925819614629174706176, "in_file"),
                )

                mask_to_epi.inputs.interp = "nearestneighbour"

                if regressor_selector["extraction_resolution"] == "Functional":
                    # apply anat2func matrix
                    mask_to_epi.inputs.apply_xfm = True
                    mask_to_epi.inputs.output_type = "NIFTI_GZ"
                    nuisance_wf.connect(
                        *(
                            pipeline_resource_pool["Functional_mean"]
                            + (mask_to_epi, "reference")
                        )
                    )
                    nuisance_wf.connect(
                        *(
                            pipeline_resource_pool["Transformations"][
                                "anat_to_func_linear_xfm"
                            ]
                            + (mask_to_epi, "in_matrix_file")
                        )
                    )

                else:
                    resolution = regressor_selector["extraction_resolution"]
                    mask_to_epi.inputs.apply_isoxfm = resolution

                    nuisance_wf.connect(
                        *(
                            pipeline_resource_pool[f"Anatomical_{resolution}mm"]
                            + (mask_to_epi, "reference")
                        )
                    )

                nuisance_wf.connect(
                    *(pipeline_resource_pool[prev_mask_key] + (mask_to_epi, "in_file"))
                )

                pipeline_resource_pool[mask_key] = (mask_to_epi, "out_file")

            if full_mask_key.startswith("CerebrospinalFluid"):
                pipeline_resource_pool = (
                    generate_summarize_tissue_mask_ventricles_masking(
                        nuisance_wf,
                        pipeline_resource_pool,
                        regressor_descriptor,
                        regressor_selector,
                        node_mask_key,
                        csf_mask_exist,
                        use_ants,
                        ventricle_mask_exist,
                    )
                )

        elif step == "erosion":
            erode_mask_node = pe.Node(
                afni.Calc(
                    args="-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k",
                    expr="a*(1-amongst(0,b,c,d,e,f,g))",
                    outputtype="NIFTI_GZ",
                ),
                name=f"{node_mask_key}",
            )

            nuisance_wf.connect(
                *(
                    pipeline_resource_pool[prev_mask_key]
                    + (erode_mask_node, "in_file_a")
                )
            )

            pipeline_resource_pool[mask_key] = (erode_mask_node, "out_file")

    return pipeline_resource_pool, full_mask_key


def generate_summarize_tissue_mask_ventricles_masking(
    nuisance_wf,
    pipeline_resource_pool: dict,
    regressor_descriptor,
    regressor_selector,
    mask_key,
    csf_mask_exist,
    use_ants=True,
    ventricle_mask_exist=True,
) -> Optional[dict]:
    """Update CSF mask to include only the lateral ventricles."""
    if not csf_mask_exist:
        IFLOGGER.warning(
            "Segmentation is Off, - therefore will be using "
            "lateral_ventricle_mask as CerebrospinalFluid_mask."
        )

    # Mask CSF with Ventricles
    if f"{mask_key}_Unmasked" not in pipeline_resource_pool:
        if ventricle_mask_exist:
            ventricles_key = "VentriclesToAnat"

            if "resolution" in regressor_descriptor:
                ventricles_key += "_{}".format(regressor_descriptor["resolution"])

            if ventricles_key not in pipeline_resource_pool:
                transforms = pipeline_resource_pool["Transformations"]

                if use_ants:
                    # perform the transform using ANTS
                    collect_linear_transforms = pe.Node(
                        util.Merge(3), name=f"{ventricles_key}_ants_transforms"
                    )

                    nuisance_wf.connect(
                        *(
                            transforms["mni_to_anat_linear_xfm"]
                            + (collect_linear_transforms, "in1")
                        )
                    )

                    # generate inverse transform flags, which depends on the number of transforms
                    inverse_transform_flags = pe.Node(
                        util.Function(
                            input_names=["transform_list"],
                            output_names=["inverse_transform_flags"],
                            function=generate_inverse_transform_flags,
                        ),
                        name=f"{ventricles_key}_inverse_transform_flags",
                    )
                    nuisance_wf.connect(
                        collect_linear_transforms,
                        "out",
                        inverse_transform_flags,
                        "transform_list",
                    )

                    lat_ven_mni_to_anat = pe.Node(
                        interface=ants.ApplyTransforms(),
                        name=f"{ventricles_key}_ants",
                        mem_gb=0.683,
                        mem_x=(
                            3811976743057169 / 302231454903657293676544,
                            "input_image",
                        ),
                    )
                    lat_ven_mni_to_anat.inputs.interpolation = "NearestNeighbor"
                    lat_ven_mni_to_anat.inputs.dimension = 3

                    nuisance_wf.connect(
                        inverse_transform_flags,
                        "inverse_transform_flags",
                        lat_ven_mni_to_anat,
                        "invert_transform_flags",
                    )
                    nuisance_wf.connect(
                        collect_linear_transforms,
                        "out",
                        lat_ven_mni_to_anat,
                        "transforms",
                    )

                    nuisance_wf.connect(
                        *(
                            pipeline_resource_pool["Ventricles"]
                            + (lat_ven_mni_to_anat, "input_image")
                        )
                    )
                    resolution = regressor_selector["extraction_resolution"]

                    if csf_mask_exist:
                        nuisance_wf.connect(
                            *(
                                pipeline_resource_pool[mask_key]
                                + (lat_ven_mni_to_anat, "reference_image")
                            )
                        )
                    elif resolution == "Functional":
                        nuisance_wf.connect(
                            *(
                                pipeline_resource_pool["Functional_mean"]
                                + (lat_ven_mni_to_anat, "reference_image")
                            )
                        )
                    else:
                        nuisance_wf.connect(
                            *(
                                pipeline_resource_pool[f"Anatomical_{resolution}mm"]
                                + (lat_ven_mni_to_anat, "reference_image")
                            )
                        )

                    pipeline_resource_pool[ventricles_key] = (
                        lat_ven_mni_to_anat,
                        "output_image",
                    )

                else:
                    # perform the transform using FLIRT
                    lat_ven_mni_to_anat = pe.Node(
                        interface=fsl.ApplyWarp(),
                        name=f"{ventricles_key}_fsl_applywarp",
                    )
                    lat_ven_mni_to_anat.inputs.interp = "nn"

                    nuisance_wf.connect(
                        *(
                            transforms["mni_to_anat_linear_xfm"]
                            + (lat_ven_mni_to_anat, "field_file")
                        )
                    )
                    nuisance_wf.connect(
                        *(
                            pipeline_resource_pool["Ventricles"]
                            + (lat_ven_mni_to_anat, "in_file")
                        )
                    )
                    nuisance_wf.connect(
                        *(
                            pipeline_resource_pool[mask_key]
                            + (lat_ven_mni_to_anat, "ref_file")
                        )
                    )

                    pipeline_resource_pool[ventricles_key] = (
                        lat_ven_mni_to_anat,
                        "out_file",
                    )

            if csf_mask_exist:
                # reduce CSF mask to the lateral ventricles
                mask_csf_with_lat_ven = pe.Node(
                    interface=afni.Calc(outputtype="NIFTI_GZ"),
                    name=f"{mask_key}_Ventricles",
                )
                mask_csf_with_lat_ven.inputs.expr = "a*b"
                mask_csf_with_lat_ven.inputs.out_file = "csf_lat_ven_mask.nii.gz"

                nuisance_wf.connect(
                    *(
                        pipeline_resource_pool[ventricles_key]
                        + (mask_csf_with_lat_ven, "in_file_a")
                    )
                )
                nuisance_wf.connect(
                    *(
                        pipeline_resource_pool[mask_key]
                        + (mask_csf_with_lat_ven, "in_file_b")
                    )
                )

                pipeline_resource_pool[f"{mask_key}_Unmasked"] = pipeline_resource_pool[
                    mask_key
                ]
                pipeline_resource_pool[mask_key] = (mask_csf_with_lat_ven, "out_file")

            else:
                pipeline_resource_pool[mask_key] = pipeline_resource_pool[
                    ventricles_key
                ]

        return pipeline_resource_pool
    return None


class NuisanceRegressor(object):
    """A nuisance regressor."""

    def __init__(self, selector):
        """Initialize the nuisance regressor."""
        self.selector = selector

        if "Bandpass" in self.selector:
            s = self.selector["Bandpass"]
            if not isinstance(s, dict) or (
                not s.get("bottom_frequency") and not s.get("top_frequency")
            ):
                del self.selector["Bandpass"]

    def get(self, key, default=None):
        """Return the value of the key in the selector."""
        return self.selector.get(key, default)

    def __contains__(self, key):
        """Return whether the key is in the selector."""
        return key in self.selector

    def __getitem__(self, key):
        """Return the value of the key in the selector."""
        return self.selector[key]

    @staticmethod
    def _derivative_params(selector):
        nr_repr = ""
        if not selector:
            return nr_repr
        if selector.get("include_squared"):
            nr_repr += "S"
        if selector.get("include_delayed"):
            nr_repr += "D"
        if selector.get("include_delayed_squared"):
            nr_repr += "B"
        if selector.get("include_backdiff"):
            nr_repr += "V"
        if selector.get("include_backdiff_squared"):
            nr_repr += "C"
        return nr_repr

    @staticmethod
    def _summary_params(selector):
        summ = selector["summary"]

        methods = {
            "PC": "PC",
            "DetrendPC": "DPC",
            "Mean": "M",
            "NormMean": "NM",
            "DetrendMean": "DM",
            "DetrendNormMean": "DNM",
        }

        if isinstance(summ, dict):
            method = summ["method"]
            rep = methods[method]
            if method in ["DetrendPC", "PC"]:
                rep += "%d" % summ["components"]
        else:
            rep = methods[summ]

        return rep

    @staticmethod
    def encode(selector: dict) -> str:
        """Return a brief string representation of the nuisance regressor."""
        regs = OrderedDict(
            [
                ("GreyMatter", "GM"),
                ("WhiteMatter", "WM"),
                ("CerebrospinalFluid", "CSF"),
                ("tCompCor", "tC"),
                ("aCompCor", "aC"),
                ("GlobalSignal", "G"),
                ("Motion", "M"),
                ("Custom", "T"),
                ("PolyOrt", "P"),
                ("Bandpass", "BP"),
                ("Censor", "C"),
            ]
        )

        tissues = ["GreyMatter", "WhiteMatter", "CerebrospinalFluid"]

        selectors_representations = []

        # tC-1.5PT-PC5S-SDB
        # aC-WC-2mmE-PC5-SDB

        # WM-2mmE-PC5-SDB
        # CSF-2mmE-M-SDB
        # GM-2mmE-DNM-SDB

        # G-PC5-SDB
        # M-SDB
        # C-S-FD1.5SD-D1.5SD
        # P-2
        # B-T0.01-B0.1

        for r in regs.keys():
            if r not in selector:
                continue

            s = selector[r]

            pieces = [regs[r]]

            if r in tissues:
                if (
                    s.get("extraction_resolution")
                    and s["extraction_resolution"] != "Functional"
                ):
                    res = "%.2gmm" % s["extraction_resolution"]
                    if s.get("erode_mask"):
                        res += "E"
                    pieces += [res]

                pieces += [NuisanceRegressor._summary_params(s)]
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == "tCompCor":
                threshold = ""
                if s.get("by_slice"):
                    threshold += "S"
                t = s.get("threshold")
                if t:
                    if not isinstance(t, str):
                        t = "%.2f" % t
                    threshold += t
                if s.get("erode_mask"):
                    threshold += "E"
                if s.get("degree"):
                    d = s.get("degree")
                    threshold += str(d)

                pieces += [threshold]
                pieces += [NuisanceRegressor._summary_params(s)]
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == "aCompCor":
                if s.get("tissues"):
                    pieces += ["+".join([regs[t] for t in sorted(s["tissues"])])]

                if s.get("extraction_resolution"):
                    res = "%.2gmm" % s["extraction_resolution"]
                    if s.get("erode_mask"):
                        res += "E"
                    pieces += [res]

                pieces += [NuisanceRegressor._summary_params(s)]
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == "Custom":
                for ss in s:
                    pieces += [
                        os.path.basename(ss["file"])[0:5] + crc_encode(ss["file"])
                    ]

            elif r == "GlobalSignal":
                pieces += [NuisanceRegressor._summary_params(s)]
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == "Motion":
                pieces += [NuisanceRegressor._derivative_params(s)]

            elif r == "PolyOrt":
                pieces += ["%d" % s["degree"]]

            elif r == "Bandpass":
                if s.get("bottom_frequency"):
                    pieces += ["B%.2g" % s["bottom_frequency"]]
                if s.get("top_frequency"):
                    pieces += ["T%.2g" % s["top_frequency"]]

            elif r == "Censor":
                censoring = {
                    "Kill": "K",
                    "Zero": "Z",
                    "Interpolate": "I",
                    "SpikeRegression": "S",
                }

                thresholds = {
                    "FD_J": "FD-J",
                    "FD_P": "FD-P",
                    "DVARS": "DV",
                }

                pieces += [censoring[s["method"]]]

                trs_range = ["0", "0"]
                if s.get("number_of_previous_trs_to_censor"):
                    trs_range[0] = "%d" % s["number_of_previous_trs_to_censor"]
                if s.get("number_of_subsequent_trs_to_censor"):
                    trs_range[1] = "%d" % s["number_of_subsequent_trs_to_censor"]

                pieces += ["+".join(trs_range)]

                threshs = sorted(s["thresholds"], reverse=True, key=lambda d: d["type"])
                for st in threshs:
                    thresh = thresholds[st["type"]]
                    if isinstance(st["value"], str):
                        thresh += st["value"]
                    else:
                        thresh += "%.2g" % st["value"]

                    pieces += [thresh]

            selectors_representations += ["-".join([_f for _f in pieces if _f])]

        return "_".join(selectors_representations)

    def __repr__(self) -> str:
        """Return a string representation of the nuisance regressor."""
        return NuisanceRegressor.encode(self.selector)
