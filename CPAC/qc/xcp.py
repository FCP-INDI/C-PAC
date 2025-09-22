# Copyright (C) 2021-2025  C-PAC Developers

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
"""
Generate XCP-stype quality control files.

.. seealso::

      `User Guide: Generate eXtensible Connectivity Pipeline-style quality control files <https://fcp-indi.github.io/docs/user/xcpqc>`_

Columns
-------
sub : str
    subject label :footcite:`BIDS21`
ses : str
    session label :footcite:`BIDS21`
task : str
    task label :footcite:`BIDS21`
run : int
    run index :footcite:`BIDS21`
desc : str
    description :footcite:`BIDS21`
regressors : str
    'Name' of regressors in the current fork
space : str
    space label :footcite:`BIDS21`
meanFD : float
    mean Jenkinson framewise displacement :footcite:`xcp_22,Jenk02` :func:`CPAC.generate_motion_statistics.calculate_FD_J` after preprocessing
relMeansRMSMotion : float
    "mean value of RMS motion" :footcite:`xcp_22,Ciri19`
relMaxRMSMotion : float
    "maximum vaue of RMS motion" :footcite:`xcp_22,Ciri19`
meanDVInit : float
    "mean DVARS" :footcite:`xcp_22,Ciri19`
meanDVFinal : float
    "mean DVARS" :footcite:`xcp_22,Ciri19`
nVolCensored : int
    "total number of volume(s) censored :footcite:`Ciri19`
nVolsRemoved : int
    number of volumes in derivative minus number of volumes in original
    functional scan
motionDVCorrInit : float
    "correlation of RMS and DVARS before regresion" :footcite:`Ciri19`
motionDVCorrFinal : float
    "correlation of RMS and DVARS after regresion" :footcite:`Ciri19`
coregDice : float
    "Coregsitration of Functional and T1w:[…] Dice index" :footcite:`xcp_22,Ciri19`
coregJaccard : float
    "Coregsitration of Functional and T1w:[…] Jaccard index" :footcite:`xcp_22,Ciri19`
coregCrossCorr : float
    "Coregsitration of Functional and T1w:[…] cross correlation" :footcite:`xcp_22,Ciri19`
coregCoverag : float
    "Coregsitration of Functional and T1w:[…] Coverage index" :footcite:`xcp_22,Ciri19`
normDice : float
    "Normalization of T1w/Functional to Template:[…] Dice index" :footcite:`xcp_22,Ciri19`
normJaccard : float
    "Normalization of T1w/Functional to Template:[…] Jaccard index" :footcite:`xcp_22,Ciri19`
normCrossCorr : float
    "Normalization of T1w/Functional to Template:[…] cross correlation" :footcite:`xcp_22,Ciri19`
normCoverage : float
    "Normalization of T1w/Functional to Template:[…] Coverage index" :footcite:`xcp_22,Ciri19`
"""  # pylint: disable=line-too-long

from io import BufferedReader
import os
import re
from typing import Any, Optional

from bids.layout import parse_file_entities
import numpy as np
import pandas as pd
from nibabel import load as nib_load  # type: ignore[reportPrivateImportUsage]
from nipype.interfaces import afni, fsl

from CPAC.generate_motion_statistics.generate_motion_statistics import (
    DVARS_strip_t0,
    ImageTo1D,
)
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.engine import ResourcePool
from CPAC.pipeline.nodeblock import nodeblock
from CPAC.qc.qcmetrics import regisQ
from CPAC.utils.interfaces.function import Function

motion_params = [
    "dvars",
    "framewise-displacement-jenkinson",
    "desc-movementParametersUnfiltered_motion",
    "desc-movementParameters_motion",
]


def _connect_motion(
    wf: pe.Workflow, strat_pool: ResourcePool, qc_file: pe.Node, pipe_num: int
) -> pe.Workflow:
    """
    Connect the motion metrics to the workflow.

    Parameters
    ----------
    wf
        The workflow to connect the motion metrics to.

    strat_pool
        The current strategy pool.

    qc_file
        A function node with the function ``generate_xcp_qc``.
    """
    # pylint: disable=invalid-name, too-many-arguments
    if strat_pool.check_rpool("censor-indices"):
        wf.connect_optional(strat_pool, "censor-indices", qc_file, "censor_indices")
    else:
        qc_file.inputs.censor_indices = []
    cal_DVARS = pe.Node(
        ImageTo1D(method="dvars"),
        name=f"cal_DVARS_{pipe_num}",
        mem_gb=0.4,
        mem_x=(739971956005215 / 151115727451828646838272, "in_file"),
        throttle=True,
    )
    cal_DVARS_strip = pe.Node(
        Function(
            input_names=["file_1D"],
            output_names=["out_file", "out_matrix"],
            function=DVARS_strip_t0,
            as_module=True,
        ),
        name=f"cal_DVARS_strip_{pipe_num}",
    )
    motion_name = "desc-movementParametersUnfiltered_motion"
    if not strat_pool.check_rpool(motion_name):
        motion_name = "desc-movementParameters_motion"
    wf.connect_optional(strat_pool, "desc-preproc_bold", cal_DVARS, "in_file")
    wf.connect_optional(strat_pool, "space-bold_desc-brain_mask", cal_DVARS, "mask")
    wf.connect_optional(strat_pool, motion_name, qc_file, "movement_parameters")
    for resource in motion_params:
        if not resource.endswith("_motion"):
            wf.connect_optional(
                strat_pool, resource, qc_file, resource.replace("-", "_")
            )
    wf.connect(
        [
            (cal_DVARS, cal_DVARS_strip, [("out_file", "file_1D")]),
            (
                cal_DVARS_strip,
                qc_file,
                [("out_file", "dvars_after_path"), ("out_matrix", "dvars_after")],
            ),
        ]
    )
    return wf


def dvcorr(dvars, fdj):
    """Correlate DVARS and FD-J."""
    dvars = np.loadtxt(dvars)
    fdj = np.loadtxt(fdj)
    if len(dvars) != len(fdj) - 1:
        msg = (
            "len(DVARS) should be 1 less than len(FDJ), but their respective "
            f"lengths are {len(dvars)} and {len(fdj)}."
        )
        raise ValueError(msg)
    return np.corrcoef(dvars, fdj[1:])[0, 1]


# This function is for a function node for which
# Nipype will connect many other nodes as inputs
def generate_xcp_qc(  # noqa: PLR0913
    sub: str,
    ses: str,
    task: str,
    run: str | int,
    desc: str,
    regressors: Optional[str],
    bold2t1w_mask: Optional[str],
    t1w_mask: Optional[str],
    bold2template_mask: Optional[str],
    template_mask: Optional[str],
    original_func: Optional[str],
    final_func: Optional[str],
    movement_parameters: Optional[str],
    dvars: Optional[str],
    censor_indices: Optional[list[int]],
    framewise_displacement_jenkinson: Optional[str],
    dvars_after: Optional[np.ndarray],
    dvars_after_path: Optional[str],
    template: Optional[str],
) -> str:
    """
    Generate an RBC-style QC CSV.

    Parameters
    ----------
    sub
        subject ID

    ses
        session ID

    task
        task ID

    run
        run ID

    desc
        description string

    regressors
        'Name' of regressors in fork

    original_func
        path to original 'bold' image

    final_bold
        path to 'space-template_desc-preproc_bold' image

    bold2t1w_mask
        path to bold-to-T1w transform applied to space-bold_desc-brain_mask
        with space-T1w_desc-brain_mask reference

    t1w_mask
        path to space-T1w_desc-brain_mask

    bold2template_mask
        path to space-template_desc-bold_mask

    template_mask
        path to T1w-brain-template-mask or EPI-template-mask

    movement_parameters
        path to movement parameters

    dvars
        path to DVARS before motion correction

    censor_indices
        list of indices of censored volumes

    framewise_displacement_jenkinson
        path to framewise displacement (Jenkinson) before motion correction

    dvars_after
        DVARS matrix for final 'bold' image

    dvars_after_path
        path to DVARS matrix for final 'bold' image

    template
        path to registration template

    Returns
    -------
    str
        path to space-template_desc-xcp_quality TSV
    """
    images = {
        key: nib_load(image)
        for key, image in [("original_func", original_func), ("final_func", final_func)]
        if image
    }

    qc_dict: dict[str, Any] = {}
    """Quality control dictionary to be converted to a DataFrame."""

    columns: list[str] = []
    """Header for the quality control DataFrame."""

    # `sub` through `space`
    from_bids: dict[str, Any] = {
        _k: _v
        for _k, _v in {
            "sub": sub,
            "ses": ses,
            "task": task,
            "run": run,
            "desc": desc,
            "regressors": regressors,
            "space": os.path.basename(template).split(".", 1)[0].split("_", 1)[0]
            if template
            else None,
        }.items()
        if _v is not None
    }
    columns.extend(["sub", "ses", "task", "run", "desc", "regressors"])
    if from_bids["space"] is not None:
        if from_bids["space"].startswith("tpl-"):
            from_bids["space"] = from_bids["space"][4:]
        columns.append("space")
    qc_dict = {**qc_dict, **from_bids}

    if framewise_displacement_jenkinson is not None:
        # `meanFD (Jenkinson)`
        power_params = {"meanFD": np.mean(np.loadtxt(framewise_displacement_jenkinson))}
        qc_dict = {**qc_dict, **power_params}
        columns.append("meanFD")

    if movement_parameters is not None:
        # `relMeansRMSMotion` & `relMaxRMSMotion`
        mot = np.genfromtxt(movement_parameters).T
        # Relative RMS of translation
        rms = np.sqrt(mot[3] ** 2 + mot[4] ** 2 + mot[5] ** 2)
        rms_params = {
            "relMeansRMSMotion": [np.mean(rms)],
            "relMaxRMSMotion": [np.max(rms)],
        }
        qc_dict = {**qc_dict, **rms_params}
        columns.extend(list(rms_params.keys()))

    if dvars is not None:
        # `meanDVInit` & `meanDVFinal`
        meanDV = {"meanDVInit": np.mean(np.loadtxt(dvars))}
        try:
            meanDV["motionDVCorrInit"] = dvcorr(dvars, framewise_displacement_jenkinson)
        except ValueError as value_error:
            meanDV["motionDVCorrInit"] = f"ValueError({value_error!s})"
        meanDV["meanDVFinal"] = np.mean(dvars_after)
        try:
            meanDV["motionDVCorrFinal"] = dvcorr(
                dvars_after_path, framewise_displacement_jenkinson
            )
        except ValueError as value_error:
            meanDV["motionDVCorrFinal"] = f"ValueError({value_error!s})"
        qc_dict = {**qc_dict, **meanDV}
        columns.extend(list(meanDV.keys()))

    if censor_indices is not None:
        # `nVolCensored` & `nVolsRemoved`
        n_vols_censored = (
            len(censor_indices) if censor_indices is not None else "unknown"
        )
        shape_params = {
            "nVolCensored": n_vols_censored,
            "nVolsRemoved": images["original_func"].shape[3]
            - images["final_func"].shape[3],
        }
        qc_dict = {**qc_dict, **shape_params}
        if "motionDVCorrFinal" in columns:
            columns.insert(-2, "meanDVInit")
            columns.insert(-2, "meanDVFinal")
            columns.extend(["motionDVCorrInit", "motionDVCorrFinal"])
        else:
            columns.extend(list(shape_params.keys()))

    if final_func is not None:
        if isinstance(final_func, BufferedReader):
            final_func = final_func.name

        desc_span = re.search(r"_desc-.*_", final_func) if final_func else None
        if desc_span:
            desc_span = desc_span.span()
            assert final_func is not None
            final_func = "_".join(
                [final_func[: desc_span[0]], final_func[desc_span[1] :]]
            )
        del desc_span

    if all(
        _var is not None
        for _var in [bold2t1w_mask, t1w_mask, bold2template_mask, template_mask]
    ):
        # `coregDice`, `coregJaccard`, `coregCrossCorr`, `coregCoverage`
        coreg_params = regisQ(
            bold2t1w_mask=bold2t1w_mask,
            t1w_mask=t1w_mask,
            bold2template_mask=bold2template_mask,
            template_mask=template_mask,
        )
        # Overlap
        overlap_params = regisQ(
            bold2t1w_mask=bold2t1w_mask,
            t1w_mask=t1w_mask,
            bold2template_mask=bold2template_mask,
            template_mask=template_mask,
        )
        qc_dict = {**qc_dict, **coreg_params, **overlap_params}
        columns.extend([*list(coreg_params.keys()), *list(overlap_params.keys())])

    df = pd.DataFrame(qc_dict, columns=columns)

    qc_filepath = os.path.join(os.getcwd(), "xcpqc.tsv")
    df.to_csv(qc_filepath, sep="\t", index=False)
    return qc_filepath


def get_bids_info(subject, scan, wf_name):
    """
    Gather BIDS information from a strat_pool.

    Parameters
    ----------
    subject : str
        subject ID

    scan : str
        scan ID

    wf_name : str
        workflow name

    Returns
    -------
    subject : str
        subject ID

    session : str
        session ID

    task : str
        task ID

    run : str or int
        run ID

    Examples
    --------
    >>> subject, session, task, run = get_bids_info(
    ...     subject='DavidBowman', scan='rest_acq-1_run-1',
    ...     wf_name='cpac_DavidBowman_3')
    >>> subject
    'DavidBowman'
    >>> session
    '3'
    >>> task
    'rest'
    >>> run
    '1'
    >>> get_bids_info(subject='sub-colornest035', scan='rest_run-01',
    ...               wf_name='cpac_sub-colornest035_ses-1')
    ('colornest035', '1', 'rest', '01')
    """
    returns = ("subject", "session", "task", "run")
    ses = wf_name.split("_")[-1]
    if not ses.startswith("ses-"):
        ses = f"ses-{ses}"
    if not subject.startswith("sub-"):
        subject = f"sub-{subject}"
    resource = "_".join([subject, ses, scan if "task" in scan else f"task-{scan}"])
    entities = parse_file_entities(resource)
    returns = {key: entities.get(key) for key in returns}
    if any(value is None for value in returns.values()):
        entity_parts = resource.split("_")

        def get_entity_part(key):
            key = f"{key}-"
            matching_parts = [part for part in entity_parts if part.startswith(key)]
            if matching_parts:
                return matching_parts[0].replace(key, "")
            return None

        for key, value in returns.items():
            if value is None:
                if key == "task":
                    returns[key] = get_entity_part(key)
                else:
                    returns[key] = get_entity_part(key[:3])
    return tuple(str(returns.get(key)) for key in ["subject", "session", "task", "run"])


@nodeblock(
    name="qc_xcp",
    config=["pipeline_setup", "output_directory", "quality_control"],
    switch=["generate_xcpqc_files"],
    inputs=[
        (
            "subject",
            "scan",
            "bold",
            "desc-preproc_bold",
            "space-T1w_sbref",
            "space-T1w_desc-brain_mask",
            "max-displacement",
            "space-template_desc-preproc_bold",
            "space-bold_desc-brain_mask",
            ["T1w-brain-template-mask", "EPI-template-mask"],
            ["space-template_desc-bold_mask", "space-EPItemplate_desc-bold_mask"],
            "desc-confounds_timeseries",
            ["T1w-brain-template-funcreg", "EPI-brain-template-funcreg"],
            [
                "desc-movementParametersUnfiltered_motion",
                "desc-movementParameters_motion",
            ],
            "dvars",
            "framewise-displacement-jenkinson",
        )
    ],
    outputs={
        "space-template_desc-xcp_quality": {"Template": "T1w-brain-template-mask"}
    },
)
def qc_xcp(wf, cfg, strat_pool, pipe_num, opt=None):
    """Generate an XCP-style quality-control file."""
    # pylint: disable=invalid-name, unused-argument
    if cfg[
        "nuisance_corrections", "2-nuisance_regression", "run"
    ] and not strat_pool.check_rpool("desc-confounds_timeseries"):
        return wf, {}
    bids_info = pe.Node(
        Function(
            input_names=["subject", "scan", "wf_name"],
            output_names=["subject", "session", "task", "run"],
            imports=["from bids.layout import parse_file_entities"],
            function=get_bids_info,
            as_module=True,
        ),
        name=f"bids_info_{pipe_num}",
    )
    bids_info.inputs.wf_name = wf.name
    qc_file = pe.Node(
        Function(
            input_names=[
                "sub",
                "ses",
                "task",
                "run",
                "desc",
                "bold2t1w_mask",
                "t1w_mask",
                "bold2template_mask",
                "template_mask",
                "original_func",
                "final_func",
                "template",
                "movement_parameters",
                "dvars",
                "censor_indices",
                "regressors",
                "framewise_displacement_jenkinson",
                "dvars_after_path",
                "dvars_after",
            ],
            output_names=["qc_file"],
            function=generate_xcp_qc,
            as_module=True,
        ),
        name=f"qcxcp_{pipe_num}",
    )
    qc_file.inputs.desc = "preproc"
    qc_file.inputs.regressors = (
        strat_pool.node_data("desc-confounds_timeseries")
        .node.name.split("desc-confounds_timeseries_")[-1][::-1]
        .split("_", 1)[-1][::-1]
    )
    bold_to_T1w_mask = pe.Node(
        interface=fsl.ImageMaths(),
        name=f"binarize_bold_to_T1w_mask_{pipe_num}",
        op_string="-bin ",
    )
    resample_bold_mask_to_template = pe.Node(
        afni.Resample(),
        name=f"resample_bold_mask_to_anat_res_{pipe_num}",
        mem_gb=0,
        mem_x=(0.0115, "in_file", "t"),
    )
    resample_bold_mask_to_template.inputs.outputtype = "NIFTI_GZ"
    wf: pe.Workflow = _connect_motion(wf, strat_pool, qc_file, pipe_num=pipe_num)
    if not hasattr(wf, "connect_optional"):
        setattr(wf, "connect_optional", pe.Workflow.connect_optional)

    for key in ["subject", "scan"]:
        wf.connect_optional(strat_pool, key, bids_info, key)
    wf.connect_optional(strat_pool, "space-T1w_sbref", bold_to_T1w_mask, "in_file")
    wf.connect_optional(strat_pool, "space-T1w_desc-brain_mask", qc_file, "t1w_mask")
    wf.connect_optional(
        strat_pool,
        ["T1w-brain-template-mask", "EPI-template-mask"],
        qc_file,
        "template_mask",
    )
    wf.connect_optional(
        strat_pool,
        ["T1w-brain-template-mask", "EPI-template-mask"],
        resample_bold_mask_to_template,
        "master",
    )
    wf.connect_optional(strat_pool, "bold", qc_file, "original_func")
    wf.connect_optional(
        strat_pool, "space-template_desc-preproc_bold", qc_file, "final_func"
    )
    wf.connect_optional(
        strat_pool,
        ["T1w-brain-template-funcreg", "EPI-brain-template-funcreg"],
        qc_file,
        "template",
    )
    wf.connect_optional(
        strat_pool,
        ["space-template_desc-bold_mask", "space-EPItemplate_desc-bold_mask"],
        resample_bold_mask_to_template,
        "in_file",
    )
    wf.connect(
        [
            (bold_to_T1w_mask, qc_file, [("out_file", "bold2t1w_mask")]),
            (
                resample_bold_mask_to_template,
                qc_file,
                [("out_file", "bold2template_mask")],
            ),
            (
                bids_info,
                qc_file,
                [
                    ("subject", "sub"),
                    ("session", "ses"),
                    ("task", "task"),
                    ("run", "run"),
                ],
            ),
        ]
    )

    return wf, {"space-template_desc-xcp_quality": (qc_file, "qc_file")}
