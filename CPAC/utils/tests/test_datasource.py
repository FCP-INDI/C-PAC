# Copyright (C) 2019-2025  C-PAC Developers

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
"""Test datasource utilities."""

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Literal, TypeAlias

from networkx.classes.digraph import DiGraph
import pytest

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.datasource import match_epi_fmaps, match_epi_fmaps_function_node
from CPAC.utils.test_resources import setup_test_wf
from CPAC.utils.utils import PE_DIRECTION
from CPAC.utils.datasource import match_epi_fmaps, match_epi_fmaps_function_node, get_fmap_type


@dataclass
class MatchEpiFmapsInputs:
    """Store test data for `match_epi_fmaps`."""

    bold_pedir: PE_DIRECTION
    epi_fmaps: list[tuple[str, dict[str, Any]]]


def match_epi_fmaps_inputs(
    generate: bool, path: Path
) -> tuple[pe.Workflow, MatchEpiFmapsInputs]:
    """Return inputs for `~CPAC.utils.datasource.match_epi_fmaps`."""
    if generate:
        # good data to use
        s3_prefix = "s3://fcp-indi/data/Projects/HBN/MRI/Site-CBIC/sub-NDARAB708LM5"
        s3_paths = [
            "func/sub-NDARAB708LM5_task-rest_run-1_bold.json",
            "fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.nii.gz",
            "fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.json",
            "fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.nii.gz",
            "fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.json",
        ]

        wf, ds, local_paths = setup_test_wf(
            s3_prefix, s3_paths, "test_match_epi_fmaps", test_dir=str(path)
        )

        opposite_pe_json = local_paths["fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.json"]
        same_pe_json = local_paths["fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.json"]
        func_json = local_paths["func/sub-NDARAB708LM5_task-rest_run-1_bold.json"]

        with open(opposite_pe_json, "r") as f:
            opposite_pe_params = json.load(f)

        with open(same_pe_json, "r") as f:
            same_pe_params = json.load(f)

        with open(func_json, "r") as f:
            func_params = json.load(f)
            bold_pedir = func_params["PhaseEncodingDirection"]

        fmap_paths_dct = {
            "epi_PA": {
                "scan": local_paths["fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.nii.gz"],
                "scan_parameters": opposite_pe_params,
            },
            "epi_AP": {
                "scan": local_paths["fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.nii.gz"],
                "scan_parameters": same_pe_params,
            },
        }
        ds.inputs.func_json = func_json
        ds.inputs.opposite_pe_json = opposite_pe_json
        ds.inputs.same_pe_json = same_pe_json
        return wf, MatchEpiFmapsInputs(
            bold_pedir,
            [
                (scan["scan"], scan["scan_parameters"])
                for scan in fmap_paths_dct.values()
            ],
        )
    _paths = [
        f"{path}/sub-NDARAB514MAJ_dir-AP_acq-fMRI_epi.nii.gz",
        f"{path}/sub-NDARAB514MAJ_dir-PA_acq-fMRI_epi.nii.gz",
    ]
    for _ in _paths:
        Path(_).touch(exist_ok=True)
    return pe.Workflow("test_match_epi_fmaps", path), MatchEpiFmapsInputs(
        "j-",
        [
            (
                _paths[0],
                {
                    "AcquisitionMatrixPE": 84,
                    "BandwidthPerPixelPhaseEncode": 23.81,
                    "BaseResolution": 84,
                    "BodyPartExamined": b"BRAIN",
                    "ConsistencyInfo": b"N4_VE11B_LATEST_20150530",
                    "ConversionSoftware": b"dcm2niix",
                    "ConversionSoftwareVersion": b"v1.0.20171215 GCC4.8.4",
                    "DerivedVendorReportedEchoSpacing": 0.00049999,
                    "DeviceSerialNumber": b"67080",
                    "DwellTime": 2.6e-06,
                    "EchoTime": 0.0512,
                    "EchoTrainLength": 84,
                    "EffectiveEchoSpacing": 0.00049999,
                    "FlipAngle": 90,
                    "ImageOrientationPatientDICOM": [1, 0, 0, 0, 1, 0],
                    "ImageType": ["ORIGINAL", "PRIMARY", "M", "ND", "MOSAIC"],
                    "InPlanePhaseEncodingDirectionDICOM": b"COL",
                    "MRAcquisitionType": b"2D",
                    "MagneticFieldStrength": 3,
                    "Manufacturer": b"Siemens",
                    "ManufacturersModelName": b"Prisma_fit",
                    "Modality": b"MR",
                    "PartialFourier": 1,
                    "PatientPosition": b"HFS",
                    "PercentPhaseFOV": 100,
                    "PhaseEncodingDirection": b"j-",
                    "PhaseEncodingSteps": 84,
                    "PhaseResolution": 1,
                    "PixelBandwidth": 2290,
                    "ProcedureStepDescription": b"CMI_HBN-CBIC",
                    "ProtocolName": b"cmrr_fMRI_DistortionMap_AP",
                    "PulseSequenceDetails": b"%CustomerSeq%_cmrr_mbep2d_se",
                    "ReceiveCoilActiveElements": b"HEA;HEP",
                    "ReceiveCoilName": b"Head_32",
                    "ReconMatrixPE": 84,
                    "RepetitionTime": 5.301,
                    "SAR": 0.364379,
                    "ScanOptions": b"FS",
                    "ScanningSequence": b"EP",
                    "SequenceName": b"epse2d1_84",
                    "SequenceVariant": b"SK",
                    "SeriesDescription": b"cmrr_fMRI_DistortionMap_AP",
                    "ShimSetting": [208, -10464, -5533, 615, -83, -88, 55, 30],
                    "SliceThickness": 2.4,
                    "SliceTiming": [
                        2.64,
                        0,
                        2.7275,
                        0.0875,
                        2.815,
                        0.175,
                        2.9025,
                        0.2625,
                        2.9925,
                        0.3525,
                        3.08,
                        0.44,
                        3.1675,
                        0.5275,
                        3.255,
                        0.615,
                        3.3425,
                        0.7025,
                        3.4325,
                        0.7925,
                        3.52,
                        0.88,
                        3.6075,
                        0.9675,
                        3.695,
                        1.055,
                        3.785,
                        1.1425,
                        3.8725,
                        1.2325,
                        3.96,
                        1.32,
                        4.0475,
                        1.4075,
                        4.135,
                        1.495,
                        4.225,
                        1.5825,
                        4.3125,
                        1.6725,
                        4.4,
                        1.76,
                        4.4875,
                        1.8475,
                        4.575,
                        1.935,
                        4.665,
                        2.0225,
                        4.7525,
                        2.1125,
                        4.84,
                        2.2,
                        4.9275,
                        2.2875,
                        5.015,
                        2.375,
                        5.105,
                        2.4625,
                        5.1925,
                        2.5525,
                    ],
                    "SoftwareVersions": b"syngo_MR_E11",
                    "SpacingBetweenSlices": 2.4,
                    "StationName": b"MRTRIO3TX72",
                    "TotalReadoutTime": 0.0414992,
                    "TxRefAmp": 209.923,
                },
            ),
            (
                _paths[1],
                {
                    "AcquisitionMatrixPE": 84,
                    "BandwidthPerPixelPhaseEncode": 23.81,
                    "BaseResolution": 84,
                    "BodyPartExamined": b"BRAIN",
                    "ConsistencyInfo": b"N4_VE11B_LATEST_20150530",
                    "ConversionSoftware": b"dcm2niix",
                    "ConversionSoftwareVersion": b"v1.0.20171215 GCC4.8.4",
                    "DerivedVendorReportedEchoSpacing": 0.00049999,
                    "DeviceSerialNumber": b"67080",
                    "DwellTime": 2.6e-06,
                    "EchoTime": 0.0512,
                    "EchoTrainLength": 84,
                    "EffectiveEchoSpacing": 0.00049999,
                    "FlipAngle": 90,
                    "ImageOrientationPatientDICOM": [1, 0, 0, 0, 1, 0],
                    "ImageType": ["ORIGINAL", "PRIMARY", "M", "ND", "MOSAIC"],
                    "InPlanePhaseEncodingDirectionDICOM": b"COL",
                    "MRAcquisitionType": b"2D",
                    "MagneticFieldStrength": 3,
                    "Manufacturer": b"Siemens",
                    "ManufacturersModelName": b"Prisma_fit",
                    "Modality": b"MR",
                    "PartialFourier": 1,
                    "PatientPosition": b"HFS",
                    "PercentPhaseFOV": 100,
                    "PhaseEncodingDirection": b"j",
                    "PhaseEncodingSteps": 84,
                    "PhaseResolution": 1,
                    "PixelBandwidth": 2290,
                    "ProcedureStepDescription": b"CMI_HBN-CBIC",
                    "ProtocolName": b"cmrr_fMRI_DistortionMap_PA",
                    "PulseSequenceDetails": b"%CustomerSeq%_cmrr_mbep2d_se",
                    "ReceiveCoilActiveElements": b"HEA;HEP",
                    "ReceiveCoilName": b"Head_32",
                    "ReconMatrixPE": 84,
                    "RepetitionTime": 5.301,
                    "SAR": 0.364379,
                    "ScanOptions": b"FS",
                    "ScanningSequence": b"EP",
                    "SequenceName": b"epse2d1_84",
                    "SequenceVariant": b"SK",
                    "SeriesDescription": b"cmrr_fMRI_DistortionMap_PA",
                    "ShimSetting": [208, -10464, -5533, 615, -83, -88, 55, 30],
                    "SliceThickness": 2.4,
                    "SliceTiming": [
                        2.64,
                        0,
                        2.73,
                        0.09,
                        2.8175,
                        0.1775,
                        2.905,
                        0.265,
                        2.9925,
                        0.3525,
                        3.08,
                        0.44,
                        3.17,
                        0.53,
                        3.2575,
                        0.6175,
                        3.345,
                        0.705,
                        3.4325,
                        0.7925,
                        3.52,
                        0.88,
                        3.61,
                        0.97,
                        3.6975,
                        1.0575,
                        3.785,
                        1.145,
                        3.8725,
                        1.2325,
                        3.9625,
                        1.32,
                        4.05,
                        1.41,
                        4.1375,
                        1.4975,
                        4.225,
                        1.585,
                        4.3125,
                        1.6725,
                        4.4025,
                        1.76,
                        4.49,
                        1.85,
                        4.5775,
                        1.9375,
                        4.665,
                        2.025,
                        4.7525,
                        2.1125,
                        4.8425,
                        2.2,
                        4.93,
                        2.29,
                        5.0175,
                        2.3775,
                        5.105,
                        2.465,
                        5.1925,
                        2.5525,
                    ],
                    "SoftwareVersions": b"syngo_MR_E11",
                    "SpacingBetweenSlices": 2.4,
                    "StationName": b"MRTRIO3TX72",
                    "TotalReadoutTime": 0.0414992,
                    "TxRefAmp": 209.923,
                },
            ),
        ],
    )


RunType: TypeAlias = Literal["nipype"] | Literal["direct"]
Direction: TypeAlias = Literal["opposite"] | Literal["same"]


@pytest.mark.parametrize("generate", [True, False])
def test_match_epi_fmaps(generate: bool, tmp_path: Path) -> None:
    """Test `~CPAC.utils.datasource.match_epi_fmaps`."""
    wf, data = match_epi_fmaps_inputs(generate, tmp_path)

    match_fmaps = match_epi_fmaps_function_node()
    match_fmaps.inputs.bold_pedir = data.bold_pedir
    match_fmaps.inputs.epi_fmap_one = data.epi_fmaps[0][0]
    match_fmaps.inputs.epi_fmap_params_one = data.epi_fmaps[0][1]
    match_fmaps.inputs.epi_fmap_two = data.epi_fmaps[1][0]
    match_fmaps.inputs.epi_fmap_params_two = data.epi_fmaps[1][1]

    wf.add_nodes([match_fmaps])

    graph: DiGraph = wf.run()
    result = list(graph.nodes)[-1].run()
    str_outputs: dict[RunType, dict[Direction, str]] = {
        "nipype": {
            "opposite": result.outputs.opposite_pe_epi,
            "same": result.outputs.same_pe_epi,
        },
        "direct": {},
    }
    path_outputs: dict[RunType, dict[Direction, Path]] = {"nipype": {}, "direct": {}}
    str_outputs["direct"]["opposite"], str_outputs["direct"]["same"] = match_epi_fmaps(
        data.bold_pedir,
        data.epi_fmaps[0][0],
        data.epi_fmaps[0][1],
        data.epi_fmaps[1][0],
        data.epi_fmaps[1][1],
    )
    directions: list[Direction] = ["opposite", "same"]
    runtypes: list[RunType] = ["nipype", "direct"]
    for direction in directions:
        for runtype in runtypes:
            path_outputs[runtype][direction] = Path(str_outputs[runtype][direction])
            assert path_outputs[runtype][direction].exists()
        assert (
            path_outputs["nipype"][direction].name
            == path_outputs["direct"][direction].name
        )


@pytest.mark.parametrize(
    "metadata, expected_type",
    [
        # Case 1: Phase-difference map (phasediff) - REQUIRED: EchoTime1 and EchoTime2
        ({"EchoTime1": 0.00600, "EchoTime2": 0.00746}, "phasediff"),
        ({"EchoTime1": 0.004, "EchoTime2": 0.006}, "phasediff"),
        
        # Case 2: Single phase map (phase) - REQUIRED: EchoTime
        ({"EchoTime": 0.00746}, "phase"),
        ({"EchoTime": 0.004}, "phase"),
        
        # Case 3: Direct field mapping (fieldmap) - REQUIRED: Units
        ({"Units": "rad/s"}, "fieldmap"),
        ({"Units": "Hz"}, "fieldmap"),
        ({"Units": "hz"}, "fieldmap"),
        ({"Units": "T"}, "fieldmap"),
        ({"Units": "Tesla"}, "fieldmap"),
        ({"Units": "hertz"}, "fieldmap"),
        
        # Case 4: EPI field maps (epi) - REQUIRED: PhaseEncodingDirection
        ({"PhaseEncodingDirection": "j-"}, "epi"),
        ({"PhaseEncodingDirection": "j"}, "epi"),
        ({"PhaseEncodingDirection": "i"}, "epi"),
        ({"PhaseEncodingDirection": "i-"}, "epi"),
        ({"PhaseEncodingDirection": "k"}, "epi"),
        ({"PhaseEncodingDirection": "k-"}, "epi"),
        
        # Edge cases and invalid inputs
        ({}, None),  # Empty metadata
        ({"SomeOtherField": "value"}, None),  # Irrelevant metadata
        ({"Units": "invalid_unit"}, None),  # Invalid units
        ({"PhaseEncodingDirection": "invalid"}, None),  # Invalid PE direction
        ({"EchoTime1": 0.006}, None),  # Only EchoTime1 without EchoTime2
        ({"EchoTime2": 0.006}, None),  # Only EchoTime2 without EchoTime1
        
        # Priority testing - phasediff should take precedence
        ({"EchoTime1": 0.006, "EchoTime2": 0.007, "EchoTime": 0.006}, "phasediff"),
        ({"EchoTime1": 0.006, "EchoTime2": 0.007, "Units": "Hz"}, "phasediff"),
        ({"EchoTime1": 0.006, "EchoTime2": 0.007, "PhaseEncodingDirection": "j-"}, "phasediff"),
        
        # Phase should take precedence over fieldmap and epi
        ({"EchoTime": 0.006, "Units": "Hz"}, "phase"),
        ({"EchoTime": 0.006, "PhaseEncodingDirection": "j-"}, "phase"),
        
        # Fieldmap should take precedence over epi
        ({"Units": "Hz", "PhaseEncodingDirection": "j-"}, "fieldmap"),
        
        # Test with optional fields that might be present (but shouldn't affect detection)
        ({"EchoTime1": 0.006, "EchoTime2": 0.007, "IntendedFor": "bids::sub-01/func/sub-01_task-motor_bold.nii.gz"}, "phasediff"),
        ({"Units": "rad/s", "IntendedFor": "bids::sub-01/func/sub-01_task-motor_bold.nii.gz"}, "fieldmap"),
        ({"PhaseEncodingDirection": "j-", "TotalReadoutTime": 0.095}, "epi"),
    ]
)
def test_get_fmap_type_dict_input(metadata: dict, expected_type: str | None) -> None:
    """Test `get_fmap_type` with dictionary input using only required BIDS fields."""
    result = get_fmap_type(metadata)
    assert result == expected_type


def test_get_fmap_type_real_world_examples() -> None:
    """Test `get_fmap_type` with realistic BIDS metadata examples (required fields only)."""
    # Real-world phasediff example (only required fields)
    phasediff_metadata = {
        "EchoTime1": 0.00600,
        "EchoTime2": 0.00746,
        # Optional fields that might be present:
        "IntendedFor": ["bids::sub-01/func/sub-01_task-motor_bold.nii.gz"]
    }
    assert get_fmap_type(phasediff_metadata) == "phasediff"
    
    # Real-world fieldmap example (only required fields)
    fieldmap_metadata = {
        "Units": "rad/s",
        # Optional fields that might be present:
        "IntendedFor": "bids::sub-01/func/sub-01_task-motor_bold.nii.gz"
    }
    assert get_fmap_type(fieldmap_metadata) == "fieldmap"
    
    # Real-world EPI example (only required fields)
    epi_metadata = {
        "PhaseEncodingDirection": "j-",
        # Optional fields that might be present:
        "TotalReadoutTime": 0.095,
        "IntendedFor": "bids::sub-01/func/sub-01_task-motor_bold.nii.gz"
    }
    assert get_fmap_type(epi_metadata) == "epi"
    
    # Real-world phase example (only required fields)
    phase_metadata = {
        "EchoTime": 0.00746
    }
    assert get_fmap_type(phase_metadata) == "phase"