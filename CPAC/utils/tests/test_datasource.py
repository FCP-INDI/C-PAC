import os
import json
from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util
import pytest

from CPAC.utils.test_resources import setup_test_wf
from CPAC.utils.datasource import match_epi_fmaps


@pytest.mark.skip(reason='needs refactoring')
def test_match_epi_fmaps():

    # good data to use
    s3_prefix = "s3://fcp-indi/data/Projects/HBN/MRI/Site-CBIC/sub-NDARAB708LM5"
    s3_paths = [
        "func/sub-NDARAB708LM5_task-rest_run-1_bold.json",
        "fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.nii.gz",
        "fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.json",
        "fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.nii.gz",
        "fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.json"
    ]

    wf, ds, local_paths = setup_test_wf(s3_prefix, s3_paths,
                                        "test_match_epi_fmaps")

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

    fmap_paths_dct = {"epi_PA":
                          {"scan": local_paths["fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.nii.gz"],
                           "scan_parameters": opposite_pe_params},
                      "epi_AP":
                          {"scan": local_paths["fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.nii.gz"],
                           "scan_parameters": same_pe_params}
                      }

    match_fmaps = \
        pe.Node(util.Function(input_names=['fmap_dct',
                                           'bold_pedir'],
                              output_names=['opposite_pe_epi',
                                            'same_pe_epi'],
                              function=match_epi_fmaps,
                              as_module=True),
                name='match_epi_fmaps')
    match_fmaps.inputs.fmap_dct = fmap_paths_dct
    match_fmaps.inputs.bold_pedir = bold_pedir

    ds.inputs.func_json = func_json
    ds.inputs.opposite_pe_json = opposite_pe_json
    ds.inputs.same_pe_json = same_pe_json

    wf.connect(match_fmaps, 'opposite_pe_epi', ds, 'should_be_dir-PA')
    wf.connect(match_fmaps, 'same_pe_epi', ds, 'should_be_dir-AP')

    wf.run()
