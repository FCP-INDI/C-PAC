import pytest
from types import SimpleNamespace
from CPAC.registration.registration import FSL_registration_connector

@pytest.mark.parametrize("sink_native_transforms", ["On", "Off"])
def test_fsl_registration_connector(sink_native_transforms):
    wf_name = "test_fsl_registration_connector"
    cfg = SimpleNamespace(
        registration_workflows=SimpleNamespace(
            sink_native_transforms=sink_native_transforms
        )
    )
    _, outputs = FSL_registration_connector(wf_name, cfg)
    expected_outputs = {
        "from-T1w_to-template_mode-image_desc-linear_xfm",
        "from-template_to-T1w_mode-image_desc-linear_xfm",
    }
    if sink_native_transforms == "On":
        assert expected_outputs.issubset(outputs.keys()), (
            f"Expected outputs {expected_outputs} not found in {outputs.keys()}"
        )
    else:
        assert not expected_outputs.intersection(outputs.keys()), (
            f"Outputs {expected_outputs} should not be present when sink_native_transforms is Off"
        )