import pytest

from CPAC.registration.registration import FSL_registration_connector


@pytest.mark.parametrize("sink_native_transforms", ["On", "Off"])
def test_fsl_registration_connector(sink_native_transforms):
    wf_name = "test_fsl_registration_connector"
    cfg = {
        "registration-workflows": {"sink_native_transforms": sink_native_transforms},
    }
    _, outputs = FSL_registration_connector(wf_name, cfg)
    if sink_native_transforms == "On":
        expected_outputs = {
            "from-T1w_to-template_mode-image_desc-linear_xfm",
            "from-template_to-T1w_mode-image_desc-linear_xfm",
        }
        assert expected_outputs.issubset(
            outputs.keys()
        ), f"Expected outputs {expected_outputs} not found in {outputs.keys()}"
    else:
        # Adjust this set based on what outputs should be present when 'Off'
        not_expected_outputs = {
            "from-T1w_to-template_mode-image_desc-linear_xfm",
            "from-template_to-T1w_mode-image_desc-linear_xfm",
        }
        assert not not_expected_outputs.intersection(
            outputs.keys()
        ), f"Outputs {not_expected_outputs} should not be present when sink_native_transforms is Off"
