import pytest

from CPAC.registration.registration import ANTs_registration_connector


@pytest.mark.parametrize("sink_native_transforms", ["On", "Off"])
def test_ants_registration_connector(sink_native_transforms):
    cfg = {
        "registration-workflows": {"sink_native_transforms": sink_native_transforms},
    }
    params = {"metric": "MI"}
    _, outputs = ANTs_registration_connector(wf_name="test", cfg=cfg, params=params)
    expected_keys = {
        "from-T1w_to-template_mode-image_desc-initial_xfm",
        "from-T1w_to-template_mode-image_desc-rigid_xfm",
        "from-T1w_to-template_mode-image_desc-affine_xfm",
    }
    if sink_native_transforms:
        assert expected_keys.issubset(
            outputs.keys()
        ), f"Expected outputs {expected_keys} not found in {outputs.keys()}"
    else:
        assert not expected_keys.intersection(
            outputs.keys()
        ), f"Outputs {expected_keys} should not be present when sink_native_transforms is Off"
