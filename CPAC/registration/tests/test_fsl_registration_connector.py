import pytest
from CPAC.registration.registration import FSL_registration_connector

class AttrDict(dict):
    def __getattr__(self, item):
        value = self[item]
        if isinstance(value, dict):
            return AttrDict(value)
        return value

@pytest.mark.parametrize("sink_native_transforms", [True, False])
def test_fsl_registration_connector(sink_native_transforms):
    wf_name = "test_fsl_registration_connector"
    cfg = AttrDict({
        "registration_workflows": {
            "sink_native_transforms": sink_native_transforms
        }
    })
    _, outputs = FSL_registration_connector(wf_name, cfg, opt="FSL")
    expected_keys = {
        "from-T1w_to-template_mode-image_desc-flirt_xfm",
        "from-template_to-T1w_mode-image_desc-flirt_xfm",
    }
    if sink_native_transforms == True:
        assert expected_keys.issubset(
            outputs.keys()
        ), f"Expected outputs {expected_keys} not found in {outputs.keys()}"
    else:
        assert not expected_keys.intersection(
            outputs.keys()
        ), f"Outputs {expected_keys} should not be present when sink_native_transforms is Off"