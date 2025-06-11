import pytest
from CPAC.registration.registration import ANTs_registration_connector

class AttrDict(dict):
    def __getattr__(self, item):
        value = self[item]
        if isinstance(value, dict):
            return AttrDict(value)
        return value

@pytest.mark.parametrize("sink_native_transforms", [True, False])
def test_ants_registration_connector(sink_native_transforms):
    wf_name = "test_ants_registration_connector"
    cfg = AttrDict({
        "registration-workflows": {"sink_native_transforms": sink_native_transforms},
        "pipeline_setup": {
            "system_config": {
                "num_ants_threads": 1
            }
        },
        "registration_workflows": {
            "sink_native_transforms": sink_native_transforms,  
            "anatomical_registration": {
                "reg_with_skull": True,
                "registration": {
                    "ANTs": {
                        "use_lesion_mask": False
                    }
                }
            }
        }
    })
    params = {"metric": "MI"}
    _, outputs = ANTs_registration_connector(wf_name, cfg=cfg, params=params)
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