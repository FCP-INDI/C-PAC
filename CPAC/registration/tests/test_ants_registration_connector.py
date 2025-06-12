import pytest
from CPAC.registration.registration import ANTs_registration_connector
from CPAC.utils.configuration import Configuration
from CPAC.utils.tests.test_utils import check_expected_keys


@pytest.mark.parametrize("sink_native_transforms", [True, False])
def test_ants_registration_connector(sink_native_transforms):
    wf_name = "test_ants_registration_connector"
    cfg = Configuration(
        {
            "pipeline_setup": {"system_config": {"num_ants_threads": 1}},
            "registration_workflows": {
                "sink_native_transforms": sink_native_transforms,
                "anatomical_registration": {
                    "reg_with_skull": True,
                    "registration": {"ANTs": {"use_lesion_mask": False}},
                },
            },
        }
    )
    params = {"metric": "MI"}
    _, outputs = ANTs_registration_connector(wf_name, cfg=cfg, params=params)
    expected_keys = {
        "from-T1w_to-template_mode-image_desc-initial_xfm",
        "from-T1w_to-template_mode-image_desc-rigid_xfm",
        "from-T1w_to-template_mode-image_desc-affine_xfm",
    }
    check_expected_keys(sink_native_transforms, outputs, expected_keys)
