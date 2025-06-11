import pytest
from CPAC.registration.registration import FSL_registration_connector
from CPAC.utils.configuration import Configuration
from CPAC.utils.tests.test_utils import check_expected_keys


@pytest.mark.parametrize("sink_native_transforms", [True, False])
def test_fsl_registration_connector(sink_native_transforms):
    wf_name = "test_fsl_registration_connector"
    cfg = Configuration(
        {"registration_workflows": {"sink_native_transforms": sink_native_transforms}}
    )
    _, outputs = FSL_registration_connector(wf_name, cfg, opt="FSL")
    expected_keys = {
        "from-T1w_to-template_mode-image_desc-flirt_xfm",
        "from-template_to-T1w_mode-image_desc-flirt_xfm",
    }
    check_expected_keys(sink_native_transforms, outputs, expected_keys)
