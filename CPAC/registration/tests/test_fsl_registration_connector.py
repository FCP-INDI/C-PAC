import pytest
from types import SimpleNamespace
from CPAC.registration.registration import FSL_registration_connector

@pytest.fixture
def dummy_module(monkeypatch):
    class DummyNode:
        def __init__(self):
            self.inputs = SimpleNamespace(inputspec=SimpleNamespace())
            self.outputspec = SimpleNamespace(
                linear_xfm="linear.mat",
                invlinear_xfm="invlinear.mat",
            )
    def dummy_create_linear(name):
        return DummyNode()
    def dummy_create_nonlinear(name):
        return DummyNode()
    monkeypatch.setattr(
        "CPAC.registration.registration.create_fsl_flirt_linear_reg", dummy_create_linear
    )
    monkeypatch.setattr(
        "CPAC.registration.registration.create_fsl_fnirt_nonlinear_reg_nhp", dummy_create_nonlinear
    )
    return FSL_registration_connector

def build_cfg(sink_native_transforms=True):
    cfg = SimpleNamespace()
    cfg.registration_workflows = {
        'sink_native_transforms': sink_native_transforms
    }
    return cfg

def test_sink_native_transforms_outputs(dummy_module):
    connector = dummy_module
    cfg = build_cfg(sink_native_transforms=True)
    _, outputs = connector(
        wf_name='test', cfg=cfg, orig="T1w", opt="FSL"
    )
    expected_keys = [
        'from-T1w_to-template_mode-image_desc-linear_xfm',
        'from-template_to-T1w_mode-image_desc-linear_xfm',
    ]
    for key in expected_keys:
        assert key in outputs