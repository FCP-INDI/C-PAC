import pytest
from types import SimpleNamespace
import CPAC.registration.ants_registration_connector as ants_registration_connector

@pytest.fixture
def dummy_module(monkeypatch):
    class DummyNode:
        def __init__(self):
            self.inputs = SimpleNamespace(inputspec=SimpleNamespace())
            self.outputspec = SimpleNamespace(
                ants_initial_xfm="initial.mat",
                ants_rigid_xfm="rigid.mat",
                ants_affine_xfm="affine.mat",
            )
    def dummy_create_wf(name):
        return DummyNode()
    monkeypatch.setattr(ants_registration_connector, 'create_wf_calculate_ants_warp', dummy_create_wf)
    monkeypatch.setattr(ants_registration_connector, 'check_transforms', lambda x: (x, len(x)))
    monkeypatch.setattr(ants_registration_connector, 'generate_inverse_transform_flags', lambda x: [True]*len(x))
    return ants_registration_connector

def build_cfg(sink_native_transforms=True):
    cfg = SimpleNamespace()
    cfg.FROM = 'default'
    cfg.registration_workflows = {
        'sink_native_transforms': sink_native_transforms
    }
    return cfg

def test_sink_native_transforms_outputs(dummy_module):
    connector = dummy_module
    cfg = build_cfg(sink_native_transforms=True)
    _, outputs = connector.ANTs_registration_connector(
        wf_name='test', cfg=cfg
    )
    expected_keys = [
        'from-T1w_to-template_mode-image_desc-initial_xfm',
        'from-T1w_to-template_mode-image_desc-rigid_xfm',
        'from-T1w_to-template_mode-image_desc-affine_xfm',
    ]
    for key in expected_keys:
        assert key in outputs

def test_no_sink_native_transforms(dummy_module):
    connector = dummy_module
    cfg = build_cfg(sink_native_transforms=False)
    _, outputs = connector.ANTs_registration_connector(
        wf_name='test', cfg=cfg
    )
    for key in [
        'from-T1w_to-template_mode-image_desc-initial_xfm',
        'from-T1w_to-template_mode-image_desc-rigid_xfm',
        'from-T1w_to-template_mode-image_desc-affine_xfm',
    ]:
        assert key not in outputs