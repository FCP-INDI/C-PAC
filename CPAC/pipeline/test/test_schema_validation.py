'''Tests for schema.py'''
from itertools import combinations
import pytest
from voluptuous.error import ExclusiveInvalid, Invalid
from CPAC.utils.configuration import Configuration


@pytest.mark.parametrize('run_value', [
    True, False, [True], [False], [True, False], [False, True], None, []])
def test_motion_estimates_and_correction(run_value):
    '''Test that any truthy forkable option for 'run' throws the custom
    human-readable exception for an invalid motion_estimate_filter.
    '''
    d = {
        'FROM': 'default',
        'functional_preproc': {'motion_estimates_and_correction': {
             'motion_estimate_filter': {'run': run_value,
                                        'filter_type': 'notch',
                                        'filter_order': 0,
                                        'breathing_rate_min': None,
                                        'breathing_rate_max': 101.5}}}
    }
    if bool(run_value) and run_value not in [[False], []]:
        with pytest.raises(Invalid) as e:
            Configuration(d)
        assert "func#motion_estimate_filter_valid_options" in str(e.value)
    else:
        Configuration(d)


@pytest.mark.parametrize('registration_using',
                         [list(combo) for _ in [list(combinations(
                          ['ANTS', 'FSL', 'FSL-linear'], i)) for i in
                          range(1, 4)] for combo in _])
def test_single_step_vs_registration(registration_using):
    '''Test that single-step resampling requires ANTS registration'''
    # pylint: disable=invalid-name
    d = {'registration_workflows': {'anatomical_registration': {
        'registration': {'using': registration_using}},
                                    'functional_registration':  {
        'func_registration_to_template': {'apply_transform': {
            'using': 'single_step_resampling_from_stc'}}}}}
    if registration_using == ['ANTS']:
        Configuration(d)  # validates without exception
    else:
        with pytest.raises(ExclusiveInvalid) as e:
            Configuration(d)
        assert "requires ANTS registration" in str(e.value)


def test_pipeline_name():
    '''Test that pipeline_name sucessfully sanitizes'''
    c = Configuration({'pipeline_setup': {'pipeline_name': ':va:lid    name'}})
    assert c['pipeline_setup', 'pipeline_name'] == 'valid_name'
