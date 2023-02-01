'''Tests for schema.py'''
import pytest

from CPAC.utils.configuration import Configuration
from voluptuous.error import Invalid


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
