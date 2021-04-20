'''Tests for schema.py'''
import pytest

from CPAC.utils.configuration import Configuration
from CPAC.pipeline.schema import valid_options
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
        assert 'func#motion_estimate_filter_valid_options' in str(e.value)
    else:
        Configuration(d)


@pytest.mark.parametrize('centrality_method', [
    'degree_centrality',
    'eigenvector_centrality',
    'local_functional_connectivity_density'])
@pytest.mark.parametrize(
    'correlation_threshold_option',
    valid_options['centrality']['threshold_options'])
@pytest.mark.parametrize(
    'correlation_threshold', [-2, -1, 0, 0.1, 1, 1.1, 100, 100.1])
def test_centrality_thresholds(
    centrality_method, correlation_threshold_option, correlation_threshold
):
    '''Test threshold range validation:
        0 <= value <= 100 for Sparsity Threshold
        -1 <= value <= 1 otherwise
    '''
    config = {
        'FROM': 'default',
        'network_centrality': {
            centrality_method: {
                'correlation_threshold_option': correlation_threshold_option,
                'correlation_threshold': correlation_threshold
            }
        }
    }

    def out_of_range_assertion(e):
        e = str(e)
        return ('correlation_threshold' in e and 'otherwise' in e) or (
            'value must be at most' in e or 'value must be at least' in e
        )

    if correlation_threshold_option == 'Sparsity threshold':
        if centrality_method == 'local_functional_connectivity_density':
            with pytest.raises(Invalid) as e:
                Configuration(config)
            assert (
                'local_functional_connectivity_density' in str(e.value) or
                out_of_range_assertion(e.value)
            )
            return
        elif not (0 <= correlation_threshold <= 100):
            with pytest.raises(Invalid) as e:
                Configuration(config)
            assert out_of_range_assertion(e.value)
            return
    elif not (-1 <= correlation_threshold <= 1):
        with pytest.raises(Invalid) as e:
            Configuration(config)
        assert out_of_range_assertion(e.value)
        return
    Configuration(config)
