import numpy as np
import pytest

from CPAC.isc.pipeline import create_isc, create_isfc


@pytest.mark.skip(reason="needs refactored")
def test_pipeline_isc():
    wf = create_isc()

    wf.inputs.inputspec.D = np.random.uniform(size=(3, 100, 20))
    wf.inputs.inputspec.permutations = 1
    wf.inputs.inputspec.collapse_subj = True
    wf.inputs.inputspec.random_state = 42
    wf.run(plugin="Linear")


@pytest.mark.skip(reason="needs refactored")
def test_pipeline_isfc():
    wf = create_isfc()

    wf.inputs.inputspec.D = np.random.uniform(size=(3, 100, 20))
    wf.inputs.inputspec.permutations = 10
    wf.inputs.inputspec.collapse_subj = True
    wf.inputs.inputspec.random_state = 42
    wf.run(plugin="Linear")
