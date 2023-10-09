# Copyright (C) 2015-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
from itertools import combinations
from pathlib import Path
import pytest
from CPAC.network_centrality.network_centrality import create_centrality_wf
from CPAC.pipeline.schema import valid_options
from CPAC.utils.interfaces.afni import AFNI_SEMVER
from CPAC.utils.typing import LIST

_DATA_DIR = Path(__file__).parent / 'data'
"""Path to test data directory"""


@pytest.mark.parametrize('method_option',
                         valid_options['centrality']['method_options'])
@pytest.mark.parametrize('weight_options',
    [*combinations(valid_options['centrality']['weight_options'], 1),
     *combinations(valid_options['centrality']['weight_options'], 2)])
@pytest.mark.parametrize('threshold_option',
                          valid_options['centrality']['threshold_options'])
@pytest.mark.parametrize('threshold', [0.001, 0.6])
@pytest.mark.skipif(AFNI_SEMVER == '0.0.0', reason='AFNI not installed')
def test_create_centrality_wf(method_option: str, weight_options: LIST[str],
                              threshold_option: str,
                              threshold: float, tmpdir: Path) -> None:
    '''Integration test of
    ~CPAC.network_centrality.network_centrality.create_centrality_wf'''
    wf_name = (f'test_{method_option[0]}'
               f'{"".join([_[0] for _ in weight_options])}'
               f'{threshold_option[0]}{threshold}'.replace('.', 'p'))
    if (method_option == 'local_functional_connectivity_density' and
        threshold_option == 'Sparsity threshold'
    ):
        with pytest.raises(ValueError):
            centrality_wf = create_centrality_wf(wf_name, method_option,
                weight_options, threshold_option, threshold, base_dir=tmpdir)
        return
    if method_option == 'eigenvector_centrality':
        centrality_wf = create_centrality_wf(wf_name, method_option,
            weight_options, threshold_option, threshold, memory_gb=3.0,
            base_dir=tmpdir)
    else:
        centrality_wf = create_centrality_wf(wf_name, method_option,
            weight_options, threshold_option, threshold, base_dir=tmpdir)
    centrality_wf.inputs.inputspec.in_file = (_DATA_DIR /
                                              'in_file.nii.gz').absolute()
    centrality_wf.inputs.inputspec.template = (_DATA_DIR /
                                               'template.nii.gz').absolute()
    centrality_wf.run()
