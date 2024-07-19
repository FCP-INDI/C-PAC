# Copyright (C) 2019-2024  C-PAC Developers

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
import os

import pytest

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.monitoring.custom_logging import getLogger
from CPAC.utils.test_mocks import configuration_strategy_mock
from CPAC.vmhc.vmhc import vmhc as create_vmhc

logger = getLogger("CPAC.utils.tests")


@pytest.mark.skip(reason="test needs refactoring")
def test_vmhc_ants():
    test_name = "test_vmhc_ants"

    # get the config and strat for the mock
    pipeline_config, strat = configuration_strategy_mock(method="ANTS")
    num_strat = 0

    workflow = pe.Workflow(name=test_name)
    workflow.base_dir = pipeline_config.workingDirectory
    workflow.config["execution"] = {
        "hash_method": "timestamp",
        "crashdump_dir": os.path.abspath(pipeline_config.crashLogDirectory),
    }

    nodes = strat.get_nodes_names()
    logger.info("nodes %s", nodes)

    workflow, strat = create_vmhc(
        workflow,
        num_strat,
        strat,
        pipeline_config,
        output_name=f"vmhc_{num_strat}",
    )

    workflow.run()
