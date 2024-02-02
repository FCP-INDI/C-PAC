import os

import pytest

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.test_mocks import configuration_strategy_mock
from CPAC.vmhc.vmhc import vmhc as create_vmhc


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

    strat.get_nodes_names()

    workflow, strat = create_vmhc(
        workflow,
        num_strat,
        strat,
        pipeline_config,
        output_name=f"vmhc_{num_strat}",
    )

    workflow.run()
