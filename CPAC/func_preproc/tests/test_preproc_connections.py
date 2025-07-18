# Copyright (C) 2023-2025  C-PAC Developers

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
"""Test graph connections for functional preprocessing."""

from itertools import product
from random import sample
import re
from typing import Literal, Optional

import pytest
from traits.trait_base import Undefined
from voluptuous.error import Invalid, RequiredFieldInvalid
from nipype.interfaces.utility import Function as NipypeFunction
from nipype.pipeline.engine import Node, Workflow as NipypeWorkflow

from CPAC.func_preproc.pipeline import stack_func_preproc_blocks
from CPAC.nuisance.nuisance import choose_nuisance_blocks
from CPAC.pipeline.cpac_pipeline import connect_pipeline
from CPAC.pipeline.engine import ResourcePool
from CPAC.pipeline.nipype_pipeline_engine import Workflow
from CPAC.utils.configuration import Configuration
from CPAC.utils.interfaces.function import Function as CpacFunction
from CPAC.utils.test_init import create_dummy_node

_FILTERS = [
    {
        "filter_type": "notch",
        "filter_order": 4,
        "center_frequency": 0.31,
        "filter_bandwidth": 0.12,
    },
    {"filter_type": "lowpass", "filter_order": 4, "lowpass_cutoff": 0.0032},
]
_PRE_RESOURCES = [
    "bold",
    "label-CSF_desc-eroded_mask",
    "label-CSF_desc-preproc_mask",
    "label-CSF_mask",
    "label-GM_desc-eroded_mask",
    "label-GM_desc-preproc_mask",
    "label-GM_mask",
    "label-WM_desc-eroded_mask",
    "label-WM_desc-preproc_mask",
    "label-WM_mask",
    "lateral-ventricles-mask",
    "space-T1w_desc-brain_mask",
    "space-T1w_desc-eroded_mask",
    "space-bold_desc-brain_mask",
    "TR",
    # "tr",
    # "tpattern",
    # "scan",
    # "subject",
    # "desc-brain_T1w",
    # "from-T1w_to-template_mode-image_desc-linear_xfm",
    # "from-bold_to-T1w_mode-image_desc-linear_xfm",
    # "from-template_to-T1w_mode-image_desc-linear_xfm",
]

NUM_TESTS = 48  # number of parameterizations to run for many-parameter tests


def _filter_assertion_message(
    subwf: NipypeWorkflow, is_filtered: list[Node], should_be_filtered: bool
) -> str:
    if is_filtered and not should_be_filtered:
        return (
            f'{subwf.name} is filtered by '
            f'{" & ".join([str(node.name) for node in is_filtered])} and should '
            'not be'
        )
    return f"{subwf.name} is not filtered and should be"


_PARAMS = {  # for test_motion_filter_connections
    "motion_estimation_timing": ["before", "after", None],
    "filters": [[_FILTERS[0]], [_FILTERS[1]], _FILTERS],
    "motion_correction": [["mcflirt"], ["3dvolreg"], ["mcflirt", "3dvolreg"]],
    "pre_resources": [
        _PRE_RESOURCES,
        ["desc-movementParameters_motion", *_PRE_RESOURCES],
    ],
    "regtool": ["ANTs", "FSL"],
    "run": [True, False, [True, False]],
}


@pytest.mark.parametrize(
    ",".join(_PARAMS.keys()),  # run n=NUM_TESTS subset
    sample(list(product(*_PARAMS.values())), NUM_TESTS),
)
def test_motion_filter_connections(
    run: bool | list[bool],
    filters: list[dict],
    regtool: list[str],
    motion_estimation_timing: Optional[Literal["before", "after"]],
    pre_resources: list[str],
    motion_correction: list[list[str]],
) -> None:
    """Test that appropriate connections occur vis-à-vis motion filters."""
    if isinstance(motion_correction, list) and len(motion_correction) != 1:
        # Until https://github.com/FCP-INDI/C-PAC/issues/1935 is resolved
        with pytest.raises(Invalid, match="FCP-INDI/C-PAC/issues/1935"):
            c = Configuration(
                {
                    "functional_preproc": {
                        "motion_estimates_and_correction": {
                            "motion_correction": {"using": motion_correction},
                            "motion_estimation_timing": motion_estimation_timing,
                        }
                    }
                }
            )
        return
    if motion_estimation_timing is None:
        with pytest.raises(
            RequiredFieldInvalid, match="motion_estimates_and_correction"
        ):
            c = Configuration(
                {
                    "functional_preproc": {
                        "motion_estimates_and_correction": {
                            "run": True,
                            "motion_correction": {"using": motion_correction},
                            "motion_estimation_timing": motion_estimation_timing,
                        }
                    }
                }
            )
        return
    # parameterized Configuration
    c = Configuration(
        {
            "functional_preproc": {
                "motion_estimates_and_correction": {
                    "motion_correction": {"using": motion_correction},
                    "motion_estimation_timing": motion_estimation_timing,
                    "motion_estimate_filter": {"run": run, "filters": filters},
                    "run": True,
                },
                "run": True,
            },
            "nuisance_corrections": {
                "2-nuisance_regression": {
                    "Regressors": [
                        {
                            "Name": "aCompCor, GSR, no censor",
                            "Motion": {
                                "include_delayed": True,
                                "include_squared": True,
                                "include_delayed_squared": True,
                            },
                            "aCompCor": {
                                "summary": {"method": "DetrendPC", "components": 5},
                                "tissues": ["WhiteMatter", "CerebrospinalFluid"],
                                "extraction_resolution": 3,
                            },
                            "GlobalSignal": {"summary": "Mean"},
                            "PolyOrt": {"degree": 2},
                            "Bandpass": {
                                "bottom_frequency": 0.01,
                                "top_frequency": 0.1,
                            },
                        }
                    ]
                }
            },
        }
    )
    # resource for intial inputs
    before_this_test = create_dummy_node("created_before_this_test", pre_resources)
    rpool = ResourcePool(cfg=c)
    for resource in pre_resources:
        if resource.endswith("xfm"):
            node_name = f"created_before_this_test_{regtool}"
        elif resource == "desc-movementParameters_motion":
            node_name = f"created_before_this_test_{motion_correction}"
        else:
            node_name = "created_before_this_test"
        rpool.set_data(resource, before_this_test, resource, {}, "", node_name)
    # set up blocks
    pipeline_blocks = stack_func_preproc_blocks({}, c, rpool)
    # Nuisance Correction
    generate_only = (
        True not in c["nuisance_corrections", "2-nuisance_regression", "run"]
    )
    if not rpool.check_rpool("desc-cleaned_bold"):
        pipeline_blocks += choose_nuisance_blocks(c, generate_only)
    wf = Workflow(re.sub(r"[\[\]\-\:\_ \'\",]", "", str(rpool)))
    connect_pipeline(wf, c, rpool, pipeline_blocks)
    # Check that filtering is happening as expected
    filter_switch_key = [
        "functional_preproc",
        "motion_estimates_and_correction",
        "motion_estimate_filter",
        "run",
    ]
    if c.switch_is_on(filter_switch_key, exclusive=True):
        assert all(
            strat.filtered_movement
            for strat in rpool.get_strats(["desc-movementParameters_motion"]).values()
        )
    elif c.switch_is_off(filter_switch_key, exclusive=True):
        assert not any(
            strat.filtered_movement
            for strat in rpool.get_strats(["desc-movementParameters_motion"]).values()
        )
    elif c.switch_is_on_off(filter_switch_key):
        assert any(
            strat.filtered_movement
            for strat in rpool.get_strats(["desc-movementParameters_motion"]).values()
        )
        if (
            "mcflirt"
            in c[
                "functional_preproc",
                "motion_estimates_and_correction",
                "motion_correction",
                "using",
            ]
            and "desc-movementParameters_motion" not in pre_resources
        ):
            # Only for [On, Off] + mcflirt, we should have at least one of each
            assert {
                wf.get_node(nodename).inputs.calc_from
                for nodename in wf.list_node_names()
                if nodename.endswith(".calculate_FDJ")
                and wf.get_node(nodename).inputs.calc_from is not Undefined
            } == {"affine", "rms"}
    regressor_subwfs = [
        wf.get_node(nodename[:-26])
        for nodename in wf.list_node_names()
        if nodename.endswith("build_nuisance_regressors")
    ]
    for subwf in regressor_subwfs:
        # a motion filter is an input to the nuisance regressor subworkflow
        is_filtered = []
        # a motion filter should be an input to the regressor subworkflow
        should_be_filtered = "_filt-" in subwf.name and "_filt-none" not in subwf.name
        for u, v in wf._graph.edges:  # pylint: disable=invalid-name,protected-access
            if (
                v == subwf
                and hasattr(u, "interface")
                and isinstance(u.interface, (NipypeFunction, CpacFunction))
                and "notch_filter_motion" in u.interface.inputs.function_str
            ):
                is_filtered.append(u)
        assert bool(is_filtered) == should_be_filtered, _filter_assertion_message(
            subwf, is_filtered, should_be_filtered
        )
