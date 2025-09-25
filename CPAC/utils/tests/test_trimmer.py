# Copyright (C) 2020-2025  C-PAC Developers

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
"""Test The Trimmer."""

from copy import copy
import tempfile

import pytest
import yaml


def accept_all(object, name, value):
    return value


@pytest.mark.skip(reason="needs refactored")
def test_trimmer():
    """Test The Trimmer."""
    import os

    from CPAC.pipeline.cpac_pipeline import build_workflow
    from CPAC.resources.configs import CONFIGS_PATH
    from CPAC.utils.configuration import Configuration
    from CPAC.utils.trimmer import (
        compute_datasink_dirs,
        expand_workflow,
        is_datasink,
        the_trimmer,
    )

    pipe_config = CONFIGS_PATH / "pipeline_config_template.yml"
    data_config = CONFIGS_PATH / "data_config_S3-BIDS-ABIDE.yml"

    data_config = yaml.safe_load(data_config.open("r"))
    sub_dict = data_config[0]

    c = Configuration(yaml.safe_load(pipe_config.open("r")))
    temp_dir = tempfile.mkdtemp()
    c.logDirectory = temp_dir
    c.workingDirectory = temp_dir
    c.outputDirectory = temp_dir
    c.crashLogDirectory = temp_dir

    # Disable functional, let only the anatomical workflow
    c_anatomical = copy(c)
    c_anatomical.runFunctional = [0]

    wf, _, _ = build_workflow(sub_dict["subject_id"], sub_dict, c_anatomical)

    # Create fake files to trick THE TRIMMER
    exec_graph = expand_workflow(wf)
    datasinks = [n for n in exec_graph.nodes() if is_datasink(n)]
    anat_derivatives = {}
    for datasink in datasinks:
        paths = compute_datasink_dirs(exec_graph, datasink)
        anat_derivatives.update(paths)
        for (node, derivative), path in paths.items():
            os.makedirs(path)
            open(os.path.join(path, "%s.txt" % derivative), "a").close()

    # Enable functional, so the workflow should only run this
    # and enable trimming
    c_functional = copy(c)
    c_functional.runFunctional = [1]

    wf, _, _ = build_workflow(sub_dict["subject_id"], sub_dict, c_functional)
    exec_wf, _ = the_trimmer(wf)
    exec_graph = exec_wf._graph

    datasinks = [n for n in exec_graph.nodes() if is_datasink(n)]
    func_derivatives = {}
    for datasink in datasinks:
        paths = compute_datasink_dirs(exec_graph, datasink)
        func_derivatives.update(paths)

    # Assert that the functional pipeline remove all the anatomical nodes,
    # as they were already computed
    assert (
        set(func_derivatives.keys()).intersection(set(anat_derivatives.keys())) == set()
    )
