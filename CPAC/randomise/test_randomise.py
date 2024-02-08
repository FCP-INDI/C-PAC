# Copyright (C) 2018-2024  C-PAC Developers

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
"""Test randomise."""
import os

import pytest
import nipype.interfaces.io as nio

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nipype_pipeline_engine.plugins import MultiProcPlugin


@pytest.mark.skip(reason="test needs refactoring")
@pytest.mark.parametrize(
    "inputs", [["subjects", "design_matrix_file", "contrast_file", "permutations"]]
)
def test_run_randomize(inputs, output_dir=None, run=True):
    """Test randomize."""
    from . import pipeline

    randomise_workflow = pe.Workflow(name="preproc")
    if output_dir is None:
        output_dir = ""

    workflow_dir = os.path.join(output_dir, "randomise_results")
    randomise_workflow.base_dir = workflow_dir

    # resource_pool = {}

    num_of_cores = 1

    t_node = pipeline.create_randomise()
    t_node.inputs.inputspec.subjects = ""
    t_node.inputs.inputspec.design_matrix_file = ""
    t_node.inputs.inputspec.constrast_file = ""
    # t_node.inputs.inputspec.f_constrast_file = ''
    t_node.inputs.inputspec.permutations = 5000
    # t_node.inputs.inputspec.mask = '' #

    dataSink = pe.Node(nio.DataSink(), name="dataSink_file")
    dataSink.inputs.base_directory = workflow_dir
    randomise_workflow.connect(t_node, "outputspec.index_file", dataSink, "index_file")
    # randomise_workflow.connect(t_node, 'outputspec.thresh_out',
    #                            dataSink,'threshold_file')
    randomise_workflow.connect(
        t_node, "outputspec.localmax_txt_file", dataSink, "localmax_txt_file"
    )
    randomise_workflow.connect(
        t_node, "outputspec.localmax_vol_file", dataSink, "localmax_vol_file"
    )
    randomise_workflow.connect(t_node, "outputspec.max_file", dataSink, "max_file")
    randomise_workflow.connect(t_node, "outputspec.mean_file", dataSink, "mean_file")
    randomise_workflow.connect(t_node, "outputspec.pval_file", dataSink, "pval_file")
    randomise_workflow.connect(t_node, "outputspec.size_file", dataSink, "size_file")
    randomise_workflow.connect(
        t_node, "outputspec.tstat_files", dataSink, "tstat_files"
    )
    randomise_workflow.connect(
        t_node, "outputspec.t_corrected_p_files", dataSink, "t_corrected_p_files"
    )
    if run:
        plugin_args = {"n_procs": num_of_cores}
        randomise_workflow.run(
            plugin=MultiProcPlugin(plugin_args), plugin_args=plugin_args
        )
        return None
    return randomise_workflow, randomise_workflow.base_dir
