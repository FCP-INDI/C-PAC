import pytest
import shutil
import tempfile
import yaml
from copy import copy
from unittest import mock


def accept_all(object, name, value):
    return value


@pytest.mark.skip(reason='needs refactored')
def test_trimmer():

    from CPAC.utils.trimmer import the_trimmer, is_datasink, expand_workflow, compute_datasink_dirs
    from CPAC.pipeline.cpac_pipeline import build_workflow
    from CPAC.utils.configuration import Configuration

    import os
    import pkg_resources as p
    
    pipe_config = \
        p.resource_filename(
            "CPAC",
            os.path.join(
                "resources",
                "configs",
                "pipeline_config_template.yml"
            )
        )
    
    data_config = \
        p.resource_filename(
            "CPAC",
            os.path.join(
                "resources",
                "configs",
                "data_config_S3-BIDS-ABIDE.yml"
            )
        )

    data_config = yaml.safe_load(open(data_config, 'r'))
    sub_dict = data_config[0]
        
    c = Configuration(yaml.safe_load(open(pipe_config, 'r')))
    temp_dir = tempfile.mkdtemp()
    c.logDirectory = temp_dir
    c.workingDirectory = temp_dir
    c.outputDirectory = temp_dir
    c.crashLogDirectory = temp_dir



    # Disable functional, let only the anatomical workflow
    c_anatomical = copy(c)
    c_anatomical.runFunctional = [0]

    wf, _, _ = build_workflow(sub_dict['subject_id'], sub_dict, c_anatomical)

    # Create fake files to trick THE TRIMMER
    exec_graph = expand_workflow(wf)
    datasinks = [
        n for n in exec_graph.nodes()
        if is_datasink(n)
    ]
    anat_derivatives = {}
    for datasink in datasinks:
        paths = compute_datasink_dirs(exec_graph, datasink)
        anat_derivatives.update(paths)
        for (node, derivative), path in paths.items():
            os.makedirs(path)
            open(os.path.join(path, '%s.txt' % derivative), 'a').close()



    # Enable functional, so the workflow should only run this
    # and enable trimming
    c_functional = copy(c)
    c_functional.runFunctional = [1]

    wf, _, _ = build_workflow(sub_dict['subject_id'], sub_dict, c_functional)
    exec_wf, _ = the_trimmer(wf)
    exec_graph = exec_wf._graph

    datasinks = [
        n for n in exec_graph.nodes()
        if is_datasink(n)
    ]
    func_derivatives = {}
    for datasink in datasinks:
        paths = compute_datasink_dirs(exec_graph, datasink)
        func_derivatives.update(paths)



    # Assert that the functional pipeline remove all the anatomical nodes,
    # as they were already computed
    assert set(func_derivatives.keys()).intersection(set(anat_derivatives.keys())) == set()
