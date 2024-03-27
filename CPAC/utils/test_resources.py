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
from CPAC.utils.monitoring import WFLOGGER


def setup_test_wf(s3_prefix, paths_list, test_name, workdirs_to_keep=None):
    """Set up a basic template Nipype workflow for testing single nodes or
    small sub-workflows.
    """
    import os
    import shutil

    from CPAC.pipeline import nipype_pipeline_engine as pe
    from CPAC.utils.datasource import check_for_s3
    from CPAC.utils.interfaces.datasink import DataSink

    test_dir = os.path.join(os.getcwd(), test_name)
    work_dir = os.path.join(test_dir, "workdir")
    out_dir = os.path.join(test_dir, "output")

    if os.path.exists(out_dir):
        try:
            shutil.rmtree(out_dir)
        except:
            pass

    if os.path.exists(work_dir):
        for dirname in os.listdir(work_dir):
            if workdirs_to_keep:
                for keepdir in workdirs_to_keep:
                    WFLOGGER.info("%s --- %s\n", dirname, keepdir)
                    if keepdir in dirname:
                        continue
            try:
                shutil.rmtree(os.path.join(work_dir, dirname))
            except:
                pass

    local_paths = {}
    for subpath in paths_list:
        s3_path = os.path.join(s3_prefix, subpath)
        local_path = check_for_s3(s3_path, dl_dir=test_dir)
        local_paths[subpath] = local_path

    wf = pe.Workflow(name=test_name)
    wf.base_dir = os.path.join(work_dir)
    wf.config["execution"] = {
        "hash_method": "timestamp",
        "crashdump_dir": os.path.abspath(test_dir),
    }

    ds = pe.Node(DataSink(), name=f"sinker_{test_name}")
    ds.inputs.base_directory = out_dir
    ds.inputs.parameterization = True

    return (wf, ds, local_paths)
