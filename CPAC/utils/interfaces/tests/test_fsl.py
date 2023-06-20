# Copyright (C) 2021-2023  C-PAC Developers

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
# pylint: disable=invalid-name
"""Tests for FSL interface"""
import os
from pathlib import Path
from ..fsl import Merge


def test_HO_ventricles_exists():
    """Make sure we have this required file
    Ref https://github.com/FCP-INDI/C-PAC/pull/1916#issuecomment-1579286459"""
    assert (Path(os.environ["FSLDIR"]) /
            'data/atlases/HarvardOxford/'
            'HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz').exists()


def test_Merge_inputs():
    input_map = dict(
        args=dict(
            argstr="%s",
        ),
        dimension=dict(
            argstr="-%s",
            mandatory=True,
            position=0,
        ),
        environ=dict(
            nohash=True,
            usedefault=True,
        ),
        in_files=dict(
            argstr="%s",
            mandatory=True,
            position=2,
        ),
        merged_file=dict(
            argstr="%s",
            extensions=None,
            hash_files=False,
            name_source="in_files",
            name_template="%s_merged",
            position=1,
        ),
        output_type=dict(),
        tr=dict(
            argstr="%.2f",
            position=-1,
        ),
    )
    inputs = Merge.input_spec()

    for key, metadata in list(input_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(inputs.traits()[key], metakey) == value


def test_Merge_outputs():
    output_map = dict(
        merged_file=dict(
            extensions=None,
        ),
    )
    outputs = Merge.output_spec()

    for key, metadata in list(output_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(outputs.traits()[key], metakey) == value
