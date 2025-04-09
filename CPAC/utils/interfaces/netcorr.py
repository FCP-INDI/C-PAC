"""Interface for AFNI 3dNetCorr."""

import subprocess

from traits.api import Bool
from nipype.interfaces.afni.preprocess import NetCorr as NipypeNetCorr, NetCorrInputSpec

NetCorrInputSpec.add_class_trait(
    "automask_off",
    Bool(
        False,
        desc="If you want to neither put in a mask *nor* have the automasking occur",
        argstr="-automask_off",
        usedefault=True,
    ),
)


class NetCorr(NipypeNetCorr):  # noqa: D101
    input_spec = NetCorrInputSpec

    def _list_outputs(self):
        """``nipype.interfaces.afni.preprocess.NetCorr._list_outputs`` with a bugfix.

        Notes
        -----
        This method can be removed once nipy/nipype#3697 is merged and a release
        including that PR is included in the C-PAC image.
        """
        # STATEMENT OF CHANGES:
        #     This function is derived from sources licensed under the Apache-2.0 terms,
        #     and this function has been changed.

        # CHANGES:
        #     * Includes changes from https://github.com/nipy/nipype/pull/3697 prior to all commits between https://github.com/nipy/nipype/tree/1.8.6 and that PR being merged and released.

        # ORIGINAL WORK'S ATTRIBUTION NOTICE:
        #    Copyright (c) 2009-2016, Nipype developers

        #    Licensed under the Apache License, Version 2.0 (the "License");
        #    you may not use this file except in compliance with the License.
        #    You may obtain a copy of the License at

        #        http://www.apache.org/licenses/LICENSE-2.0

        #    Unless required by applicable law or agreed to in writing, software
        #    distributed under the License is distributed on an "AS IS" BASIS,
        #    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
        #    See the License for the specific language governing permissions and
        #    limitations under the License.

        #    Prior to release 0.12, Nipype was licensed under a BSD license.

        # Modifications copyright (C) 2024  C-PAC Developers
        import glob
        import os

        from nipype.interfaces.base.traits_extension import isdefined

        outputs = self.output_spec().get()

        if not isdefined(self.inputs.out_file):
            prefix = self._gen_fname(self.inputs.in_file, suffix="_netcorr")
        else:
            prefix = self.inputs.out_file

        # All outputs should be in the same directory as the prefix
        odir = os.path.dirname(os.path.abspath(prefix))
        outputs["out_corr_matrix"] = glob.glob(os.path.join(odir, "*.netcc"))[0]

        if self.inputs.ts_wb_corr or self.inputs.ts_wb_Z:
            corrdir = os.path.join(odir, prefix + "_000_INDIV")
            outputs["out_corr_maps"] = glob.glob(os.path.join(corrdir, "*.nii.gz"))

        return outputs


NetCorr.__doc__ = f"""{NipypeNetCorr.__doc__}
`CPAC.utils.interfaces.netcorr.NetCorr` adds an additional optional input, `automask_off`

Examples
--------
>>> ncorr.inputs.automask_off = True
>>> ncorr.cmdline
'3dNetCorr -prefix sub0.tp1.ncorr -automask_off -fish_z -inset functional.nii -in_rois maps.nii -mask mask.nii -ts_wb_Z -ts_wb_corr'
"""


def strip_afni_output_header(in_file, out_file):
    """Rewrite a file with all but the first 6 lines."""
    subprocess.run(f"tail -n +7 {in_file} > {out_file}", shell=True, check=True)
    return out_file


__all__ = ["NetCorr", "strip_afni_output_header"]
