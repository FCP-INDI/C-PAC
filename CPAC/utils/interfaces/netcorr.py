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
