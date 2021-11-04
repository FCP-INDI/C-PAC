"""Interface for AFNI 3dNetCorr"""
import os
from nipype.interfaces.afni.base import AFNICommand, AFNICommandInputSpec, \
                                        AFNICommandOutputSpec
from nipype.interfaces.base import File, Str
from traits.api import Bool


class NetCorrInputSpec(AFNICommandInputSpec):
    """Input specification for 3dNetCorr. Although 3dNetCorr "will internally
    'automask', based on when non-uniformly-zero time series are", this
    interface requires an ROI mask."""
    out_file = File(name_template='%s_connectome', position=0,
                    desc='output image file name',
                    argstr=f'-prefix $PWD/%s',
                    name_source=['timeseries'])
    measure = Str(desc='if not Pearson', argstr='%s', position=1,
                  mandatory=False) or None
    timeseries = File(desc='time series file (4D data set)', position=2,
                      exists=True, argstr='-inset %s', mandatory=True)
    parcellation = File(desc='a set of ROIs, each labelled with distinct '
                             'integers', position=3, exists=True,
                        argstr='-in_rois %s', mandatory=True)
    automask_off = Bool(True, desc='If you want to neither put in a mask '
                                   '*nor* have the automasking occur',
                        argstr='-automask_off', usedefault=True)


class NetCorr(AFNICommand):
    """3dNetCorr interface. For complete details, see the `3dNetCorr Documentation.
    https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dNetCorr.html>`_
    """
    _cmd = '3dNetCorr'
    input_spec = NetCorrInputSpec
    output_spec = AFNICommandOutputSpec
    # def _gen_filename(self, name):  # this doesn't seem to do anything?
    #     return '_'.join([*self.inputs.timeseries.split('_')[:-1],
    #                       'connectome.netcc'])
