"""Interface for AFNI 3dNetCorr"""
import glob
import os
import subprocess
from nipype.interfaces.afni.base import AFNICommand, AFNICommandInputSpec, \
                                        AFNICommandOutputSpec
from nipype.interfaces.base import File, isdefined
from traits.api import Bool, List


class NetCorrInputSpec(AFNICommandInputSpec):
    in_file = File(
        desc="input time series file (4D data set)",
        exists=True,
        argstr="-inset %s",
        mandatory=True,
    )
    in_rois = File(
        desc="input set of ROIs, each labelled with distinct integers",
        exists=True,
        argstr="-in_rois %s",
        mandatory=True,
    )
    mask = File(
        desc="can include a whole brain mask within which to "
        "calculate correlation. Otherwise, data should be "
        "masked already",
        exists=True,
        argstr="-mask %s",
    )
    automask_off = Bool(
        False,
        desc='If you want to neither put in a mask '
             '*nor* have the automasking occur',
        argstr='-automask_off', usedefault=True
    )
    weight_ts = File(
        desc="input a 1D file WTS of weights that will be applied "
        "multiplicatively to each ROI's average time series. "
        "WTS can be a column- or row-file of values, but it "
        "must have the same length as the input time series "
        "volume. "
        "If the initial average time series was A[n] for "
        "n=0,..,(N-1) time points, then applying a set of "
        "weights W[n] of the same length from WTS would "
        "produce a new time series:  B[n] = A[n] * W[n]",
        exists=True,
        argstr="-weight_ts %s",
    )
    fish_z = Bool(
        desc="switch to also output a matrix of Fisher Z-transform "
        "values for the corr coefs (r): "
        "Z = atanh(r) , "
        "(with Z=4 being output along matrix diagonals where "
        "r=1, as the r-to-Z conversion is ceilinged at "
        "Z = atanh(r=0.999329) = 4, which is still *quite* a "
        "high Pearson-r value",
        argstr="-fish_z",
    )
    part_corr = Bool(
        desc="output the partial correlation matrix", argstr="-part_corr"
    )
    ts_out = Bool(
        desc="switch to output the mean time series of the ROIs that "
        "have been used to generate the correlation matrices. "
        "Output filenames mirror those of the correlation "
        "matrix files, with a '.netts' postfix",
        argstr="-ts_out",
    )
    ts_label = Bool(
        desc="additional switch when using '-ts_out'. Using this "
        "option will insert the integer ROI label at the start "
        "of each line of the *.netts file created. Thus, for "
        "a time series of length N, each line will have N+1 "
        "numbers, where the first is the integer ROI label "
        "and the subsequent N are scientific notation values",
        argstr="-ts_label",
    )
    ts_indiv = Bool(
        desc="switch to create a directory for each network that "
        "contains the average time series for each ROI in "
        "individual files (each file has one line). "
        "The directories are labelled PREFIX_000_INDIV/, "
        "PREFIX_001_INDIV/, etc. (one per network). Within each "
        "directory, the files are labelled ROI_001.netts, "
        "ROI_002.netts, etc., with the numbers given by the "
        "actual ROI integer labels",
        argstr="-ts_indiv",
    )
    ts_wb_corr = Bool(
        desc="switch to create a set of whole brain correlation maps. "
        "Performs whole brain correlation for each "
        "ROI's average time series; this will automatically "
        "create a directory for each network that contains the "
        "set of whole brain correlation maps (Pearson 'r's). "
        "The directories are labelled as above for '-ts_indiv' "
        "Within each directory, the files are labelled "
        "WB_CORR_ROI_001+orig, WB_CORR_ROI_002+orig, etc., with "
        "the numbers given by the actual ROI integer labels",
        argstr="-ts_wb_corr",
    )
    ts_wb_Z = Bool(
        desc="same as above in '-ts_wb_corr', except that the maps "
        "have been Fisher transformed to Z-scores the relation: "
        "Z=atanh(r). "
        "To avoid infinities in the transform, Pearson values "
        "are effectively capped at |r| = 0.999329 (where |Z| = 4.0). "
        "Files are labelled WB_Z_ROI_001+orig, etc",
        argstr="-ts_wb_Z",
    )
    ts_wb_strlabel = Bool(
        desc="by default, '-ts_wb_{corr,Z}' output files are named "
        "using the int number of a given ROI, such as: "
        "WB_Z_ROI_001+orig. "
        "With this option, one can replace the int (such as '001') "
        "with the string label (such as 'L-thalamus') "
        "*if* one has a labeltable attached to the file",
        argstr="-ts_wb_strlabel",
    )
    nifti = Bool(
        desc="output any correlation map files as NIFTI files "
        "(default is BRIK/HEAD). Only useful if using "
        "'-ts_wb_corr' and/or '-ts_wb_Z'",
        argstr="-nifti",
    )
    output_mask_nonnull = Bool(
        desc="internally, this program checks for where there are "
        "nonnull time series, because we don't like those, in "
        "general.  With this flag, the user can output the "
        "determined mask of non-null time series.",
        argstr="-output_mask_nonnull",
    )
    push_thru_many_zeros = Bool(
        desc="by default, this program will grind to a halt and "
        "refuse to calculate if any ROI contains >10 percent "
        "of voxels with null times series (i.e., each point is "
        "0), as of April, 2017.  This is because it seems most "
        "likely that hidden badness is responsible. However, "
        "if the user still wants to carry on the calculation "
        "anyways, then this option will allow one to push on "
        "through.  However, if any ROI *only* has null time "
        "series, then the program will not calculate and the "
        "user will really, really, really need to address their masking",
        argstr="-push_thru_many_zeros",
    )
    ignore_LT = Bool(
        desc="switch to ignore any label table labels in the "
        "'-in_rois' file, if there are any labels attached",
        argstr="-ignore_LT",
    )
    out_file = File(
        desc="output file name part",
        name_template="%s_netcorr",
        argstr="-prefix %s",
        position=1,
        name_source="in_file",
    )


class NetCorrOutputSpec(AFNICommandOutputSpec):
    out_corr_matrix = File(
        desc="output correlation matrix between ROIs written to a text "
             "file with .netcc suffix"
    )
    out_corr_maps = List(
        File(), desc="output correlation maps in Pearson and/or Z-scores"
    )


class NetCorr(AFNICommand):
    """Calculate correlation matrix of a set of ROIs (using mean time series of
    each). Several networks may be analyzed simultaneously, one per brick.

    For complete details, see the `3dNetCorr Documentation
    <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dNetCorr.html>`_.

    Examples
    --------
    >>> from CPAC.utils.interfaces.netcorr import NetCorr
    >>> ncorr = NetCorr()
    >>> ncorr.inputs.in_file = 'functional.nii'  # doctest: +SKIP
    >>> ncorr.inputs.mask = 'mask.nii'  # doctest: +SKIP
    >>> ncorr.inputs.in_rois = 'maps.nii'  # doctest: +SKIP
    >>> ncorr.inputs.ts_wb_corr = True
    >>> ncorr.inputs.ts_wb_Z = True
    >>> ncorr.inputs.fish_z = True
    >>> ncorr.inputs.out_file = 'sub0.tp1.ncorr'
    >>> ncorr.cmdline  # doctest: +SKIP
    '3dNetCorr -prefix sub0.tp1.ncorr -fish_z -inset functional.nii -in_rois maps.nii -mask mask.nii -ts_wb_Z -ts_wb_corr'
    >>> res = ncorr.run()  # doctest: +SKIP

    """  # noqa: E501  # pylint: disable=line-too-long
    _cmd = "3dNetCorr"
    input_spec = NetCorrInputSpec
    output_spec = NetCorrOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()

        if not isdefined(self.inputs.out_file):
            prefix = self._gen_fname(self.inputs.in_file, suffix="_netcorr")
        else:
            prefix = self.inputs.out_file

        # All outputs should be in the same directory as the prefix
        odir = os.path.dirname(os.path.abspath(prefix))
        outputs["out_corr_matrix"] = glob.glob(
            os.path.join(odir, "*.netcc"))[0]

        if (
            isdefined(self.inputs.ts_wb_corr) or
            isdefined(self.inputs.ts_wb_Z)
        ):
            corrdir = os.path.join(odir, prefix + "_000_INDIV")
            outputs["out_corr_maps"] = glob.glob(
                os.path.join(corrdir, "*.nii.gz"))

        return outputs


def strip_afni_output_header(in_file, out_file):
    """Function to rewrite a file with all but the first 6 lines"""
    subprocess.run(f'tail -n +7 {in_file} > {out_file}', shell=True,
                   check=True)
    return out_file
