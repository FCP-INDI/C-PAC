from __future__ import print_function, division, unicode_literals, absolute_import

from nipype.interfaces.base import (traits, isdefined, File, InputMultiPath)

from nipype.interfaces.afni.base import (AFNICommand, AFNICommandInputSpec, AFNICommandOutputSpec)


class TprojectInputSpec(AFNICommandInputSpec):

    in_file = File(
        desc='input file to 3dTproject',
        argstr='-input %s',
        mandatory=True,
        exists=True,
        copyfile=False
    )

    out_file = File(
        name_template='%s_tproject_residual',
        desc='output file from 3dTproject',
        argstr='-prefix %s',
        name_source='in_file',
        genfile=True
    )

    censor_file = File(
        desc='censor file, single column with a row for every volume in dataset, 1 indicates the volume should be '
             'retained, 0 indicates it should be censored',
        argstr='-censor %s',
        exists=True
    )

    censor_idx = traits.List(
        traits.Int(),
        desc='List of indices of volumes to be censored, the first volume\'s index = 0',
        argstr='-CENSORTR %s',
        minlen=1
    )

    _censor_modes = [
        "KILL",
        "ZERO",
        "NTRP"
    ]

    censor_mode = traits.Enum(
        *_censor_modes,
        argstr='-cenmode %s',
        desc='Defines the method to use for censoring. \'KILL\' (the default) removes the volume from the time series, '
             '\'ZERO\' replaces censored volumes with zeros, and \'NTRP\' uses interpolation to replace censored '
             'volumes. Censoring is performed before orthogonalization.'
    )

    catenation_file = File(
        desc='If the input dataset was formed by concatenating other files, it must be handled differently to make '
             'sure that the filtering is correct. The catenation file contains the indices of the _start_ of each of '
             'the concatenated files in the input dataset. If the input dataset is automatically catenated from a '
             'collection of datasets, then the run start indexes are determined directly, and this parameter is not '
             'needed',
        argstr='-concat %s',
        exists=True
    )

    noblock = traits.Bool(
        argstr='-noblock',
        desc='As in 3dDeconvolve, if you want the program to treat an auto-catenated dataset as one long run, use '
             'this option. This flag will be ignored if used with the -concat option.'
    )

    orthogonalize_file = InputMultiPath(
        File(exists=True),
        argstr='-ort %s',
        desc='Also orthogonalize input to columns in f.1D. Multiple \'-ort\' '
             'options are allowed.'
    )

    orthogonalize_polynomial = traits.Int(
        argstr='-polort %s',
        desc='Include Legendre polynomials up to and including degree corresponding to this value. Default is 2. It '
             'makes no sense to use a value of pp greater than 2, if you are filtering out the lower frequencies! '
             'For catenated datasets, each run gets a separate set of pp+1 Legendre polynomial regressors.'
    )

    orthogonalize_dataset = File(
        exists=True,
        argstr='-dsort %s',
        desc='Orthogonalize each voxel to the corresponding voxel time series '
             'in dataset \'fset\', which must have the same spatial and '
             'temporal grid structure as the main input dataset. At present, '
             'only one \'-dsort\' option is allowed.'
    )

    bandpass = traits.List(
        traits.Float(),
        desc='Remove all frequencies EXCEPT those in the range fbot..ftop. ++ Only one -bandpass option is allowed.',
        argstr='-bandpass %s',
        minlen=2,
        maxlen=2
    )

    stopband = traits.List(
        traits.Float(),
        argstr='-stopband %s',
        desc='Remove all frequencies in the range sbot..stop. More than one -stopband option is allowed. To specify '
             'an additional add a string containing one or more \'-stopband sbot stop\' sequences to the extra_args '
             'parameter',
        minlen=2,
        maxlen=2
    )

    tr = traits.Float(
        argstr='-dt %f',
        desc='Set time step (TR) in sec [default=from dataset header].'
    )

    mask = File(
        desc='mask file',
        position=2,
        argstr='-mask %s',
        exists=True
    )

    automask = traits.Bool(
        argstr='-automask',
        desc='Create a mask from the input dataset.'
    )

    blur = traits.Float(
        argstr='-blur %f',
        desc='Blur (inside the mask only) with a filter width (FWHM) of '
             '\'fff\' millimeters.'
    )

    normalize = traits.Bool(
        argstr='-norm',
        desc='Make all output time series have L2 norm = 1 (i.e., sum of '
             'squares = 1).'
    )

    verb = traits.Bool(
        argstr='-verb',
        desc='The program will save the fixed ort matrix and its singular values into .1D files, for post-mortems. It '
             'will also print out more progress messages, which might help with figuring out what\'s happening when '
             'problems occur.'
    )

    extra_args = traits.Str(
        argstr='%s',
        desc='Extra arguments to pass to 3dTproject, particularly useful for adding additional stop bands'
    )


class Tproject(AFNICommand):
    """This program projects (regresses) out various 'nuisance' time series from each
     voxel in the input dataset.  Note that all the projections are done via linear
     regression, including the frequency-based options such as '-passband'.  In this
     way, you can bandpass time-censored data, and at the same time, remove other
     time series of no interest (e.g., physiological estimates, motion parameters).

     For complete details, see the `3dTproject Documentation.
     <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTproject.html>`_

     Examples
     ========
     >>> from nipype.interfaces import afni
     >>> from nipype.testing import  example_data
     >>> nuisance_regression = afni.Tproject()
     >>> nuisance_regression.inputs.in_file = 'functional.nii'
     >>> nuisance_regression.inputs.bandpass = [0.005 0.1]
     >>> nuisance_regression.cmdline  # doctest: +ALLOW_UNICODE
     # '3dTproject -prefix functional_bp -input functional.nii -passband 0.005000 0.100000'
     >>> res = nuisance_regression.run()  # doctest: +SKIP
    """

    _cmd = '3dTproject'
    input_spec = TprojectInputSpec
    output_spec = AFNICommandOutputSpec


class LocalstatInputSpec(AFNICommandInputSpec):

    in_file = File(
        desc='input file to 3dTproject',
        argstr='-input %s',
        mandatory=True,
        exists=True,
        copyfile=False,
        position=-1
    )

    neighborhood = traits.Str(
        desc="""   
            String that defines the region around each
            voxel that will be extracted for the statistics
            calculation.  The format of the 'nnn' string are:
            * 'SPHERE(r)' where 'r' is the radius in mm;
              the neighborhood is all voxels whose center-to-
              center distance is less than or equal to 'r'.
              ** The distances are computed in 3 dimensions,
                 so a SPHERE(1) on a 1mm3 grid gives a 7 voxel-
                 neighborhood - the center voxel and the six
                 facing voxels, 4 in plane and 2 above and below.
                 A SPHERE(1.42) contains 19 voxels, the center voxel
                 with the 8 others in plane, and the 5 above and
                 below (all voxels sharing an edge with the center)
                 A SPHERE(1.74) contains 27 voxels, all voxels
                 sharing a face, edge or corner with the center
              ** A negative value for 'r' means that the region
                 is calculated using voxel indexes rather than
                 voxel dimensions; that is, the neighborhood
                 region is a "sphere" in voxel indexes of
                 "radius" abs(r).
            * 'RECT(a,b,c)' is a rectangular block which
              proceeds plus-or-minus 'a' mm in the x-direction,
              'b' mm in the y-direction, and 'c' mm in the
              z-direction.  The correspondence between the
              dataset xyz axes and the actual spatial orientation
              can be determined by using program 3dinfo.
              ** Note the a,b,c are not the full dimensions of
                 of the block. They are radially used - effectively
                 half the dimension of a side. So if one wanted to
                 compute a 5-slice projection on a 1mm3 volume,
                 then a RECT(0,0,2) would be appropriate, and 
                 the program would report 5 voxels used in the mask
                 Any dimension less than a voxel will avoid
                 voxels in that direction.
              ** A negative value for 'a' means that the region
                 extends plus-and-minus abs(a) voxels in the
                 x-direction, rather than plus-and-minus a mm.
                 Mutatis mutandum for negative 'b' and/or 'c'.
            * 'RHDD(a)' where 'a' is the size parameter in mm;
              this is Kepler's rhombic dodecahedron [volume=2*a^3].
            * 'TOHD(a)' where 'a' is the size parameter in mm;
              this is a truncated octahedron. [volume=4*a^3]
              ** This is the polyhedral shape that tiles space
                 and is the most 'sphere-like'.
            * If no '-nbhd' option is given, the region extracted
              will just be the voxel and its 6 nearest neighbors.
            * Voxels not in the mask (if any) or outside the
              dataset volume will not be used.  This means that
              different output voxels will have different numbers
              of input voxels that went into calculating their
              statistics.  The 'num' statistic can be used to
              get this count on a per-voxel basis, if you need it.
        """,
        argstr='-nbhd %s'
    )

    statistic = traits.Str(
        desc="""   
            Compute the statistic named 'sss' on the values
            extracted from the region around each voxel:
            * mean   = average of the values
            * stdev  = standard deviation
            * var    = variance (stdev*stdev)
            * cvar   = coefficient of variation = stdev/fabs(mean)
            * median = median of the values
            * MAD    = median absolute deviation
            * min    = minimum
            * max    = maximum
            * absmax = maximum of the absolute values
            * mode   = mode
            * nzmode = non-zero mode
            * num    = number of the values in the region:
                       with the use of -mask or -automask,
                       the size of the region around any given
                       voxel will vary; this option lets you
                       map that size.  It may be useful if you
                       plan to compute a t-statistic (say) from
                       the mean and stdev outputs.
            * sum    = sum of the values in the region:
            * FWHM   = compute (like 3dFWHM) image smoothness
                       inside each voxel's neighborhood.  Results
                       are in 3 sub-bricks: FWHMx, FWHMy, and FWHMz.
                       Places where an output is -1 are locations
                       where the FWHM value could not be computed
                       (e.g., outside the mask).
            * FWHMbar= Compute just the average of the 3 FWHM values
                       (normally would NOT do this with FWHM also).
            * perc:P0:P1:Pstep = 
                       Compute percentiles between P0 and P1 with a 
                       step of Pstep.
                       Default P1 is equal to P0 and default P2 = 1
            * rank   = rank of the voxel's intensity
            * frank  = rank / number of voxels in neighborhood
            * P2skew = Pearson's second skewness coefficient
                        3 * (mean - median) / stdev 
            * ALL    = all of the above, in that order 
                      (except for FWHMbar and perc).
            * mMP2s  = Exactly the same output as:
                       -stat median -stat MAD -stat P2skew
                       but it a little faster
            * mmMP2s  = Exactly the same output as:
                    -stat mean -stat median -stat MAD -stat P2skew
            * diffs   = Compute differences between central voxel
                        and all neighbors. Values output are the 
                        average difference, followed by the min and max
                        differences.
            * list    = Just output the voxel values in the neighborhood
                        The order in which the neighbors are listed 
                        depends on the neighborhood selected. Only
                        SPHERE results in a neighborhood list sorted by
                        the distance from the center.
                        Regardless of the neighborhood however, the first
                        value should always be that of the central voxel.
            * hist:MIN:MAX:N[:IGN] = Compute the histogram in the voxel's
                        neighborhood. You must specify the min, max, and 
                        the number of bins in the histogram. You can also
                        ignore values outside the [min max] range by 
                        setting IGN to 1. IGN = 0 by default.
                        The histograms are scaled by the number 
                        of values that went into the histogram.
                        That would be the number of non-masked voxels
                        in the neighborhood if outliers are NOT
                        ignored (default).
                    For histograms of labeled datasets, use 3dLocalHistog

            More than one '-stat' option can be used. This can be accomplished
            by passing the additional arguments to 'extra_args'. 
        """,
        argstr='-stat %s'
    )

    mask = File(
        desc='mask file',
        position=2,
        argstr='-mask %s',
        exists=True
    )

    automask = traits.Bool(
        argstr='-automask',
        desc='Create a mask from the input dataset.'
    )

    use_nonmask = traits.Bool(
        argstr='-use_nonmask',
        desc="Voxels NOT in the mask will not have their local statistics computed.  This option will make it so that "
             "voxels not in the mask WILL have their local statistics computed from all voxels in their neighborhood "
             "that ARE in the mask. * You could use '-use_nonmask' to compute the average local white matter time "
             "series, for example, even at non-WM voxels."
    )

    output_data_type = traits.Str(
        desc="Coerce the output data to be stored as the given type, which may be byte, short, or float. "
             "Default is float",
        argstr='-datum %s'
    )

    label_ext = traits.Str(
        desc="Append '.LABEXT' to each sub-brick label",
        argstr='-datum %s'
    )

    reduce_grid = traits.List(
        traits.Int(),
        desc="<Rx [Ry Rz]>. Compute output on a grid that is reduced by a factor of Rx Ry Rz in the X, Y, and Z "
             "directions of the input data set. This option speeds up computations at the expense of resolution. You "
             "should only use it when the nbhd is quite large with respect to the input's resolution, and the "
             "resultant stats are expected to be smooth. You can either set Rx, or Rx Ry and Rz. If you only specify "
             "Rx the same value is applied to Ry and Rz.",
        argstr='-reduce_grid %s'
    )

    reduce_restore_grid = traits.List(
        traits.Int(),
        desc="<Rx [Ry Rz]>. Like reduce_grid, but also resample output back to input grid.",
        argstr='-reduce_restore_grid %s'
    )

    reduce_max_vox = traits.Int(
        desc="<MAX_VOX> Like -reduce_restore_grid, but automatically set Rx Ry Rz so that the computation grid is at "
             "a resolution of nbhd/MAX_VOX voxels.",
        argstr="-reduce_max_vox %s"
    )

    _resample_modes = ['NN', 'Li', 'Cu', 'Bk']

    grid_resample_mode = traits.Enum(
        *_resample_modes,
        argstr='-grid_rmode %s',
        desc="Interpolation to use when resampling the output with reduce_restore_grid option. The resampling method "
             "string RESAM should come from the set {'NN', 'Li', 'Cu', 'Bk'}.  These stand for 'Nearest Neighbor', "
             "'Linear', 'Cubic' and 'Blocky' interpolation, respectively. Default is Linear"
    )

    proceed_small_N = traits.Bool(
        argstr='-proceed_small_N',
        desc='Do not crash if neighborhood is too small for certain estimates.'
    )

    verb = traits.Bool(
        argstr='-verb',
        desc='The program will save the fixed ort matrix and its singular values into .1D files, for post-mortems. It '
             'will also print out more progress messages, which might help with figuring out what\'s happening when '
             'problems occur.'
    )

    extra_args = traits.Str(
        argstr='%s',
        desc='Extra arguments to pass to 3dLocalstat, particularly useful for adding additional stats'
    )


class Localstat(AFNICommand):
    """This program computes statistics at each voxel, based on a local neighborhood of that voxel.
          - The neighborhood is defined by the '-nbhd' option.
          - Statistics to be calculated are defined by the '-stat' option(s).

     For complete details, see the `3dTproject Documentation.
     <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dLocalstat.html>`_

    """
    _cmd = '3dLocalstat'
    input_spec = LocalstatInputSpec
    output_spec = AFNICommandOutputSpec
