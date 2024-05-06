# Copyright (C) 2012-2023  C-PAC Developers

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
"""Functions to calculate motion statistics"""
import os
import sys
from typing import Optional
import nibabel as nb
from nipype.interfaces.afni.base import (AFNICommand, AFNICommandInputSpec)
from nipype.interfaces.base import (TraitedSpec, traits, File)
from nipype.interfaces import utility as util
import numpy as np
import pandas as pd
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.interfaces.function import Function
from CPAC.utils.pytest import skipif
from CPAC.utils.typing import LITERAL, TUPLE


def motion_power_statistics(name='motion_stats',
                            motion_correct_tool='3dvolreg',
                            filtered=False):
    """
    The main purpose of this workflow is to get various statistical measures
     from the movement/motion parameters obtained in functional preprocessing.

    Parameters
    ----------
    :param str name: Name of the workflow, defaults to 'motion_stats'
    :return: Nuisance workflow.
    :rtype: nipype.pipeline.engine.Workflow

    Notes
    -----

    Workflow Inputs::

        inputspec.motion_correct : string (func/rest file or a list of func/rest nifti file)
            Path to motion corrected functional data

        inputspec.max_displacement : string (Mat file)
            maximum displacement (in mm) vector for brain voxels in each volume.
            This file is obtained in functional preprocessing step

        inputspec.movement_parameters : string (Mat file)
            1D file containing six movement/motion parameters(3 Translation, 3 Rotations)
            in different columns (roll pitch yaw dS  dL  dP), obtained in functional preprocessing step


    Workflow Outputs::

        outputspec.FDP_1D : 1D file
            mean Framewise Displacement (FD)

        outputspec.power_params : txt file
            Text file containing various power parameters for scrubbing

        outputspec.motion_params : txt file
            Text file containing various movement parameters


    Order of commands:

    - Calculate Framewise Displacement FD as per power et al., 2012

      Differentiating head realignment parameters across frames yields a six dimensional timeseries that represents instantaneous head motion.
      Rotational displacements are converted from degrees to millimeters by calculating displacement on the surface of a sphere of radius 50 mm.[R5]

    - Calculate Framewise Displacement FD as per jenkinson et al., 2002

    - Calculate DVARS

      DVARS (D temporal derivative of timecourses, VARS referring to RMS variance over voxels) indexes
      the rate of change of BOLD signal across the entire brain at each frame of data.To calculate
      DVARS, the volumetric timeseries is differentiated (by backwards differences) and RMS signal
      change is calculated over the whole brain.DVARS is thus a measure of how much the intensity
      of a brain image changes in comparison to the previous timepoint (as opposed to the global
      signal, which is the average value of a brain image at a timepoint).[R5]

    - Calculate Power parameters::

        MeanFD : Mean (across time/frames) of the absolute values for Framewise Displacement (FD),
        computed as described in Power et al., Neuroimage, 2012)

        rootMeanSquareFD : Root mean square (RMS; across time/frames) of the absolute values for FD

        rmsFD : Root mean square (RMS; across time/frames) of the absolute values for FD

        FDquartile(top 1/4th FD) : Mean of the top 25% highest FD values

        MeanDVARS : Mean of voxel DVARS

    - Calculate Motion Parameters

      Following motion parameters are calculated::

        Scan
        Mean Relative RMS Displacement
        Max Relative RMS Displacement
        Movements > threshold
        Mean Relative Mean Rotation
        Mean Relative Maxdisp
        Max Relative Maxdisp
        Max Abs Maxdisp
        Max Relative Roll
        Max Relative Pitch
        Max Relative Yaw
        Max Relative dS-I
        Max Relative dL-R
        Max Relative dP-A
        Mean Relative Roll
        Mean Relative Pitch
        Mean Relative Yaw
        Mean Relative dS-I
        Mean Relative dL-R
        Mean Relative dP-A
        Max Abs Roll
        Max Abs Pitch
        Max Abs Yaw
        Max Abs dS-I
        Max Abs dL-R
        Max Abs dP-A
        Mean Abs Roll
        Mean Abs Pitch
        Mean Abs Yaw
        Mean Abs dS-I
        Mean Abs dL-R
        Mean Abs dP-A

    .. exec::
        from CPAC.generate_motion_statistics import motion_power_statistics
        wf = motion_power_statistics()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/motion_statistics.dot'
        )

    High Level Workflow Graph:

    .. image:: ../../images/generated/motion_statistics.png
       :width: 1000

    Detailed Workflow Graph:

    .. image:: ../../images/generated/motion_statistics_detailed.png
       :width: 1000

    Examples
    --------
    >>> from CPAC import generate_motion_statistics
    >>> wf = generate_motion_statistics.motion_power_statistics("generate_statistics")
    >>> wf.inputs.inputspec.movement_parameters = 'CPAC_outupts/sub01/func/movement_parameteres/rest_mc.1D'  # doctest: +SKIP
    >>> wf.inputs.inputspec.max_displacement = 'CPAC_outputs/sub01/func/max_dispalcement/max_disp.1D'  # doctest: +SKIP
    >>> wf.inputs.inputspec.motion_correct = 'CPAC_outputs/sub01/func/motion_correct/rest_mc.nii.gz'  # doctest: +SKIP
    >>> wf.inputs.inputspec.mask = 'CPAC_outputs/sub01/func/func_mask/rest_mask.nii.gz'  # doctest: +SKIP
    >>> wf.inputs.inputspec.transformations = 'CPAC_outputs/sub01/func/coordinate_transformation/rest_mc.aff12.1D'  # doctest: +SKIP
    >>> wf.base_dir = './working_dir'  # doctest: +SKIP
    >>> wf.run()  # doctest: +SKIP

    References
    ----------

    .. [1] Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious
           but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3),
           2142-2154. doi:10.1016/j.neuroimage.2011.10.018

    .. [2] Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Steps
           toward optimizing motion artifact removal in functional connectivity MRI; a reply to Carp.
           NeuroImage. doi:10.1016/j.neuroimage.2012.03.017

    .. [3] Jenkinson, M., Bannister, P., Brady, M., Smith, S., 2002. Improved optimization for the robust
           and accurate linear registration and motion correction of brain images. Neuroimage 17, 825-841.
    """
    wf = pe.Workflow(name=name)
    input_node = pe.Node(util.IdentityInterface(fields=['movement_parameters',
                                                        'max_displacement',
                                                        'rels_displacement',
                                                        'motion_correct',
                                                        'mask',
                                                        'transformations']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['FDP_1D',
                                                         'FDJ_1D',
                                                         'DVARS_1D',
                                                         'power_params',
                                                         'motion_params',
                                                         'motion',
                                                         'desc-summary_motion']
                                                ),
                          name='outputspec')

    cal_DVARS = pe.Node(ImageTo1D(method='dvars'),
                        name='cal_DVARS',
                        mem_gb=0.4,
                        mem_x=(739971956005215 / 151115727451828646838272,
                               'in_file'),
                        throttle=True)

    cal_DVARS_strip = pe.Node(Function(input_names=['file_1D'],
                                       output_names=['out_file', 'DVARS_val'],
                                       function=DVARS_strip_t0,
                                       as_module=True),
                              name='cal_DVARS_strip')

    # calculate mean DVARS
    wf.connect(input_node, 'motion_correct', cal_DVARS, 'in_file')
    wf.connect(input_node, 'mask', cal_DVARS, 'mask')
    wf.connect(cal_DVARS, 'out_file', cal_DVARS_strip, 'file_1D')
    wf.connect(cal_DVARS_strip, 'out_file', output_node, 'DVARS_1D')

    # Calculating mean Framewise Displacement as per power et al., 2012
    calculate_FDP = pe.Node(Function(input_names=['in_file'],
                                     output_names=['out_file', 'fd'],
                                     function=calculate_FD_P,
                                     as_module=True),
                            name='calculate_FD')

    wf.connect(input_node, 'movement_parameters', calculate_FDP, 'in_file')
    wf.connect(calculate_FDP, 'out_file', output_node, 'FDP_1D')

    # Calculating mean Framewise Displacement as per jenkinson et al., 2002
    calculate_FDJ = pe.Node(Function(input_names=['in_file',
                                                  'calc_from',
                                                  'center'],
                                     output_names=['out_file', 'fd'],
                                     function=calculate_FD_J,
                                     as_module=True),
                            name='calculate_FDJ')

    if filtered or motion_correct_tool == '3dvolreg':
        wf.connect(input_node, 'transformations', calculate_FDJ, 'in_file')
        calculate_FDJ.inputs.calc_from = 'affine'
    elif motion_correct_tool == 'mcflirt':
        calculate_FDJ.inputs.calc_from = 'rms'
        wf.connect(input_node, 'rels_displacement', calculate_FDJ, 'in_file')

    wf.connect(calculate_FDJ, 'out_file', output_node, 'FDJ_1D')

    calc_motion_parameters = pe.Node(Function(input_names=[
                                                  'movement_parameters',
                                                  'max_displacement',
                                                  'motion_correct_tool',
                                                  'rels_displacement'],
                                              output_names=['out_file',
                                                            'info',
                                                            'maxdisp',
                                                            'relsdisp'],
                                              function=gen_motion_parameters,
                                              as_module=True),
                                     name='calc_motion_parameters')

    get_all_motion_parameters = pe.Node(Function(input_names=[
                                                     'fdj',
                                                     'fdp',
                                                     'maxdisp',
                                                     'motion',
                                                     'power',
                                                     'relsdisp',
                                                     'dvars'],
                                                 output_names=[
                                                     'all_motion_val',
                                                     'summary_motion_power'],
                                                 function=get_allmotion,
                                                 as_module=True),
                                        name='get_all_motion_parameters')

    calc_motion_parameters.inputs.motion_correct_tool = motion_correct_tool
    wf.connect(calculate_FDJ, 'fd', get_all_motion_parameters, 'fdj')
    wf.connect(calculate_FDP, 'fd', get_all_motion_parameters, 'fdp')
    wf.connect(calc_motion_parameters, 'maxdisp',
               get_all_motion_parameters, 'maxdisp')
    wf.connect(calc_motion_parameters, 'relsdisp',
               get_all_motion_parameters, 'relsdisp')
    wf.connect(calc_motion_parameters, 'info',
               get_all_motion_parameters, 'motion')
    wf.connect(cal_DVARS_strip, 'DVARS_val',
               get_all_motion_parameters, 'dvars')
    wf.connect(input_node, 'movement_parameters',
               calc_motion_parameters, 'movement_parameters')
    wf.connect(input_node, 'max_displacement',
               calc_motion_parameters, 'max_displacement')
    wf.connect(input_node, 'rels_displacement',
               calc_motion_parameters, 'rels_displacement')
    wf.connect(calc_motion_parameters, 'out_file',
               output_node, 'motion_params')
    wf.connect(get_all_motion_parameters, 'all_motion_val',
               output_node, 'motion')
    wf.connect(get_all_motion_parameters, 'summary_motion_power',
               output_node, 'desc-summary_motion')

    calc_power_parameters = pe.Node(Function(input_names=['fdp',
                                                          'fdj',
                                                          'dvars',
                                                          'motion_correct_tool'],
                                             output_names=['out_file', 'info'],
                                             function=gen_power_parameters,
                                             as_module=True),
                                    name='calc_power_parameters')

    calc_power_parameters.inputs.motion_correct_tool = motion_correct_tool
    wf.connect(calc_power_parameters, 'info',
               get_all_motion_parameters, 'power')
    wf.connect(cal_DVARS, 'out_file',
               calc_power_parameters, 'dvars')

    wf.connect(calculate_FDP, 'out_file',
               calc_power_parameters, 'fdp')

    if motion_correct_tool == '3dvolreg':
        wf.connect(calculate_FDJ, 'out_file', calc_power_parameters, 'fdj')

    wf.connect(calc_power_parameters, 'out_file',
               output_node, 'power_params')

    return wf


def calculate_FD_P(in_file):
    """
    Method to calculate Framewise Displacement (FD)  as per Power et al., 2012

    Parameters
    ----------
    in_file : string
        movement parameters vector file path

    Returns
    -------
    out_file : string
        Frame-wise displacement mat
        file path
    fd : array
        Frame-wise displacement mat
    """

    motion_params = np.genfromtxt(in_file).T

    rotations = np.transpose(np.abs(np.diff(motion_params[0:3, :])))
    translations = np.transpose(np.abs(np.diff(motion_params[3:6, :])))

    fd = np.sum(translations, axis=1) + \
         (50 * np.pi / 180) * np.sum(rotations, axis=1)

    fd = np.insert(fd, 0, 0)

    out_file = os.path.join(os.getcwd(), 'FD.1D')
    np.savetxt(out_file, fd)

    return out_file, fd


@Function.sig_imports(['import os', 'import sys',
                       'from typing import Optional',
                       'import numpy as np',
                       'from CPAC.utils.pytest import skipif',
                       'from CPAC.utils.typing import LITERAL, TUPLE'])
@skipif(sys.version_info < (3, 10),
        reason="Test requires Python 3.10 or higher")
def calculate_FD_J(in_file: str, calc_from: LITERAL['affine', 'rms'],
                   center: Optional[np.ndarray] = None
                  ) -> TUPLE[str, np.ndarray]:
    """
    Method to calculate framewise displacement as per Jenkinson et al. 2002

    Parameters
    ----------
    in_file : string
        matrix transformations from volume alignment file path if
        calc_from is 'affine', or FDRMS (*_rel.rms) output if
        calc_from is 'rms'.
    calc_from : string
        one of {'affine', 'rms'}
    center : ~numpy.ndarray, optional
        optional volume center for the from-affine calculation

    Returns
    -------
    out_file : string
        Framewise displacement file path

    fdj : ~numpy.ndarray
        Framewise displacement array

    Examples
    --------
    The file and array output by this function and the "rels_rms"
    property of the pickled test data (offset by a leading zero)
    should all be equal (rounded to the neareast 0.001):
    >>> import gzip, os, pickle
    >>> from unittest import mock
    >>> import numpy as np
    >>> with gzip.open('/code/CPAC/generate_motion_statistics/test/'
    ...                'fdj_test_data.pklz') as _pickle:
    ...     test_data = pickle.load(_pickle)
    >>> with mock.patch('nibabel.load',
    ...                 return_value=test_data.img), mock.patch(
    ...        'numpy.genfromtxt', return_value=test_data.affine):
    ...     fdj_file, fdj = calculate_FD_J(
    ...         test_data.affine, calc_from='affine',
    ...         center=find_volume_center(test_data.img))
    >>> fdj_from_file = np.genfromtxt(fdj_file)
    >>> fdj_test_data = np.insert(test_data.rels_rms, 0, 0)
    >>> all(np.isclose(fdj, fdj_from_file, atol=0.001))
    True
    >>> all(np.isclose(fdj, fdj_test_data, atol=0.001))
    True
    >>> all(np.isclose(fdj_from_file, fdj_test_data, atol=0.001))
    True
    >>> os.unlink(fdj_file)
    """
    if calc_from == 'affine':
        if center is None:
            center = np.zeros((3, 1))
        else:
            center = np.asarray(center).reshape((3, 1))
        pm_ = np.genfromtxt(in_file)

        pm = np.zeros((pm_.shape[0], pm_.shape[1] + 4))
        pm[:, :12] = pm_
        pm[:, 12:] = [0.0, 0.0, 0.0, 1.0]

        # The default radius (as in FSL) of a sphere represents the brain
        rmax = 80.0

        T_rb_prev = pm[0].reshape(4, 4)

        fd = np.zeros(pm.shape[0])

        for i in range(1, pm.shape[0]):
            T_rb = pm[i].reshape(4, 4)

            M = np.dot(T_rb, np.linalg.inv(T_rb_prev)) - np.eye(4)
            A = M[0:3, 0:3]
            b = M[0:3, 3:4] + A @ center

            fd[i] = np.sqrt(
                (rmax * rmax / 5) * np.trace(np.dot(A.T, A)) + np.dot(b.T, b)
            )

            T_rb_prev = T_rb

    elif calc_from == 'rms':
        rel_rms = np.loadtxt(in_file)
        fd = np.append(0, rel_rms)

    else:
        raise ValueError(f"calc_from {calc_from} not supported")

    out_file = os.path.join(os.getcwd(), 'FD_J.1D')
    np.savetxt(out_file, fd, fmt='%.8f')

    return out_file, fd


def find_volume_center(img_file : str) -> np.ndarray:
    """
    Find the center of mass of a Nifti image volume

    Parameters
    ----------
    img_file : string (nifti file)
        path to nifti volume image

    Returns
    -------
    center : ndarray 
        volume center of mass vector
    """
    img = nb.load(img_file)
    dim = np.array(img.header["dim"][1:4])
    pixdim = np.array(img.header["pixdim"][1:4])
    # Calculation follows MCFLIRT
    # https://github.com/fithisux/FSL/blob/7aa2932949129f5c61af912ea677d4dbda843895/src/mcflirt/mcflirt.cc#L479
    center = 0.5 * (dim - 1) * pixdim
    return center


def gen_motion_parameters(movement_parameters, max_displacement,
                          motion_correct_tool, rels_displacement=None):
    """
    Method to calculate all the movement parameters

    Parameters
    ----------
    max_displacement : string
        path of file with maximum displacement (in mm) for brain voxels
        in each volume

    movement_parameters : string
        path of 1D file containing six movement/motion parameters
        (3 Translation, 3 Rotations) in different columns
        (roll pitch yaw dS  dL  dP)

    Returns
    -------
    out_file : string
        path to csv file containing various motion parameters

    info : text
        contains information about motion parameters

    maxdisp : array
        max displacement value

    relsdisp : array
        rels displacement value
    """
    mot = np.genfromtxt(movement_parameters).T

    # Relative RMS of translation
    rms = np.sqrt(mot[3] ** 2 + mot[4] ** 2 + mot[5] ** 2)

    # remove any other information other than matrix from
    # max displacement file. AFNI adds information to the file
    if motion_correct_tool == '3dvolreg':
        maxdisp = np.loadtxt(max_displacement)
        relsdisp = []
        relsdisp = pd.DataFrame(relsdisp)

    elif motion_correct_tool == 'mcflirt':
        # TODO: mcflirt outputs absdisp, instead of maxdisp
        maxdisp = np.loadtxt(max_displacement)

        # rels_disp output only for mcflirt
        relsdisp = np.loadtxt(rels_displacement)

    abs_relative = lambda v: np.abs(np.diff(v))
    max_relative = lambda v: np.max(abs_relative(v))
    avg_relative = lambda v: np.mean(abs_relative(v))
    max_abs = lambda v: np.max(np.abs(v))
    avg_abs = lambda v: np.mean(np.abs(v))

    info = [
        ('Mean_Relative_RMS_Displacement', avg_relative(rms)),
        ('Max_Relative_RMS_Displacement', max_relative(rms)),
        ('Movements_gt_threshold', np.sum(abs_relative(rms) > 0.1)),
        ('Mean_Relative_Mean_Rotation',
         avg_relative(np.abs(mot[0:3]).mean(axis=0))),
        ('Mean_Relative_Maxdisp', avg_relative(maxdisp)),  # to be updated
        ('Max_Relative_Maxdisp', max_relative(maxdisp)),  # to be updated
        ('Max_Abs_Maxdisp', max_abs(maxdisp)),  # to be updated
        ('Max Relative_Roll', max_relative(mot[0])),
        ('Max_Relative_Pitch', max_relative(mot[1])),
        ('Max_Relative_Yaw', max_relative(mot[2])),
        ('Max_Relative_dS-I', max_relative(mot[3])),
        ('Max_Relative_dL-R', max_relative(mot[4])),
        ('Max_Relative_dP-A', max_relative(mot[5])),
        ('Mean_Relative_Roll', avg_relative(mot[0])),
        ('Mean_Relative_Pitch', avg_relative(mot[1])),
        ('Mean_Relative_Yaw', avg_relative(mot[2])),
        ('Mean_Relative_dS-I', avg_relative(mot[3])),
        ('Mean_Relative_dL-R', avg_relative(mot[4])),
        ('Mean_Relative_dP-A', avg_relative(mot[5])),
        ('Max_Abs_Roll', max_abs(mot[0])),
        ('Max_Abs_Pitch', max_abs(mot[1])),
        ('Max_Abs_Yaw', max_abs(mot[2])),
        ('Max_Abs_dS-I', max_abs(mot[3])),
        ('Max_Abs_dL-R', max_abs(mot[4])),
        ('Max_Abs_dP-A', max_abs(mot[5])),
        ('Mean_Abs_Roll', avg_abs(mot[0])),
        ('Mean_Abs_Pitch', avg_abs(mot[1])),
        ('Mean_Abs_Yaw', avg_abs(mot[2])),
        ('Mean_Abs_dS-I', avg_abs(mot[3])),
        ('Mean_Abs_dL-R', avg_abs(mot[4])),
        ('Mean_Abs_dP-A', avg_abs(mot[5])),
    ]

    out_file = os.path.join(os.getcwd(), 'motion_parameters.txt')
    with open(out_file, 'w') as f:
        f.write(','.join(t for t, v in info))
        f.write('\n')
        f.write(','.join(
            v if type(v) == str else '{0:.6f}'.format(v) for t, v in info))
        f.write('\n')

    return out_file, info, maxdisp, relsdisp


def gen_power_parameters(fdp=None, fdj=None, dvars=None,
                         motion_correct_tool='3dvolreg'):
    """
    Method to generate Power parameters for scrubbing

    Parameters
    ----------
    fdp : string
        framewise displacement(FD as per power et al., 2012) file path
    fdj : string
        framewise displacement(FD as per jenkinson et al., 2002) file path
    dvars : string
        path to numpy file containing DVARS

    Returns
    -------
    out_file : string (csv file)
        path to csv file containing all the pow parameters
    info : text
        contains information about power parameters
    """
    import numpy as np

    meanFD_Power = []
    meanDVARS = []
    meanFD_Jenkinson = []
    rmsFDJ = []
    FDJquartile = []

    if fdp:
        fdp_data = np.loadtxt(fdp)
        dvars_data = np.loadtxt(dvars)

        # Mean (across time/frames) of the absolute values
        # for Framewise Displacement (FD)
        meanFD_Power = np.mean(fdp_data)

        # Mean DVARS
        meanDVARS = np.mean(dvars_data)

        if motion_correct_tool == '3dvolreg':
            if fdj:
                fdj_data = np.loadtxt(fdj)

                # Mean FD Jenkinson
                meanFD_Jenkinson = np.mean(fdj_data)

                # Root mean square (RMS; across time/frames)
                # of the absolute values for FD
                rmsFDJ = np.sqrt(np.mean(fdj_data))

                # Mean of the top quartile of FD is $FDquartile
                quat = int(len(fdj_data) / 4)
                FDJquartile = np.mean(np.sort(fdj_data)[::-1][:quat])

                info = [
                    ('MeanFD_Power', meanFD_Power),
                    ('MeanFD_Jenkinson', meanFD_Jenkinson),
                    ('rootMeanSquareFD', rmsFDJ),
                    ('FDquartile(top1/4thFD)', FDJquartile),
                    ('MeanDVARS', meanDVARS)]

        elif motion_correct_tool == 'mcflirt':
            info = [
                ('MeanFD_Power', meanFD_Power),
                ('MeanDVARS', meanDVARS)]

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    with open(out_file, 'a') as f:
        f.write(','.join(t for t, v in info))
        f.write('\n')
        f.write(','.join(
            v if type(v) == str else '{0:.4f}'.format(v) for t, v in info))
        f.write('\n')

    return out_file, info


def DVARS_strip_t0(file_1D):
    x = np.loadtxt(file_1D)
    x = x[1:]
    x = np.insert(x, 0, 0)
    np.savetxt('dvars_strip.1D', x)
    return os.path.abspath('dvars_strip.1D'), x


class ImageTo1DInputSpec(AFNICommandInputSpec):
    in_file = File(desc='input file to 3dTto1D',
                   argstr='-input %s',
                   position=1,
                   mandatory=True,
                   exists=True,
                   copyfile=False)

    mask = File(desc='-mask dset = use dset as mask to include/exclude voxels',
                argstr='-mask %s',
                position=2,
                exists=True)

    out_file = File(name_template="%s_3DtoT1.1D", desc='output 1D file name',
                    argstr='-prefix %s', name_source="in_file", keep_extension=True)

    _methods = [
        'enorm', 'dvars',
        'rms', 'srms', 's_srms',
        'mdiff', 'smdiff',
        '4095_count', '4095_frac', '4095_warn',
    ]

    method = traits.Enum(
        *_methods,
        argstr='-method %s'
    )


class ImageTo1DOutputSpec(TraitedSpec):
    out_file = File(desc='output 1D file name')


class ImageTo1D(AFNICommand):
    _cmd = '3dTto1D'
    input_spec = ImageTo1DInputSpec
    output_spec = ImageTo1DOutputSpec


def calculate_DVARS(func_brain, mask):
    """
    Method to calculate DVARS as per power's method

    Parameters
    ----------
    func_brain : string (nifti file)
        path to motion correct functional data
    mask : string (nifti file)
        path to brain only mask for functional data

    Returns
    -------
    out_file : string (numpy mat file)
        path to file containing array of DVARS calculation for each voxel
    dvars: array
        file containing array of DVARS calculation for each voxel
    """
    import numpy as np
    import nibabel as nb
    rest_data = nb.load(func_brain).get_fdata().astype(np.float32)
    mask_data = nb.load(mask).get_fdata().astype('bool')

    # square of relative intensity value for each voxel across every timepoint
    data = np.square(np.diff(rest_data, axis=3))

    # applying mask, getting the data in the brain only
    data = data[mask_data]

    # square root and mean across all timepoints inside mask
    dvars = np.sqrt(np.mean(data, axis=0))

    out_file = os.path.join(os.getcwd(), 'DVARS.txt')
    np.savetxt(out_file, dvars)

    dvars = np.insert(dvars, 0, 0)

    return out_file, dvars


def get_allmotion(fdj, fdp, maxdisp, motion, power, relsdisp=None, dvars=None):
    """
    Method to append all the motion and power parameters into 2 files

    Parameters
    ----------
    fdj
        framewise displacement (Jenkinson)
    fdp
        framewise displacement (Power)
    maxdisp
        maximum displacement value
    motion
        motion info
    power
        Power values
    relsdisp
        rels displacement value
    dvars
        DVARS value

    Returns
    -------
    all_motion_val : str
        path to file containing all motion parameters appended

    summary_motion_power : str
        path to file containing all motion parameters appended
    """
    all_motion_val = os.path.join(os.getcwd(), 'motion.tsv')
    summary_motion_power = os.path.join(os.getcwd(), 'desc-summary_motion.tsv')

    df_fdj = pd.DataFrame(fdj)
    df_fdj.columns = ['Framewise displacement Jenkinson']
    df_fdp = pd.DataFrame(fdp)
    df_fdp.columns = ['Framewise displacement Power']
    df_dvars = pd.DataFrame(dvars)
    df_dvars.columns = ['DVARS']
    df_maxdisp = pd.DataFrame(maxdisp)
    df_maxdisp.columns = ['Max Displacement']
    df_relsdisp = pd.DataFrame(relsdisp)
    df_maxdisp.columns = ['Rels Displacement']
    data_frames = [df_fdj, df_fdp, df_dvars, df_maxdisp, df_relsdisp]
    all_motion_val_df = pd.concat(data_frames, axis=1)

    if len(all_motion_val_df.columns) == 5:
        np.savetxt(all_motion_val, all_motion_val_df, delimiter="\t",
                   header="Framewise displacement Jenkinson\tFramewise "
                          "displacement power\tDVARS\tMax Displacement"
                          "\tRels Displacement", comments='')
    if len(all_motion_val_df.columns) == 4:
        np.savetxt(all_motion_val, all_motion_val_df, delimiter="\t",
                   header="Framewise displacement Jenkinson\tFramewise "
                          "displacement power\tDVARS\tMax Displacement",
                   comments='')

    df_motion = pd.DataFrame(motion)
    df_power = pd.DataFrame(power)
    data_frames_motionpower = [df_motion, df_power]
    summary_motion_pow_df = pd.concat(data_frames_motionpower).T
    summary_motion_pow_df.columns = summary_motion_pow_df.iloc[0]
    summary_motion_pow_df.drop(summary_motion_pow_df.index[0], inplace=True)
    summary_motion_pow_df.to_csv(summary_motion_power, sep='\t', header=True,
                                 index=False)

    return all_motion_val, summary_motion_power
