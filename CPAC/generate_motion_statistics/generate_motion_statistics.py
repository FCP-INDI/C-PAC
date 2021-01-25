import os
import numpy as np
import nibabel as nb
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from CPAC.utils.interfaces.function import Function


def motion_power_statistics(name='motion_stats',
                            motion_correct_tool='3dvolreg'):
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

        inputspec.subject_id : string
            Subject name or id

        inputspec.scan_id : string
            Functional Scan id or name

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

        Subject
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
    >>> import generate_motion_statistics
    >>> wf = generate_motion_statistics.motion_power_statistics("generate_statistics")
    >>> wf.inputs.inputspec.movement_parameters = 'CPAC_outupts/sub01/func/movement_parameteres/rest_mc.1D'
    >>> wf.inputs.inputspec.max_displacement = 'CPAC_outputs/sub01/func/max_dispalcement/max_disp.1D'
    >>> wf.inputs.inputspec.motion_correct = 'CPAC_outputs/sub01/func/motion_correct/rest_mc.nii.gz'
    >>> wf.inputs.inputspec.mask = 'CPAC_outputs/sub01/func/func_mask/rest_mask.nii.gz'
    >>> wf.inputs.inputspec.transformations = 'CPAC_outputs/sub01/func/coordinate_transformation/rest_mc.aff12.1D'
    >>> wf.inputs.inputspec.subject_id = 'sub01'
    >>> wf.inputs.inputspec.scan_id = 'rest_1'
    >>> wf.base_dir = './working_dir'
    >>> wf.run()

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
    input_node = pe.Node(util.IdentityInterface(fields=['subject_id',
                                                        'scan_id',
                                                        'movement_parameters',
                                                        'max_displacement',
                                                        'motion_correct',
                                                        'mask',
                                                        'transformations']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['FDP_1D',
                                                         'FDJ_1D',
                                                         'DVARS_1D',
                                                         'power_params',
                                                         'motion_params']),
                          name='outputspec')

    cal_DVARS = pe.Node(Function(input_names=['func_brain', 'mask'],
                                 output_names=['out_file'],
                                 function=calculate_DVARS,
                                 as_module=True),
                        name='cal_DVARS')

    # calculate mean DVARS
    wf.connect(input_node, 'motion_correct', cal_DVARS, 'func_brain')
    wf.connect(input_node, 'mask', cal_DVARS, 'mask')
    wf.connect(cal_DVARS, 'out_file', output_node, 'DVARS_1D')

    # Calculating mean Framewise Displacement as per power et al., 2012
    calculate_FDP = pe.Node(Function(input_names=['in_file'],
                                     output_names=['out_file'],
                                     function=calculate_FD_P,
                                     as_module=True),
                            name='calculate_FD')

    wf.connect(input_node, 'movement_parameters', calculate_FDP, 'in_file')
    wf.connect(calculate_FDP, 'out_file', output_node, 'FDP_1D')

    # Calculating mean Framewise Displacement as per jenkinson et al., 2002
    calculate_FDJ = pe.Node(Function(input_names=['in_file',
                                                  'motion_correct_tool'],
                                     output_names=['out_file'],
                                     function=calculate_FD_J,
                                     as_module=True),
                            name='calculate_FDJ')

    calculate_FDJ.inputs.motion_correct_tool = motion_correct_tool
    if motion_correct_tool == '3dvolreg':
        wf.connect(input_node, 'transformations', calculate_FDJ, 'in_file')
    elif motion_correct_tool == 'mcflirt':
        wf.connect(input_node, 'max_displacement', calculate_FDJ, 'in_file')

    wf.connect(calculate_FDJ, 'out_file', output_node, 'FDJ_1D')

    calc_motion_parameters = pe.Node(Function(input_names=['subject_id',
                                                           'scan_id',
                                                           'movement_parameters',
                                                           'max_displacement',
                                                           'motion_correct_tool'],
                                              output_names=['out_file'],
                                              function=gen_motion_parameters,
                                              as_module=True),
                                     name='calc_motion_parameters')

    calc_motion_parameters.inputs.motion_correct_tool = motion_correct_tool
    wf.connect(input_node, 'subject_id',
               calc_motion_parameters, 'subject_id')
    wf.connect(input_node, 'scan_id',
               calc_motion_parameters, 'scan_id')
    wf.connect(input_node, 'movement_parameters',
               calc_motion_parameters, 'movement_parameters')
    wf.connect(input_node, 'max_displacement',
               calc_motion_parameters, 'max_displacement')

    wf.connect(calc_motion_parameters, 'out_file',
               output_node, 'motion_params')

    calc_power_parameters = pe.Node(Function(input_names=['subject_id',
                                                          'scan_id',
                                                          'fdp',
                                                          'fdj',
                                                          'dvars',
                                                          'motion_correct_tool'],
                                             output_names=['out_file'],
                                             function=gen_power_parameters,
                                             as_module=True),
                                    name='calc_power_parameters')

    calc_power_parameters.inputs.motion_correct_tool = motion_correct_tool
    wf.connect(input_node, 'subject_id',
               calc_power_parameters, 'subject_id')
    wf.connect(input_node, 'scan_id',
               calc_power_parameters, 'scan_id')

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

    """

    motion_params = np.genfromtxt(in_file).T

    rotations = np.transpose(np.abs(np.diff(motion_params[0:3, :])))
    translations = np.transpose(np.abs(np.diff(motion_params[3:6, :])))

    fd = np.sum(translations, axis=1) + \
         (50 * np.pi / 180) * np.sum(rotations, axis=1)

    fd = np.insert(fd, 0, 0)

    out_file = os.path.join(os.getcwd(), 'FD.1D')
    np.savetxt(out_file, fd)

    return out_file


def calculate_FD_J(in_file, motion_correct_tool='3dvolreg'):
    """
    Method to calculate framewise displacement as per Jenkinson et al. 2002

    Parameters
    ----------
    in_file : string
        matrix transformations from volume alignment file path

    Returns
    -------
    out_file : string
        Frame-wise displacement file path

    """

    if motion_correct_tool == '3dvolreg':
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
            b = M[0:3, 3]

            fd[i] = np.sqrt(
                (rmax * rmax / 5) * np.trace(np.dot(A.T, A)) + np.dot(b.T, b)
            )

            T_rb_prev = T_rb

    elif motion_correct_tool == 'mcflirt':
        rel_rms = np.loadtxt(in_file[1])
        fd = np.append(0, rel_rms)

    out_file = os.path.join(os.getcwd(), 'FD_J.1D')
    np.savetxt(out_file, fd, fmt='%.8f')

    return out_file


def gen_motion_parameters(subject_id, scan_id, movement_parameters,
                          max_displacement, motion_correct_tool):
    """
    Method to calculate all the movement parameters

    Parameters
    ----------
    subject_id : string
        subject name or id
    scan_id : string
        scan name or id
    max_displacement : string
        path of file with maximum displacement (in mm) for brain voxels in each volume
    movement_parameters : string
        path of 1D file containing six movement/motion parameters(3 Translation,
        3 Rotations) in different columns (roll pitch yaw dS  dL  dP)

    Returns
    -------
    out_file : string
        path to csv file containing various motion parameters

    """

    mot = np.genfromtxt(movement_parameters).T

    # Relative RMS of translation
    rms = np.sqrt(mot[3] ** 2 + mot[4] ** 2 + mot[5] ** 2)

    # remove any other information other than matrix from
    # max displacement file. AFNI adds information to the file
    if motion_correct_tool == '3dvolreg':
        maxdisp = np.loadtxt(max_displacement)

    elif motion_correct_tool == 'mcflirt':
        maxdisp = np.loadtxt(max_displacement[
                                 0])  # TODO: mcflirt outputs absdisp, instead of maxdisp
        reldisp = np.loadtxt(max_displacement[1])

    abs_relative = lambda v: np.abs(np.diff(v))
    max_relative = lambda v: np.max(abs_relative(v))
    avg_relative = lambda v: np.mean(abs_relative(v))
    max_abs = lambda v: np.max(np.abs(v))
    avg_abs = lambda v: np.mean(np.abs(v))

    info = [
        ('Subject', subject_id),
        ('Scan', scan_id),
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

    return out_file


def gen_power_parameters(subject_id, scan_id, fdp=None, fdj=None, dvars=None,
                         motion_correct_tool='3dvolreg'):
    """
    Method to generate Power parameters for scrubbing

    Parameters
    ----------
    subject_id : string
        subject name or id
    scan_id : string
        scan name or id
    FDP_1D: string
        framewise displacement(FD as per power et al., 2012) file path
    FDJ_1D: string
        framewise displacement(FD as per jenkinson et al., 2002) file path
    threshold : float
        scrubbing threshold set in the configuration
        by default the value is set to 1.0
    DVARS : string
        path to numpy file containing DVARS

    Returns
    -------
    out_file : string (csv file)
        path to csv file containing all the pow parameters
    """

    fdp_data = np.loadtxt(fdp)
    dvars_data = np.loadtxt(dvars)

    # Mean (across time/frames) of the absolute values
    # for Framewise Displacement (FD)
    meanFD_Power = np.mean(fdp_data)

    # Mean DVARS
    meanDVARS = np.mean(dvars_data)

    if motion_correct_tool == '3dvolreg':

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
            ('Subject', subject_id),
            ('Scan', scan_id),
            ('MeanFD_Power', meanFD_Power),
            ('MeanFD_Jenkinson', meanFD_Jenkinson),
            ('rootMeanSquareFD', rmsFDJ),
            ('FDquartile(top1/4thFD)', FDJquartile),
            ('MeanDVARS', meanDVARS),
        ]

    elif motion_correct_tool == 'mcflirt':
        info = [
            ('Subject', subject_id),
            ('Scan', scan_id),
            ('MeanFD_Power', meanFD_Power),
            ('MeanDVARS', meanDVARS),
        ]

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')
    with open(out_file, 'w') as f:
        f.write(','.join(t for t, v in info))
        f.write('\n')
        f.write(','.join(
            v if type(v) == str else '{0:.4f}'.format(v) for t, v in info))
        f.write('\n')

    return out_file


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
    """

    rest_data = nb.load(func_brain).get_data().astype(np.float32)
    mask_data = nb.load(mask).get_data().astype('bool')

    # square of relative intensity value for each voxel across every timepoint
    data = np.square(np.diff(rest_data, axis=3))

    # applying mask, getting the data in the brain only
    data = data[mask_data]

    # square root and mean across all timepoints inside mask
    dvars = np.sqrt(np.mean(data, axis=0))

    out_file = os.path.join(os.getcwd(), 'DVARS.txt')
    np.savetxt(out_file, dvars)
    return out_file
