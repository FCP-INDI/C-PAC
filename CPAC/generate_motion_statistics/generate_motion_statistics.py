import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def motion_power_statistics(wf_name = 'gen_motion_stats'):

    """
    The main purpose of this workflow is to get various statistical measures from the 
    movement/motion parameters obtained in functional preprocessing. These parameters
    (FD calculations) are also required to carry out scrubbing.
    
    Parameters
    ----------
    wf_name : workflow object
        Workflow name
    
    Returns 
    -------
    param_wf : workflow object
          Workflow object containing various movement/motion and power parameters estimates.  
    
    Notes
    -----
    
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/generate_parmeters/generate_parmeters.py>`_
    
    Workflow Inputs::
    
        inputspec.subject_id : string 
            Subject name or id
            
        inputspec.scan_id : string
            Functional Scan id or name
                    
        inputspec.motion_correct : string (func/rest file or a list of func/rest nifti file) 
            Path to motion corrected functional data
            
        inputspec.mask : string (nifti file)
            Path to field contianing brain-only mask for the functional data
                
        inputspec.max_displacement : string (Mat file)
            maximum displacement (in mm) vector for brain voxels in each volume.
            This file is obtained in functional preprocessing step
        
        inputspec.movement_parameters : string (Mat file)
            1D file containing six movement/motion parameters(3 Translation, 3 Rotations) 
            in different columns (roll pitch yaw dS  dL  dP), obtained in functional preprocessing step
        
        threshold_input.threshold : string (float)
            scrubbing threshold
        
    Workflow Outputs::
        
        outputspec.FD_1D : 1D file
            mean Framewise Displacement (FD)
            
        outputspec.frames_ex_1D : 1D file
            Number of frames that would be censored ("scrubbed")
            also removing the offending time frames (i.e., those exceeding the threshold), 
            the preceeding frame, and the two subsequent frames
        
        outputspec.frames_in_1D : 1d file
            Number of frames left after removing for scrubbing
        
        outputspec.power_params : txt file
            Text file various power parameters for scrubbing.
        
        outputspec.motion_params : txt file
           Text file containing various movement parameters
        
    
    Order of commands:
    
    - Calculate Frame Wise Displacement FD
    
      Differentiating head realignment parameters across frames yields a six dimensional timeseries that represents instantaneous head motion.   
      Rotational displacements are converted from degrees to millimeters by calculating displacement on the surface of a sphere of radius 50 mm.[R5]
      
    - Calculate Frames to exclude
    
      Remove all frames which are below the threshold
    
    - Calculate Frames to include
    
      Include all the frames which are above the threshold
    
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
        
        NumFD >=threshold : Number of frames (time points) where movement (FD) exceeded threshold
        
        rmsFD : Root mean square (RMS; across time/frames) of the absolute values for FD
        
        FDquartile(top 1/4th FD) : Mean of the top 25% highest FD values
        
        PercentFD( > threshold) : Number of frames (time points) where movement (FD) exceeded threshold 
                                  expressed as a percentage of the total number of frames (time points)
        
        MeanDVARS : Mean of voxel DVARS
            
    - Calculate Motion Parameters
        
      Following motion parameters are calculated::
         
        Subject, Scan, Mean Relative RMS Displacement, Max Relative RMS Displacement,
        Movements >threshold, Mean Relative Mean Rotation, Mean Relative Maxdisp,
        Max Relative Maxdisp, Max Abs Maxdisp, Max Relative Roll,Max Relative Pitch,
        Max Relative Yaw, Max Relative dS-I, Max Relative dL-R,Max Relative dP-A,
        Mean Relative Roll, Mean Relative Pitch,Mean Relative Yaw, Mean Relative dS-I,
        Mean Relative dL-R, Mean Relative dP-A, Max Abs Roll, Max Abs Pitch, Max Abs Yaw,
        Max Abs dS-I, Max Abs dL-R, Max Abs dP-A, Mean Abs Roll,Mean Abs Pitch,Mean Abs Yaw,
        Mean Abs dS-I,Mean Abs dL-R,Mean Abs dP-A

    
    High Level Workflow Graph:
    
    .. image:: ../images/parameters.dot.png
       :width: 1000
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/parameters_detailed.dot.png
       :width: 1000

    Examples
    --------
    
    >>> import generate_motion_statistics
    >>> wf = generate_motion_statistics.motion_power_statistics()
    >>> wf.inputs.inputspec.movement_parameters = 'CPAC_outupts/sub01/func/movement_parameteres/rest_mc.1D'
    >>> wf.inputs.inputspec.max_displacement = 'CPAC_outputs/sub01/func/max_dispalcement/max_disp.1D'
    >>> wf.inputs.inputspec.motion_correct = 'CPAC_outputs/sub01/func/motion_correct/rest_mc.nii.gz'
    >>> wf.inputs.inputspec.mask = 'CPAC_outputs/sub01/func/func_mask/rest_mask.nii.gz'
    >>> wf.inputs.inputspec.subject_id = 'sub01'
    >>> wf.inputs.inputspec.scan_id = 'rest_1'
    >>> wf.inputs.threshold_input.threshold = 0.5
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
    
     
    """
    pm = pe.Workflow(name=wf_name)
    inputNode = pe.Node(util.IdentityInterface(fields=['subject_id',
                                                       'scan_id',
                                                       'movement_parameters',
                                                       'max_displacement',
                                                       'motion_correct',
                                                       'mask'
                                                    ]),
                        name='inputspec')

    inputnode_threshold = pe.Node(util.IdentityInterface(fields=['threshold']),
                             name='threshold_input')

    outputNode = pe.Node(util.IdentityInterface(fields=['FD_1D',
                                                        'frames_ex_1D',
                                                        'frames_in_1D',
                                                        'power_params',
                                                        'motion_params']),
                        name='outputspec')


    cal_DVARS = pe.Node(util.Function(input_names=['rest', 
                                                   'mask'],
                                           output_names=['out_file'],
                                           function=calculate_DVARS),
                             name='cal_DVARS')
    ##calculate mean DVARS
    pm.connect(inputNode, 'motion_correct', cal_DVARS, 'rest')
    pm.connect(inputNode, 'mask', cal_DVARS, 'mask')
    
    ###Calculating mean Framewise Displacement    
    calculate_FD = pe.Node(util.Function(input_names=['in_file'],
                                         output_names=['out_file'],
                                           function=set_FD),
                             name='calculate_FD')
    
    pm.connect(inputNode, 'movement_parameters', 
               calculate_FD, 'in_file' )
    
    pm.connect(calculate_FD, 'out_file', 
               outputNode, 'FD_1D')

    ##calculating frames to exclude and include after scrubbing
    exclude_frames = pe.Node(util.Function(input_names=['in_file', 
                                                        'threshold'],
                                           output_names=['out_file'],
                                           function=set_frames_ex),
                             name='exclude_frames')

    pm.connect(calculate_FD, 'out_file', 
               exclude_frames, 'in_file')
    pm.connect(inputnode_threshold, 'threshold', 
               exclude_frames, 'threshold')
    
    pm.connect(exclude_frames, 'out_file', 
               outputNode, 'frames_ex_1D')
    

    include_frames = pe.Node(util.Function(input_names=['in_file', 
                                                        'threshold', 
                                                        'exclude_list'],
                                           output_names=['out_file'],
                                           function=set_frames_in),
                             name='include_frames')
    pm.connect(calculate_FD, 'out_file', 
               include_frames, 'in_file')
    pm.connect(inputnode_threshold, 'threshold', 
               include_frames, 'threshold')
    pm.connect(exclude_frames, 'out_file', 
               include_frames, 'exclude_list')

    pm.connect(include_frames, 'out_file', 
               outputNode, 'frames_in_1D')

    
    calc_motion_parameters = pe.Node(util.Function(input_names=["subject_id", 
                                                                "scan_id", 
                                                                "movement_parameters",
                                                                "max_displacement"],
                                                   output_names=['out_file'],
                                                   function=gen_motion_parameters),
                                     name='calc_motion_parameters')
    pm.connect(inputNode, 'subject_id',
               calc_motion_parameters, 'subject_id')
    pm.connect(inputNode, 'scan_id',
               calc_motion_parameters, 'scan_id')
    pm.connect(inputNode, 'movement_parameters',
                calc_motion_parameters, 'movement_parameters')
    pm.connect(inputNode, 'max_displacement',
               calc_motion_parameters, 'max_displacement')
    
    pm.connect(calc_motion_parameters, 'out_file', 
               outputNode, 'motion_params')


    calc_power_parameters = pe.Node(util.Function(input_names=["subject_id", 
                                                                "scan_id", 
                                                                "FD_1D", 
                                                                "threshold",
                                                                "DVARS"],
                                                   output_names=['out_file'],
                                                   function=gen_power_parameters),
                                     name='calc_power_parameters')
    pm.connect(inputNode, 'subject_id',
               calc_power_parameters, 'subject_id')
    pm.connect(inputNode, 'scan_id',
               calc_power_parameters, 'scan_id')
    pm.connect(cal_DVARS, 'out_file',
               calc_power_parameters, 'DVARS')
    pm.connect(calculate_FD, 'out_file',
               calc_power_parameters, 'FD_1D')
    pm.connect(inputnode_threshold, 'threshold',
               calc_power_parameters, 'threshold')

    pm.connect(calc_power_parameters, 'out_file', 
               outputNode, 'power_params')


    return pm


def set_FD(in_file):
    """
    Method to calculate Framewise Displacement (FD) calculations
    
    Parameters
    ----------
    in_file : string
        movement parameters vector file path
    
    Returns
    -------
    out_file : string
        Frame -wise displalcement mat 
        file path
    
    """
    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'FD.1D')
    
    #Rotational displacements are converted from degrees to millimeters
    #by calculating displacement on the surface of a sphere of radius 50 mm
    cmd1 = sb.Popen(
        ['awk', '{x=$4} {y=$5} {z=$6} {a=$1} {b=$2} {c=$3} ' +
        '{print 2*3.142*50*(a/360),2*3.142*50*(b/360),' +
        ' 2*3.142*50*(c/360), x, y, z}',
        in_file], stdin=sb.PIPE, stdout=sb.PIPE,)

    #calculating relative displacement
    cmd2 = sb.Popen(
        ['awk', '{a=$1} {b=$2} {c=$3} {x=$4} {y=$5} {z=$6} ' +
        'NR>=1{print a-d, b-e, c-f, x-u, y-v, z-w}' +
        '{d=a} {e=b} {f=c} {u=x} {v=y} {w=z}'],
        stdin=cmd1.stdout, stdout=sb.PIPE,)

    #taking the absolute
    cmd3 = sb.Popen(
        ['awk', '{ for (i=1; i<=NF; i=i+1) {' +
        'if ($i < 0) $i = -$i} print}'],
        stdin=cmd2.stdout, stdout=sb.PIPE,)

    #summing the columns up
    cmd4 = sb.Popen(
        ['awk', '{a=$1+$2+$3+$4+$5+$6} {print a}'],
        stdin=cmd3.stdout, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd4.communicate()
    output = stdout_value.split()
    f = open(out_file, 'w')
    for out in output:
        print >>f, float(out)
    f.close()
    return out_file


def set_frames_ex(in_file, threshold):
    """
    Method to calculate Number of frames that would be censored 
    ("scrubbed") by removing the offending time frames 
    (i.e., those exceeding FD threshold), the preceeding frame, 
    and the two subsequent frames
    
    Parameters
    ----------
    in_file : string 
        framewise displacement(FD) file path
    threshold : float
         scrubbing thereshold value set in configuration file
    
    Returns
    -------
    out_file : string
        path to file containing offending time frames
    """

    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'frames_ex.1D')
    data = loadtxt(in_file)
    #masking zeroth timepoint value as 0, since the mean displacment value for
    #zeroth timepoint cannot be calculated, as there is no timepoint before it
    data[0] = 0

    extra_indices = []

    indices = [i[0] for i in (np.argwhere(data >= threshold)).tolist()]

    #adding addtional 2 after and one before timepoints
    for i in indices:
        if i > 0:
            extra_indices.append(i-1)
            if i+1 < data.size:
                extra_indices.append(i+1)
            if i+2 < data.size:
                extra_indices.append(i+2)

    indices = list(set(indices) | set(extra_indices))


    f = open(out_file, 'a')

    print "indices->", indices
    for idx in indices:
        f.write('%s,' % int(idx))

    f.close()

    return out_file

def set_frames_in(in_file, threshold, exclude_list):

    """
     Method to Calculate  the frames that are left
     after censoring for scrubbing.
     
     Parameters
     ----------
     in_file : string
        framewise displacement(FD) file path
     threshold : float
        scrubbing thereshold set in configuration file
     exclude_list : string
        path of file containing sensored timepoints
    
     Returns
     -------
     out_file : string 
        path of file containing remaining uncensored timepoints 
    """

    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'frames_in.1D')

    data = loadtxt(in_file)
    #masking zeroth timepoint value as 0, since the mean displacment value for
    #zeroth timepoint cannot be calculated, as there is no timepoint before it
    data[0] = 0

    indices = [i[0] for i in (np.argwhere(data < threshold)).tolist()]

    indx = []
    f = open(exclude_list, 'r')
    line = f.readline()
    if line:
        line = line.strip(',')
        indx = map(int, line.split(","))
    f.close()
    print indx

    if indx:
        indices = list(set(indices) - set(indx))

    f = open(out_file, 'a')

    for idx in indices:
        f.write('%s,' % int(idx))

    f.close()

    return out_file


def gen_motion_parameters(subject_id, scan_id, movement_parameters, max_displacement):
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
        path to csv file contianing various motion parameters

    """

    import os
    import numpy as np
    import re

    out_file = os.path.join(os.getcwd(), 'motion_parameters.txt')

    f = open(out_file, 'w')
    print >>f, "Subject,Scan,Mean_Relative_RMS_Displacement," \
        "Max_Relative_RMS_Displacement,Movements_gt_threshold,"\
        "Mean_Relative_Mean_Rotation,Mean_Relative_Maxdisp,Max_Relative_Maxdisp," \
        "Max_Abs_Maxdisp,Max Relative_Roll,Max_Relative_Pitch," \
        "Max_Relative_Yaw,Max_Relative_dS-I,Max_Relative_dL-R," \
        "Max_Relative_dP-A,Mean_Relative_Roll,Mean_Relative_Pitch,Mean_Relative_Yaw," \
        "Mean_Relative_dS-I,Mean_Relative_dL-R,Mean_Relative_dP-A,Max_Abs_Roll," \
        "Max_Abs_Pitch,Max_Abs_Yaw,Max_Abs_dS-I,Max_Abs_dL-R,Max_Abs_dP-A," \
        "Mean_Abs_Roll,Mean_Abs_Pitch,Mean_Abs_Yaw,Mean_Abs_dS-I,Mean_Abs_dL-R,Mean_Abs_dP-A"


    f.write("%s," % (subject_id))
    f.write("%s," % (scan_id))

    arr = np.genfromtxt(movement_parameters)
    arr = arr.T

    ##Relative RMS of translation
    rms = np.sqrt(arr[3]*arr[3] + arr[4]*arr[4] + arr[5]*arr[5])
    diff = np.diff(rms)
    MEANrms = np.mean(abs(diff))
    f.write("%.3f," % (MEANrms))

    #Max Relative RMS Displacement
    MAXrms = np.max(abs(diff))
    f.write("%.3f," % (MAXrms))

    ##NUMBER OF relative RMS movements >0.1mm
    NUMmove = np.sum(abs(diff) > 0.1)
    f.write("%.3f," % (NUMmove))

    ##Mean of mean relative rotation (params 1-3)
    MEANrot = np.mean(np.abs(np.diff((abs(arr[0])+ abs(arr[1])+ abs(arr[2]))/3 ) ) )
    f.write("%.3f," % (MEANrot))

    file = open(max_displacement, 'r')
    lines = file.readlines()
    file.close()
    list1 = []

    #remove any other information aother than matrix from
    #max displacement file. afni adds infomration to the file
    for l in lines:
        if re.match("^\d+?\.\d+?$", l.strip()):
            list1.append(float(l.strip()))

    arr2 = np.array(list1, dtype='float')

    #Mean Relative Maxdisp
    mean = np.mean(np.diff(arr2))
    f.write("%.3f," % (mean))

    #Max Relative Maxdisp
    relMAX = np.max(abs(np.diff(arr2)))
    f.write("%.3f," % (relMAX))

    #Max Abs Maxdisp
    MAX= np.max(arr2)
    f.write("%.3f," %(MAX))

    #Max Relative Roll,Max Relative Pitch,
    #Max Relative Yaw,Max Relative dS-I,
    #Max Relative dL-R,Max Relative dP-A
    for i in range(6):
        f.write("%.6f," %(np.max(abs(np.diff(arr[i])))))

    #Mean Relative Roll,Mean Relative Pitch,
    #Mean Relative Yaw,Mean Relative dS-I,
    #Mean Relative dL-R,Mean Relative dP-A
    for i in range(6):
        f.write("%.6f," %(np.mean(np.diff(arr[i]))))
        
    #Max Abs Roll,Max Abs Pitch,Max Abs Yaw,
    #Max Abs dS-I,Max Abs dL-R,Max Abs dP-A
    for i in range(6):
        f.write("%.6f," %(np.max(abs(arr[i]))))
    
    #Mean Abs Roll,Mean Abs Pitch,Mean Abs Yaw,
    #Mean Abs dS-I,Mean Abs dL-R,Mean Abs dP-A 
    for i in range(6):
        f.write("%.6f," %(np.mean(abs(arr[i]))))

    f.close()
    return out_file


def gen_power_parameters(subject_id, scan_id, FD_1D, DVARS, threshold = 1.0):
    
    """
    Method to generate Power parameters for scrubbing
    
    Parameters
    ----------
    subject_id : string
        subject name or id
    scan_id : string
        scan name or id
    FD_ID: string 
        framewise displacement file path
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

    import os
    import numpy as np
    from numpy import loadtxt

    out_file = os.path.join(os.getcwd(), 'pow_params.txt')

    f= open(out_file,'w')
    print >>f, "Subject,Scan,MeanFD,NumFD_greater_than_%.2f," \
    "rootMeanSquareFD,FDquartile(top1/4thFD),PercentFD_greater_than_%.2f," \
     "MeanDVARS"%(threshold,threshold)

    f.write("%s," % subject_id)
    f.write("%s," % scan_id)

    data= loadtxt(FD_1D)
    #Mean (across time/frames) of the absolute values 
    #for Framewise Displacement (FD)
    meanFD  = np.mean(data)
    f.write('%.4f,' % meanFD)
    
    #Number of frames (time points) where movement 
    #(FD) exceeded threshold
    numFD = float(data[data >threshold].size)
    f.write('%.4f,' % numFD)
    
    #Root mean square (RMS; across time/frames) 
    #of the absolute values for FD
    rmsFD = np.sqrt(np.mean(data))
    f.write('%.4f,' % rmsFD)

    #Mean of the top quartile of FD is $FDquartile
    quat=int(len(data)/4)
    FDquartile=np.mean(np.sort(data)[::-1][:quat])
    f.write('%.4f,' % FDquartile)

    ##NUMBER OF FRAMES >threshold FD as percentage of total num frames
    count = np.float(data[data>threshold].size)
    percentFD = (count*100/(len(data)+1))
    f.write('%.4f,' %percentFD)

    #Mean DVARS 
    meanDVARS = np.mean(np.load(DVARS))
    f.write('%.4f,' % meanDVARS)

    f.close()
    return out_file


def calculate_DVARS(rest, mask):
    """
    Method to calculate DVARS as per
    power's method
    
    Parameters
    ----------
    rest : string (nifti file)
        path to motion correct functional data
    mask : string (nifti file)
        path to brain only mask for functional data
        
    Returns
    -------
    out_file : string (numpy mat file)
        path to file containing  array of DVARS 
        calculation for each voxel
    """
    
    import numpy as np
    import nibabel as nib
    import os
    
    out_file = os.path.join(os.getcwd(), 'DVARS.npy')
    
    rest_data = nib.load(rest).get_data().astype(np.float32)
    mask_data = nib.load(mask).get_data().astype('bool')
    
    #square of relative intensity value for each voxel across
    #every timepoint 
    data = np.square(np.diff(rest_data, axis = 3))
    #applying mask, getting the data in the brain only
    data = data[mask_data]
    #square root and mean across all timepoints inside mask
    DVARS = np.sqrt(np.mean(data, axis=0))
    
    np.save(out_file, DVARS)
    
    return out_file
    
    
    
