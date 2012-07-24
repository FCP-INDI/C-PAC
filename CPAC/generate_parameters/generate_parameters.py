import CPAC.interfaces.afni.preprocess as e_afni
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def mov_pow_parameters():

    """
    The main purpose of this workflow is to get different statistical measures from the movement parameters obtained 
    in functional preprocessing and also generate power parameters required for scrubbing.
    
    Returns 
    -------
    param_wf : workflow object
          Workflow object containing various movement/motion and power parameters estimates.  
    
    Notes
    -----
    
    `Source <https://github.com/openconnectome/C-PAC/blob/master/CPAC/generate_parmeters/generate_parmeters.py>`_
    
    Workflow Inputs::
    
        inputspec.rest : func/rest file or a list of func/rest nifti file 
            User input functional(T2) Image, in any of the 8 orientations

        inputspec.subject_id : string 
            Subject name or id
            
        inputspec.scan_id : string
            Functional Scan id or name
            
        inputspec.max_displacement : string (Mat file)
            maximum displacement (in mm) vector for brain voxels in each volume.
            This file is obtained in functional preprocessing step
        
        inputspec.movement_parameters : string (Mat file)
            1D file containing six movement/motion parameters(3 Translation, 3 Rotations) 
            in different columns (roll pitch yaw dS  dL  dP), obtained in functional preprocessing step
        
        threshold_input.threshold : float
            scrubbing threshold
        
    Workflow Outputs::
    
        outputspec.mean_deriv_sq_1D : string (1D file)
          Mean over the masked brain region
          of the derivative of the raw data  
        
        outputspec.mean_raw_sq_1D: string (1D file)
          Mean over the masked brain region 
          for the raw data 
        
        outputspec.FD_1D : 1D file
            mean Framewise Displacement (FD)
        
        outputspec.sqrt_mean_deriv_sq_1D : 1D file
            RMS values for each volume of the 
            derivative file of the raw data
        
        outputspec.sqrt_mean_raw_sq_1D : 1D file
            RMS values for each volume of 
            the raw data
        
        outputspec.frames_ex_1D : 1D file
            Number of frames that would be censored ("scrubbed")
            also removing the offending time frames (i.e., those exceeding the threshold), 
            the preceeding frame, and the two subsequent frames
        
        outputspec.frames_in_1D :
            Number of frames left after removing for scrubbing
        
        outputspec.power_params : txt file
            Text file various power parameters for scrubbing.
        
        outputspec.motion_params : txt file
           Text file containing various movement parameters
        
        outputspec.ftof_percent_change_1D : 1D file
            The frame-wise RMS of the temporal derivative of the RAW (unprocessed) 
            data as a percentage of the RMS of the first of the two time points.
    
    Order of commands:
    
    - Calculate Frame Wise Displacement FD
    
      Differentiating head realignment parameters across frames yields a six dimensional timeseries that represents instantaneous head motion.   
      
      .. math:: FD_i = | \Delta d_{ix} | + | \Delta d_{iy} | + | \Delta d_{iz} | + | \Delta \\alpha_{i} | + | \Delta \\beta_{i} | + | \Delta \gamma_{i} |
      
      where  :math:`\Delta d_{ix} = d_{(i-1)x} - d_{ix}` 
    
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
      signal, which is the average value of a brain image at a timepoint).
      
      .. math:: \sqrt { \\left \\langle \\left [ \Delta  I_i \\left ( \\vec {x} \\right ) \\right ]^2 \\right \\rangle } = \sqrt { \\left \\langle \\left [ I_i \\left ( \\vec {x} \\right )  - I_{i-1} \\left ( \\vec {x} \\right ) \\right ]^2  \\right \\rangle}
      
      where, :math:`I_i \\left ( x \\right )` is image intensity at locus x on frame i and angle brackets denote the spatial average over the whole brain [R5]
                   
        - Create a dilated brain mask::
             
            3dAutomask -dilate 1  
                       -prefix mask 
                        rest.nii.gz
        
        - Calculate the difference in Image intensities between I and I-1 frame::
        
            3dcalc -a rest.nii.gz[4..299]  
                   -b rest.nii.gz[3..298] 
                   -expr '(a-b)' 
                   -prefix temp_deriv
        
        - Square the Intensity values of the derivative file obtained in above step::
         
            3dcalc -a temp_deriv 
                   -expr 'a*a' 
                   -prefix temp_deriv_sq
        
        - Calculate the mean intensity statistics over masked(whole brain) region
          for each volume of the squared derivative file::
        
            3dROIstats -quiet 
                       -mask mask
                       temp_deriv_sq
        
        - Square the intensity values of raw image::

            3dcalc -a rest.nii.gz[3..298] 
                   -expr 'a*a' 
                   -prefix raw_sq
        
        - Get the mean intensity over masked region for each volume of the
          squared raw file::
        
            3dROIstats -quiet 
                       -mask mask
                        raw_sq
    
        - calculate Square root of Mean value for derivative Image 
            
        - calculate square root of Mean for the Raw Image
    
      
    - Calculate Power parameters::
        
        MeanFD : Mean (across time/frames) of the absolute values for Framewise Displacement (FD), 
        computed as described in Power et al., Neuroimage, 2012)
        
        rootMeanSquareFD : Root mean square (RMS; across time/frames) of the absolute values for FD
        
        NumFD >=threshold : Number of frames (time points) where movement (FD) exceeded threshold
        
        rmsFD : Root mean square (RMS; across time/frames) of the absolute values for FD
        
        FDquartile(top 1/4th FD) : Mean of the top 25% highest FD values
        
        PercentFD( > threshold) : Number of frames (time points) where movement (FD) exceeded threshold 
                                  expressed as a percentage of the total number of frames (time points)
        
        Num5 : Number of frames (time points) where modDVARS exceeded 5%, expressed as a percentage of 
               the total number of frames (time points)
        
        Num10 : Number of frames (time points) where modDVARS exceeded 10%, expressed as a percentage of 
                the total number of frames (time points)
        
        MeanDVARS : The mean (across time) of a measure similar to DVARS as described by Power et al., 2012a. 
                    The frame-wise RMS of the temporal derivative of the RAW (unprocessed) data as a percentage of the 
                    RMS of the first of the two time points
        
        MeanDVARS_POW : MeanDVARS as mean of frame wise RMS of the RAW data. 
            
    - Calculate Motion Parameters
        
      Following motion parameters are calculated::
         
        Subject, Scan, Mean Relative RMS Displacement, Max Relative RSM Displacement,
        Movements >threshold, Mean Relative Mean Rotation, Mean Relative Maxdisp,
        Max Relative Maxdisp, Max Abs Maxdisp, Max Relative Roll,Max Relative Pitch,
        Max Relative Yaw, Max Relative dS-I, Max Relative dL-R,Max Relative dP-A,
        Mean Relative Roll, Mean Relative Pitch,Mean Relative Yaw, Mean Relative dS-I,
        Mean Relative dL-R, Mean Relative dP-A, Max Abs Roll, Max Abs Pitch, Max Abs Yaw,
        Max Abs dS-I, Max Abs dL-R, Max Abs dP-A

    
    High Level Workflow Graph:
    
    .. image:: ../images/parameters.dot.png
       :width: 500
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/parameters_detailed.dot.png
       :width: 500

    Examples
    --------
    
    >>> import generate_parameters
    >>> wf = generate_parameters.mov_pow_parameters()
    >>> wf.inputs.inputspec.movement_parameters = 'CPAC_outupts/sym_link/sub01/rest_1/func/rest_mc.1D'
    >>> wf.inputs.inputspec.max_displacement = 'CPAC_outputs/sym_link/sub01/rest_1/func/max_disp.1D'
    >>> wf.inputs.inputspec.subject_id = 'sub01'
    >>> wf.inputs.inputspec.scan_id = 'rest_1'
    >>> wf.inputs.threshold_input.threshold = 0.5
    >>> wf.base_dir = './working_dir'
    >>> wf.run()
    
    
    See Also
    --------
        

    
    References
    ----------
    
    .. [1] Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious 
           but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3),
           2142-2154. doi:10.1016/j.neuroimage.2011.10.018
           
    .. [2] Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Steps
           toward optimizing motion artifact removal in functional connectivity MRI; a reply to Carp.
           NeuroImage. doi:10.1016/j.neuroimage.2012.03.017
    
 

     
    """
    pm = pe.Workflow(name='param_wf')
    inputNode = pe.Node(util.IdentityInterface(fields=['subject_id',
                                                       'scan_id',
                                                       'rest',
                                                       'movement_parameters',
                                                       'max_displacement'
                                                    ]),
                        name='inputspec')

    inputnode_threshold = pe.Node(util.IdentityInterface(fields=['threshold']),
                             name='threshold_input')

    outputNode = pe.Node(util.IdentityInterface(fields=['mean_deriv_sq_1D',
                                                        'mean_raw_sq_1D',
                                                        'FD_1D',
                                                        'sqrt_mean_deriv_sq_1D',
                                                        'sqrt_mean_raw_sq_1D',
                                                        'frames_ex_1D',
                                                        'frames_in_1D',
                                                        'power_params',
                                                        'motion_params',
                                                        'ftof_percent_change_1D']),
                        name='outputspec')

    NVOLS = pe.Node(util.Function(input_names=['in_files'],
                                  output_names=['nvols'],
                                  function=get_img_nvols),
                    name='NVOLS')

    last = pe.MapNode(util.Function(input_names=['nvols'], 
                                              output_names=['last_volume'],
                                              function=last_vol), 
                                name='last', 
                                iterfield=['nvols'])

    last_minus_one = pe.MapNode(util.Function(input_names=['nvols', 'stopIdx', 'startIdx'], 
                                              output_names=['last_vol_minus_one'],
                                              function=trend_minus1), 
                                name='last_minus_one', 
                                iterfield=['nvols'])

    calculate_FD = pe.MapNode(util.Function(input_names=['in_file'], 
                                           output_names=['out_file'],
                                           function=set_FD), 
                             name='calculate_FD', 
                             iterfield=["in_file"])


    exclude_frames = pe.MapNode(util.Function(input_names=['in_file', 'threshold'], 
                                           output_names=['out_file'],
                                           function=set_frames_ex), 
                             name='exclude_frames', 
                             iterfield=["in_file"])

    include_frames = pe.MapNode(util.Function(input_names=['in_file', 'threshold', 'exclude_list'], 
                                           output_names=['out_file'],
                                           function=set_frames_in), 
                             name='include_frames', 
                             iterfield=["in_file", "exclude_list"])

    meanDVARS_perc_change = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'], 
                                                    output_names=['out_file'],
                                                    function=set_ftof_percent_change), 
                                      name='meanDVARS_perc_change', 
                                      iterfield=["infile_a", "infile_b"])

    DVARS_deriv = pe.MapNode(util.Function(input_names=['in_file'], output_names=['out_file'],
                                                function=set_sqrtmean_deriv), 
                                  name='DVARS_deriv',
                                  iterfield=["in_file"])

    DVARS_raw = pe.MapNode(util.Function(input_names=['in_file'], 
                                              output_names=['out_file'],
                                              function=set_sqrtmean_raw), 
                                name='DVARS_raw', 
                                iterfield=["in_file"])

    DVARS_calc1 = pe.MapNode(interface=e_afni.Threedcalc(), 
                          name='DVARS_calc1',
                          iterfield=["infile_a", "stop_idx", "infile_b", "stop_idx2"])
    DVARS_calc1.inputs.start_idx = 4
    DVARS_calc1.inputs.start_idx2 = 3
    DVARS_calc1.inputs.expr = '\'(a-b)\''
    DVARS_calc1.inputs.out_file = 'temp_deriv'

    DVARS_calc2 = pe.MapNode(interface=e_afni.Threedcalc(), 
                          name='DVARS_calc2', 
                          iterfield=["infile_a"])
    DVARS_calc2.inputs.expr = '\'a*a\''
    DVARS_calc2.inputs.out_file = 'temp_deriv_sq'

    DVARS_calc3 = pe.MapNode(interface=e_afni.Threedcalc(), 
                          name='DVARS_calc3', 
                          iterfield=["infile_a", "stop_idx"])
    DVARS_calc3.inputs.start_idx = 3
    DVARS_calc3.inputs.expr = '\'a*a\''
    DVARS_calc3.inputs.out_file = 'raw_sq'


    get_mask = pe.MapNode(interface=e_afni.ThreedAutomask(), 
                             name='get_mask', 
                             iterfield=["in_file"])
    get_mask.inputs.dilate = 1
    get_mask.inputs.genbrickhead = True
    get_mask.inputs.out_file = './mask'


    mean_derv = pe.MapNode(interface=e_afni.ThreedROIstats(), 
                                 name='mean_derv',
                                 iterfield=["in_file", "mask"])
    mean_derv.inputs.quiet = True

    mean_raw = pe.MapNode(interface=e_afni.ThreedROIstats(), 
                                 name='mean_raw',
                                 iterfield=["in_file", "mask"])
    mean_raw.inputs.quiet = True


    calc_motion_parameters = pe.MapNode(util.Function(input_names=["subject_id","scan_id", "rest", "movement_parameters", 
                                                                "max_displacement"],
                                                   output_names=['out_file'],
                                                   function=gen_motion_parameters),
                                     name='calc_motion_parameters',
                                     iterfield=["rest", 
                                                "movement_parameters", 
                                                "max_displacement"])

    calc_power_parameters = pe.MapNode(util.Function(input_names=["subject_id","scan_id","rest", "FD_1D", "threshold",
                                                               "ftof_percent", "sqrt_mean_raw"],
                                                   output_names=['out_file'],
                                                   function=gen_power_parameters),
                                     name='calc_power_parameters',
                                     iterfield=["rest", "FD_1D", 
                                                "ftof_percent", 
                                                "sqrt_mean_raw"])

    pm.connect(inputNode, 'rest', NVOLS, 'in_files')

    pm.connect(inputNode, 'rest', DVARS_calc1, 'infile_a')
    pm.connect(inputNode, 'rest', DVARS_calc1, 'infile_b')
    pm.connect(NVOLS, 'nvols', last, 'nvols')
    pm.connect(last, 'last_volume', DVARS_calc1, 'stop_idx')
    pm.connect(NVOLS, 'nvols', last_minus_one, 'nvols')
    pm.connect(last_minus_one, 'last_vol_minus_one', DVARS_calc1, 'stop_idx2')

    pm.connect(DVARS_calc1, 'brik_file', DVARS_calc2, 'infile_a')

    pm.connect(inputNode, 'rest', DVARS_calc3, 'infile_a')
    pm.connect(last_minus_one, 'last_vol_minus_one', DVARS_calc3, 'stop_idx')

    pm.connect(inputNode, 'rest', get_mask, 'in_file')

    pm.connect(get_mask, 'brik_file', mean_derv, 'mask')
    pm.connect(DVARS_calc2, 'brik_file', mean_derv, 'in_file')

    pm.connect(get_mask, 'brik_file', mean_raw, 'mask')
    pm.connect(DVARS_calc3, 'brik_file', mean_raw, 'in_file')

    pm.connect(mean_derv, 'stats', DVARS_deriv, 'in_file' )
    pm.connect(mean_raw, 'stats', DVARS_raw, 'in_file')

    pm.connect(DVARS_deriv, 'out_file', meanDVARS_perc_change, 'infile_a' )
    pm.connect(DVARS_raw, 'out_file', meanDVARS_perc_change, 'infile_b')

    ###Calculating mean Framewise Displacement
    pm.connect(inputNode, 'movement_parameters', calculate_FD, 'in_file' )
    
    ##calculating frames to exclude and include after scrubbing
    pm.connect(calculate_FD, 'out_file', exclude_frames, 'in_file')
    pm.connect(inputnode_threshold, 'threshold', exclude_frames, 'threshold')


    pm.connect(calculate_FD, 'out_file', include_frames, 'in_file')
    pm.connect(inputnode_threshold, 'threshold', include_frames, 'threshold')
    pm.connect(exclude_frames, 'out_file', include_frames, 'exclude_list')

    pm.connect(inputNode, 'subject_id', 
               calc_motion_parameters, 'subject_id')
    pm.connect(inputNode, 'scan_id', 
               calc_motion_parameters, 'scan_id')
    pm.connect(inputNode, 'rest', 
               calc_motion_parameters, 'rest')
    pm.connect(inputNode, 'movement_parameters', 
                calc_motion_parameters, 'movement_parameters')
    pm.connect(inputNode, 'max_displacement',
               calc_motion_parameters, 'max_displacement')

    pm.connect(inputNode, 'subject_id', 
               calc_power_parameters, 'subject_id')
    pm.connect(inputNode, 'scan_id', 
               calc_power_parameters, 'scan_id')
    pm.connect(inputNode, 'rest',  
               calc_power_parameters, 'rest')
    pm.connect(meanDVARS_perc_change, 'out_file', 
               calc_power_parameters, 'ftof_percent')
    pm.connect(calculate_FD, 'out_file',
               calc_power_parameters, 'FD_1D')
    pm.connect(inputnode_threshold, 'threshold',
               calc_power_parameters, 'threshold')
    pm.connect(DVARS_raw, 'out_file',
                calc_power_parameters, 'sqrt_mean_raw')

    pm.connect(mean_derv, 'stats', outputNode, 'mean_deriv_sq_1D')
    pm.connect(mean_raw, 'stats', outputNode, 'mean_raw_sq_1D')

    pm.connect(calculate_FD, 'out_file', outputNode, 'FD_1D')
    pm.connect(DVARS_deriv, 'out_file', outputNode, 'sqrt_mean_deriv_sq_1D')
    pm.connect(DVARS_raw, 'out_file', outputNode, 'sqrt_mean_raw_sq_1D')
    pm.connect(exclude_frames, 'out_file', outputNode, 'frames_ex_1D')
    pm.connect(include_frames, 'out_file', outputNode, 'frames_in_1D')
    pm.connect(meanDVARS_perc_change, 'out_file', outputNode, 'ftof_percent_change_1D')
    pm.connect(calc_motion_parameters, 'out_file', outputNode, 'motion_params')
    pm.connect(calc_power_parameters, 'out_file', outputNode, 'power_params')


    return pm


def last_vol(nvols):
    """
    Get the last volume of the image 
    
    Parameters
    ----------
    nvols : int
        no of scans
    
    Returns
    -------
    last_volume : int
        last scan
    """
    
    vol = nvols
    last_volume = (int(vol) - 1)
    return last_volume


def trend_minus1(nvols):
    """
    Get the second to last volume of the image
    
    Parameters
    ----------
    nvols : int
        no of scans
        
    Returns 
    -------
    last_vol_minus_one : int
        second to last scan
    
    """

    vol = nvols
    last_vol_minus_one = None

    last_volume = (int(vol) - 1)
    last_vol_minus_one = last_volume - 1

    return last_vol_minus_one


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


def set_sqrtmean_deriv(in_file):
    
    """
    Method to calculate mean of 
    rms of temporal derivative of raw data
    
    Parameters
    ----------
    in_file : string
        path to 1D file containing rms of temporal 
        derivative of raw data
    
    Returns
    -------
    out_file : string
        path to 1D file containing mean of rms of
        temporal derivative of raw data
    """
    import os
    import numpy 

    out_file = os.path.join(os.getcwd(), 'sqrt_mean_deriv_sq.1D')

    data = numpy.loadtxt(in_file)
    data = numpy.around(numpy.sqrt(data), decimals =6)
    numpy.savetxt(out_file, data)
    
    return out_file


def set_sqrtmean_raw(in_file):
    
    """
    Method to calculate mean of 
    rms of raw data
    
    Parameters
    ----------
    in_file : string
        path to 1D file containing rms 
        of raw data
    
    Returns
    -------
    out_file : string
        path to 1D file containing mean of rms 
        of raw data
    """

    import os
    import numpy
    
    out_file = os.path.join(os.getcwd(), 'sqrt_mean_raw_sq.1D')
    
    data = numpy.loadtxt(in_file)
    data = numpy.around(numpy.sqrt(data), decimals =6)
    numpy.savetxt(out_file, data)

    return out_file


def set_ftof_percent_change(infile_a, infile_b):
    
    """
    Method to calculate percentage change in the 
    RMS value of the timepoints
    
    Parameters
    ----------
    infile_a : string
        path to 1D file containing rms of temporal 
        derivative of raw data
    infile_b : string
         path to 1D file containing mean of rms 
         of raw data
    
    Returns
    -------
    out_file : string
        path to file containing percentage change in
        RMS value of timepoints
    
    """
    import subprocess as sb
    import os

    out_file = os.path.join(os.getcwd(), 'ftof_percent_change.1D')

    cmd1 = sb.Popen(['awk', 'NR==FNR{a[NR]=$1; next} {print a[FNR], $1}',
                     infile_a, infile_b], stdin=sb.PIPE, stdout=sb.PIPE,)

    cmd2 = sb.Popen(['awk', '{x=$1} {y=$2} {printf( "%.6f ", ((x/y)*100))}'],
                    stdin=cmd1.stdout, stdout=sb.PIPE,)

    stdout_value, stderr_value = cmd2.communicate()

    output = stdout_value.split()
    f = open(out_file, 'a')

    for out in output:
        print >>f, float(out)

    f.close()

    return out_file


def gen_motion_parameters(subject_id, scan_id, rest, movement_parameters, max_displacement):
    """
    Method to calculate all the movement parameters
    
    Parameters
    ----------
    subject_id : string
        subject name or id
    scan_id : string
        scan name or id
    rest : string 
        functional file path
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

    f= open(out_file,'w')
    print >>f, "Subject,Scan,Mean Relative RMS Displacement,"\
    "Max Relative RSM Displacement,Movements >threshold,Mean " \
    "Relative Mean Rotation,Mean Relative Maxdisp,Max Relative Maxdisp," \
    "Max Abs Maxdisp,Max Relative Roll,Max Relative Pitch," \
    "Max Relative Yaw,Max Relative dS-I,Max Relative dL-R," \
    "Max Relative dP-A,Mean Relative Roll,Mean Relative Pitch,Mean Relative Yaw," \
    "Mean Relative dS-I,Mean Relative dL-R,Mean Relative dP-A,Max Abs Roll," \
    "Max Abs Pitch,Max Abs Yaw,Max Abs dS-I,Max Abs dL-R,Max Abs dP-A"


    f.write("%s," %(subject_id))
    f.write("%s," %(scan_id))

    arr = np.genfromtxt(movement_parameters)
    arr = arr.T

    ##Relative RMS of translation
    rms= np.sqrt(arr[3]*arr[3] + arr[4]*arr[4] + arr[5]*arr[5])
    diff = np.diff(rms)
    MEANrms = np.mean(abs(diff))
    f.write("%.3f," %(MEANrms))
    
    #Max Relative RSM Displacement
    MAXrms= np.max(abs(diff))
    f.write("%.3f," %(MAXrms))

    ##NUMBER OF relative RMS movements >0.1mm
    NUMmove= np.sum(abs(diff)>0.1)
    f.write("%.3f," %(NUMmove))

    ##Mean of mean relative rotation (params 1-3)
    MEANrot= np.mean(np.abs(np.diff( (abs(arr[0])+ abs(arr[1])+ abs(arr[2]))/3 ) ) )
    f.write("%.3f," %(MEANrot))

    file = open(max_displacement, 'r')
    lines =file.readlines()
    file.close()
    list1=[]

    #remove any other information aother than matrix from
    #max displacement file. afni adds infomration to the file
    for l in lines:
        if re.match("^\d+?\.\d+?$", l.strip()):
            list1.append(float(l.strip()))

    arr2=np.array(list1, dtype='float')
    
    #Mean Relative Maxdisp
    mean=np.mean(np.diff(arr2))
    f.write("%.3f," %(mean))

    #Max Relative Maxdisp
    relMAX=np.max(abs(np.diff(arr2)))
    f.write("%.3f," %(relMAX))
    
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

    f.close()
    return out_file


def gen_power_parameters(subject_id, scan_id, rest, FD_1D, threshold, ftof_percent, sqrt_mean_raw):
    
    """
    Method to generate Power parameters for scrubbing
    
    Parameters
    ----------
    subject_id : string
        subject name or id
    scan_id : string
        scan name or id
    rest : string 
        functional file path
    FD_ID: string 
        framewise displacement file path
    threshold : float
        scrubbing threshold set in the configuration
    ftof_percent : string 
        path to 1D file containing percentage change in RMS value of the time points.
    sqrt_mean_raw : string 
        path to 1D file containing mean rms value of raw data for each scan
    
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
    print >>f, "Subject,Scan, MeanFD, rootMeanSquareFD, NumFD > %.2f," \
    "rmsFD, FDquartile(top 1/4th FD), PercentFD( >%.2f), Num5, \
    Num10, MeanDVARS, MeanDVARS_POW"%(threshold,threshold)

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


    data=loadtxt(ftof_percent)
    ###NUMBER OF relative FRAMES >5%
    num5= np.sum(data>=5)*100/len(data)
    f.write('%.4f,' % num5)
    
    #Number of frames (time points) where modDVARS exceeded 10%, 
    #expressed as a percentage of the total number of frames (time points)
    num10= np.sum(data>=10)*100/len(data)
    f.write('%.4f,' % num10)

    #Mean DVARS 
    meanDVARS  = np.mean(data)
    f.write('%.4f,' % meanDVARS)

    #mean DVARS as mean of rms value of raw data
    # as per powers paper
    data = loadtxt(sqrt_mean_raw)
    meanDVARS_POW = np.mean(data)
    f.write('%.4f,' % meanDVARS_POW)

    f.close()
    return out_file


def get_img_nvols(in_files):
    
    """
    Method to get number of volumes
    in the image 
    
    Parameters
    ----------
    in_files: string, list
        file path or list of file paths for 
        functional runs
    
    Returns 
    -------
    out : list
        list containing no of scans of 
        each functional run
    
    """

    out = []
    from nibabel import load
    if(isinstance(in_files, list)):
        for in_file in in_files:
            img = load(in_file)
            hdr = img.get_header()

            if len(hdr.get_data_shape()) > 3:
                nvols = int(hdr.get_data_shape()[3])
            else:
                nvols = 1
            out.append(nvols)
        return out

    else:
        img = load(in_files)
        hdr = img.get_header()
        if len(hdr.get_data_shape()) > 3:
            nvols = int(hdr.get_data_shape()[3])
        else:
            nvols = 1
        return [nvols]