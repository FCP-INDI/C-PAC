from nipype import logging
logger = logging.getLogger('workflow')

import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
#import CPAC.interfaces.afni.preprocess as preprocess
from nipype.interfaces.afni import preprocess
from CPAC.utils import dbg_file_lineno

# workflow to edit the scan to the proscribed TRs
def create_wf_edit_func( wf_name = "edit_func" ):
    """
    Workflow Inputs::
    
        inputspec.func : func file or a list of func/rest nifti file 
            User input functional(T2*) Image
            
        inputspec.start_idx : string 
            Starting volume/slice of the functional image (optional)
            
        inputspec.stop_idx : string
            Last volume/slice of the functional image (optional)
            
    Workflow Outputs::
    
        outputspec.edited_func : string (nifti file)
            Path to Output image with the initial few slices dropped
          
           
    Order of commands:
    
    - Get the start and the end volume index of the functional run. If not defined by the user, return the first and last volume.
    
        get_idx(in_files, stop_idx, start_idx)
        
    - Dropping the initial TRs. For details see `3dcalc <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_::
        
        3dcalc -a rest.nii.gz[4..299] 
               -expr 'a' 
               -prefix rest_3dc.nii.gz
               
    """

    # allocate a workflow object
    try:
        preproc = pe.Workflow(name=wf_name)
    except:
        logger.info( "Error allocating workflow %s."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

    # configure the workflow's input spec
    try:
        inputNode = pe.Node(util.IdentityInterface(fields=['func',
                                                           'start_idx',
                                                           'stop_idx']),
                            name='inputspec')
    except:
        logger.info( "Error allocating inputspec (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

    # configure the workflow's output spec
    try:
        outputNode = pe.Node(util.IdentityInterface(fields=['edited_func']),
                          name='outputspec')
    except:
        logger.info( "Error allocating output spec (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

    # allocate a node to check that the requested edits are
    # reasonable given the data
    try:
        func_get_idx = pe.Node(util.Function(input_names=['in_files', 
                                                          'stop_idx', 
                                                          'start_idx'],
                                   output_names=['stopidx', 
                                                 'startidx'],
                     function=get_idx), name='func_get_idx')
    except:
        logger.info( "Error allocating get_idx function node (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

   
    # wire in the func_get_idx node
    try: 
        preproc.connect(inputNode, 'func',
                        func_get_idx, 'in_files')
    except:
        logger.info( "Error connecting 'in_files' input to get_idx function node (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

    try:
        preproc.connect(inputNode, 'start_idx',
                        func_get_idx, 'start_idx')
    except:
        logger.info( "Error connecting 'start_idx' input to get_idx function node (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

    try:
        preproc.connect(inputNode, 'stop_idx',
                     func_get_idx, 'stop_idx')
    except:
        logger.info( "Error connecting 'stop_idx' input to get_idx function node (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

    try:
        # allocate a node to edit the functional file 
        func_drop_trs = pe.Node(interface=preprocess.Calc(),
                               name='func_drop_trs')
        func_drop_trs.inputs.expr = 'a'
        func_drop_trs.inputs.outputtype = 'NIFTI_GZ'
    except:
        logger.info( "Error allocating afni Calc node (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

   
    # wire in the inpus
    try: 
        preproc.connect(inputNode, 'func',
                        func_drop_trs, 'in_file_a')
    except:
        logger.info( "Error connecting 'in_file_a' input to afni Calc node (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

    try: 
        preproc.connect(func_get_idx, 'startidx',
                        func_drop_trs, 'start_idx')
    except:
        logger.info( "Error connecting 'start_idx' input to afni Calc node (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise

    try: 
        preproc.connect(func_get_idx, 'stopidx',
                        func_drop_trs, 'stop_idx')
    except:
        logger.info( "Error connecting 'stop_idx' input to afni Calc node (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise
  
    try: 
        # wire the output 
        preproc.connect(func_drop_trs, 'out_file',
                        outputNode, 'edited_func')
    except:
        logger.info( "Error connecting output (wflow %s)."+\
                     " (%s:%d)" % (wf_name, dbg_file_lineno() ))
        raise
  
    return preproc
    


#functional preprocessing
def create_func_preproc(use_bet = False, wf_name = 'func_preproc'):
    """
    
    The main purpose of this workflow is to process functional data. Raw rest file is deobliqued and reoriented 
    into RPI. Then take the mean intensity values over all time points for each voxel and use this image 
    to calculate motion parameters. The image is then skullstripped, normalized and a processed mask is 
    obtained to use it further in Image analysis.
    
    Parameters
    ----------
    
    wf_name : string
        Workflow name
    
    Returns 
    -------
    func_preproc : workflow object
        Functional Preprocessing workflow object
    
    Notes
    -----
    
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/func_preproc/func_preproc.py>`_
    
    Workflow Inputs::
    
        inputspec.rest : func/rest file or a list of func/rest nifti file 
            User input functional(T2) Image, in any of the 8 orientations
            
        scan_params.tr : string
            Subject TR
        
        scan_params.acquistion : string
            Acquisition pattern (interleaved/sequential, ascending/descending)
    
        scan_params.ref_slice : integer
            Reference slice for slice timing correction
            
    Workflow Outputs::
    
        outputspec.refit : string (nifti file)
            Path to deobliqued anatomical data 
        
        outputspec.reorient : string (nifti file)
            Path to RPI oriented anatomical data 
        
        outputspec.motion_correct_ref : string (nifti file)
             Path to Mean intensity Motion corrected image 
             (base reference image for the second motion correction run)
        
        outputspec.motion_correct : string (nifti file)
            Path to motion corrected output file
        
        outputspec.max_displacement : string (Mat file)
            Path to maximum displacement (in mm) for brain voxels in each volume
        
        outputspec.movement_parameters : string (Mat file)
            Path to 1D file containing six movement/motion parameters(3 Translation, 3 Rotations) 
            in different columns (roll pitch yaw dS  dL  dP)
        
        outputspec.skullstrip : string (nifti file)
            Path to skull stripped Motion Corrected Image 
        
        outputspec.mask : string (nifti file)
            Path to brain-only mask
            
        outputspec.example_func : string (nifti file)
            Mean, Skull Stripped, Motion Corrected output T2 Image path
            (Image with mean intensity values across voxels) 
        
        outputpsec.preprocessed : string (nifti file)
            output skull stripped, motion corrected T2 image 
            with normalized intensity values 

        outputspec.preprocessed_mask : string (nifti file)
           Mask obtained from normalized preprocessed image
           
    Order of commands:
               
    - Deobliqing the scans.  For details see `3drefit <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3drefit.html>`_::
    
        3drefit -deoblique rest_3dc.nii.gz
        
    - Re-orienting the Image into Right-to-Left Posterior-to-Anterior Inferior-to-Superior (RPI) orientation. For details see `3dresample <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dresample.html>`_::
    
        3dresample -orient RPI 
                   -prefix rest_3dc_RPI.nii.gz 
                   -inset rest_3dc.nii.gz
        
    - Calculate voxel wise statistics. Get the RPI Image with mean intensity values over all timepoints for each voxel. For details see `3dTstat <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTstat.html>`_::
    
        3dTstat -mean 
                -prefix rest_3dc_RPI_3dT.nii.gz 
                rest_3dc_RPI.nii.gz
    
    - Motion Correction. For details see `3dvolreg <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dvolreg.html>`_::  
       
        3dvolreg -Fourier 
                 -twopass 
                 -base rest_3dc_RPI_3dT.nii.gz/
                 -zpad 4 
                 -maxdisp1D rest_3dc_RPI_3dvmd1D.1D 
                 -1Dfile rest_3dc_RPI_3dv1D.1D 
                 -prefix rest_3dc_RPI_3dv.nii.gz 
                 rest_3dc_RPI.nii.gz
                 
      The base image or the reference image is the mean intensity RPI image obtained in the above the step.For each volume 
      in RPI-oriented T2 image, the command, aligns the image with the base mean image and calculates the motion, displacement 
      and movement parameters. It also outputs the aligned 4D volume and movement and displacement parameters for each volume.
                 
    - Calculate voxel wise statistics. Get the motion corrected output Image from the above step, with mean intensity values over all timepoints for each voxel. 
      For details see `3dTstat <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTstat.html>`_::
    
        3dTstat -mean 
                -prefix rest_3dc_RPI_3dv_3dT.nii.gz 
                rest_3dc_RPI_3dv.nii.gz
    
    - Motion Correction and get motion, movement and displacement parameters. For details see `3dvolreg <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dvolreg.html>`_::   

        3dvolreg -Fourier 
                 -twopass 
                 -base rest_3dc_RPI_3dv_3dT.nii.gz 
                 -zpad 4 
                 -maxdisp1D rest_3dc_RPI_3dvmd1D.1D 
                 -1Dfile rest_3dc_RPI_3dv1D.1D 
                 -prefix rest_3dc_RPI_3dv.nii.gz 
                 rest_3dc_RPI.nii.gz
        
      The base image or the reference image is the mean intensity motion corrected image obtained from the above the step (first 3dvolreg run). 
      For each volume in RPI-oriented T2 image, the command, aligns the image with the base mean image and calculates the motion, displacement 
      and movement parameters. It also outputs the aligned 4D volume and movement and displacement parameters for each volume.
    
    - Create a  brain-only mask. For details see `3dautomask <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutomask.html>`_::
    
        3dAutomask  
                   -prefix rest_3dc_RPI_3dv_automask.nii.gz 
                   rest_3dc_RPI_3dv.nii.gz

    - Edge Detect(remove skull) and get the brain only. For details see `3dcalc <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_::
    
        3dcalc -a rest_3dc_RPI_3dv.nii.gz 
               -b rest_3dc_RPI_3dv_automask.nii.gz 
               -expr 'a*b' 
               -prefix rest_3dc_RPI_3dv_3dc.nii.gz
    
    - Normalizing the image intensity values. For details see `fslmaths <http://www.fmrib.ox.ac.uk/fsl/avwutils/index.html>`_::
      
        fslmaths rest_3dc_RPI_3dv_3dc.nii.gz 
                 -ing 10000 rest_3dc_RPI_3dv_3dc_maths.nii.gz 
                 -odt float
                 
      Normalized intensity = (TrueValue*10000)/global4Dmean
                 
    - Calculate mean of skull stripped image. For details see `3dTstat <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTstat.html>`_::
        
        3dTstat -mean -prefix rest_3dc_RPI_3dv_3dc_3dT.nii.gz rest_3dc_RPI_3dv_3dc.nii.gz
        
    - Create Mask (Generate mask from Normalized data). For details see `fslmaths <http://www.fmrib.ox.ac.uk/fsl/avwutils/index.html>`_::
        
        fslmaths rest_3dc_RPI_3dv_3dc_maths.nii.gz 
               -Tmin -bin rest_3dc_RPI_3dv_3dc_maths_maths.nii.gz 
               -odt char

    High Level Workflow Graph:
    
    .. image:: ../images/func_preproc.dot.png
       :width: 1000
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/func_preproc_detailed.dot.png
       :width: 1000

    Examples
    --------
    
    >>> import func_preproc
    >>> preproc = create_func_preproc(bet=True)
    >>> preproc.inputs.inputspec.func='sub1/func/rest.nii.gz'
    >>> preproc.run() #doctest: +SKIP


    >>> import func_preproc
    >>> preproc = create_func_preproc(bet=False)
    >>> preproc.inputs.inputspec.func='sub1/func/rest.nii.gz'
    >>> preproc.run() #doctest: +SKIP
    
    """

    preproc = pe.Workflow(name=wf_name)
    inputNode = pe.Node(util.IdentityInterface(fields=['func']),
                        name='inputspec')
    
    outputNode = pe.Node(util.IdentityInterface(fields=['refit',
                                                        'reorient',
                                                        'reorient_mean',
                                                        'motion_correct',
                                                        'motion_correct_ref',
                                                        'movement_parameters',
                                                        'max_displacement',
                                                        #'xform_matrix',
                                                        'mask',
                                                        'skullstrip',
                                                        'example_func',
                                                        'preprocessed',
                                                        'preprocessed_mask',
                                                        'slice_time_corrected',
                                                        'oned_matrix_save']),

                          name='outputspec')

    
    func_deoblique = pe.Node(interface=preprocess.Refit(),
                            name='func_deoblique')
    func_deoblique.inputs.deoblique = True
    
    preproc.connect(inputNode, 'func',
                        func_deoblique, 'in_file')
    
    func_reorient = pe.Node(interface=preprocess.Resample(),
                               name='func_reorient')
    func_reorient.inputs.orientation = 'RPI'
    func_reorient.inputs.outputtype = 'NIFTI_GZ'

    preproc.connect(func_deoblique, 'out_file',
                    func_reorient, 'in_file')
    
    preproc.connect(func_reorient, 'out_file',
                    outputNode, 'reorient')
    
    func_get_mean_RPI = pe.Node(interface=preprocess.TStat(),
                            name='func_get_mean_RPI')
    func_get_mean_RPI.inputs.options = '-mean'
    func_get_mean_RPI.inputs.outputtype = 'NIFTI_GZ'
    
    preproc.connect(func_reorient, 'out_file',
                    func_get_mean_RPI, 'in_file')
        
    #calculate motion parameters
    func_motion_correct = pe.Node(interface=preprocess.Volreg(),
                             name='func_motion_correct')
    func_motion_correct.inputs.args = '-Fourier -twopass'
    func_motion_correct.inputs.zpad = 4
    func_motion_correct.inputs.outputtype = 'NIFTI_GZ'

    
    preproc.connect(func_reorient, 'out_file',
                    func_motion_correct, 'in_file')
    preproc.connect(func_get_mean_RPI, 'out_file',
                    func_motion_correct, 'basefile')


    func_get_mean_motion = func_get_mean_RPI.clone('func_get_mean_motion')
    preproc.connect(func_motion_correct, 'out_file',
                    func_get_mean_motion, 'in_file')
    
    preproc.connect(func_get_mean_motion, 'out_file',
                    outputNode, 'motion_correct_ref')
    
    
    func_motion_correct_A = func_motion_correct.clone('func_motion_correct_A')
    func_motion_correct_A.inputs.md1d_file = 'max_displacement.1D'
    
    preproc.connect(func_reorient, 'out_file',
                    func_motion_correct_A, 'in_file')
    preproc.connect(func_get_mean_motion, 'out_file',
                    func_motion_correct_A, 'basefile')
    
    preproc.connect(func_motion_correct_A, 'out_file',
                    outputNode, 'motion_correct')
    preproc.connect(func_motion_correct_A, 'md1d_file',
                    outputNode, 'max_displacement')
    preproc.connect(func_motion_correct_A, 'oned_file',
                    outputNode, 'movement_parameters')
    preproc.connect(func_motion_correct_A, 'oned_matrix_save',
                    outputNode, 'oned_matrix_save')


    if use_bet == False:
    
        func_get_brain_mask = pe.Node(interface=preprocess.Automask(),
                                   name='func_get_brain_mask')

        func_get_brain_mask.inputs.outputtype = 'NIFTI_GZ'

        preproc.connect(func_motion_correct_A, 'out_file',
                        func_get_brain_mask, 'in_file')
    
        preproc.connect(func_get_brain_mask, 'out_file',
                        outputNode, 'mask')

    else:

        func_get_brain_mask = pe.Node(interface=fsl.BET(),
                                      name='func_get_brain_mask_BET')

        func_get_brain_mask.inputs.mask = True
        func_get_brain_mask.inputs.functional = True

        erode_one_voxel = pe.Node(interface=fsl.ErodeImage(),
                                  name='erode_one_voxel')

        erode_one_voxel.inputs.kernel_shape = 'box'
        erode_one_voxel.inputs.kernel_size = 1.0

        preproc.connect(func_motion_correct_A, 'out_file',
                        func_get_brain_mask, 'in_file')

        preproc.connect(func_get_brain_mask, 'mask_file',
                        erode_one_voxel, 'in_file')

        preproc.connect(erode_one_voxel, 'out_file',
                        outputNode, 'mask')

        
    
    func_edge_detect = pe.Node(interface=preprocess.Calc(),
                            name='func_edge_detect')
    func_edge_detect.inputs.expr = 'a*b'
    func_edge_detect.inputs.outputtype = 'NIFTI_GZ'

    preproc.connect(func_motion_correct_A, 'out_file',
                    func_edge_detect, 'in_file_a')


    if use_bet == False:

        preproc.connect(func_get_brain_mask, 'out_file',
                        func_edge_detect, 'in_file_b')

    else:

        preproc.connect(erode_one_voxel, 'out_file',
                        func_edge_detect, 'in_file_b')



    preproc.connect(func_edge_detect, 'out_file',
                    outputNode, 'skullstrip')

    
    func_mean_skullstrip = pe.Node(interface=preprocess.TStat(),
                           name='func_mean_skullstrip')
    func_mean_skullstrip.inputs.options = '-mean'
    func_mean_skullstrip.inputs.outputtype = 'NIFTI_GZ'

    preproc.connect(func_edge_detect, 'out_file',
                    func_mean_skullstrip, 'in_file')
    
    preproc.connect(func_mean_skullstrip, 'out_file',
                    outputNode, 'example_func')
    
    
    func_normalize = pe.Node(interface=fsl.ImageMaths(),
                            name='func_normalize')
    func_normalize.inputs.op_string = '-ing 10000'
    func_normalize.inputs.out_data_type = 'float'

    preproc.connect(func_edge_detect, 'out_file',
                    func_normalize, 'in_file')
    
    preproc.connect(func_normalize, 'out_file',
                    outputNode, 'preprocessed')
    
    
    func_mask_normalize = pe.Node(interface=fsl.ImageMaths(),
                           name='func_mask_normalize')
    func_mask_normalize.inputs.op_string = '-Tmin -bin'
    func_mask_normalize.inputs.out_data_type = 'char'

    preproc.connect(func_normalize, 'out_file',
                    func_mask_normalize, 'in_file')
    
    preproc.connect(func_mask_normalize, 'out_file',
                    outputNode, 'preprocessed_mask')

    return preproc



def get_idx(in_files, stop_idx=None, start_idx=None):

    """
    Method to get the first and the last slice for
    the functional run. It verifies the user specified
    first and last slice. If the values are not valid, it 
    calculates and returns the very first and the last slice 
    
    Parameters
    ----------
    in_file : string (nifti file)
       Path to input functional run
        
    stop_idx : int
        Last volume to be considered, specified by user
        in the configuration file 
    
    stop_idx : int
        First volume to be considered, specified by user 
        in the configuration file 
    
    Returns
    -------
    stop_idx :  int
        Value of first slice to consider for the functional run 
        
    start_idx : int 
        Value of last slice to consider for the functional run
        
    """

    #stopidx = None
    #startidx = None
    # Import packages
    from nibabel import load

    # Init variables
    img = load(in_files)
    hdr = img.get_header()
    shape = hdr.get_data_shape()
    
    # Check to make sure the input file is 4-dimensional
    if len(shape) != 4:
        raise TypeError('Input nifti file: %s is not a 4D file' % in_files)
    # Grab the number of volumes
    nvols = int(hdr.get_data_shape()[3])

    if (start_idx == None) or (start_idx < 0) or (start_idx > (nvols - 1)):
        startidx = 0
    else:
        startidx = start_idx

    if (stop_idx == None) or (stop_idx > (nvols - 1)):
        stopidx = nvols - 1
    else:
        stopidx = stop_idx

    return stopidx, startidx


