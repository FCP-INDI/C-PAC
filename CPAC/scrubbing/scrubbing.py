#import e_afni
import CPAC.interfaces.afni.preprocess as e_afni
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def create_scrubbing_preproc():
    """
    This workflow essentially takes the list of offending timepoints that are to be removed
    and removes it from the motion corrected input image. Also, it removes the information
    of discarded time points from the movement parameters file obtained during motion correction.

    
    Returns
    -------
    scrub: object
        Scrubbing workfow object

    
    Notes
    -----
    `Source <https://github.com/openconnectome/C-PAC/blob/master/CPAC/scrubbing/scrubbing.py>`_

    
    Workflow Inputs::
        
        inputspec.frames_in_ID : file
            List of time points for which FD > threshold
        inputspec.movement_parameters : mat file
            1D file containing six movement/motion parameters
            (3 Translation, 3 Rotations) in different columns 
        inputspec.preprocessed : nifti file
            Preprocessed input image

            
    Workflow Outputs::
        
        outputspec.preprocessed : nifti file
            Preprocessed scrubbed output image 
        outputspec.scrubbed_movement_parameters : mat file
            1D file containing six movement/motion parameters
            for the timepoints which are not discarded by scrubbing

        
    Order of Commands:
    
    - Remove the movement parameters for all the time frames other than 
      those that are present in frames_in_1D file
      
    - Remove the discarded timepoints from the input image
        
        3dcalc -a bandpassed_demeaned_filtered.nii.gz[0,1,5,6,7,8,9,10,15,16,17,18,19,20,24,25,287,288,289,290,291,292,293,294,295] 
               -expr 'a' 
	       -prefix bandpassed_demeaned_filtered_3dc.nii.gz

               
    High Level Workflow Graph:
    
    .. image:: ../images/scrubbing.dot.png
       :width: 500
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/scrubbing_detailed.dot.png
       :width: 500
       
    Example
    -------
    >>> import scrubbing
    >>> sc = scrubbing.create_scrubbing_preproc()
    >>> sc.inputs.inputspec.frames_in_ID = 'frames_in.1D'
    >>> sc.inputs.inputpsec.movement_parameters = 'rest_mc.1D'
    >>> sc.inputs.inputpsec.preprocessed = 'rest_pp.nii.gz'
    >>> sc.run()  -- SKIP doctest
    
    """

    scrub = pe.Workflow(name='sc_preproc')

    inputNode = pe.Node(util.IdentityInterface(fields=['frames_in_1D',
                                                       'movement_parameters',
                                                       'preprocessed'
                                                    ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=['',
                                                         'scrubbed_movement_parameters']),
                        name='outputspec')


    scrubbed_movement_parameters = pe.MapNode(util.Function(input_names=['infile_a', 'infile_b'], 
                                                 output_names=['out_file'],
                                                 function=get_mov_parameters), 
                                   name='scrubbed_movement_parameters',
                                   iterfield=["infile_a", "infile_b"])

    scrubbed_preprocessed = pe.MapNode(interface=e_afni.Threedcalc(), 
                               name='scrubbed_preprocessed',
                               iterfield=["infile_a", "list_idx"] )
    scrubbed_preprocessed.inputs.expr = '\'a\''

    scrub.connect(inputNode, 'preprocessed', scrubbed_preprocessed, 'infile_a')
    scrub.connect(inputNode, ('frames_in_1D', get_indx), scrubbed_preprocessed, 'list_idx')

    scrub.connect(inputNode, 'movement_parameters', scrubbed_movement_parameters, 'infile_b')
    scrub.connect(inputNode, 'frames_in_1D', scrubbed_movement_parameters, 'infile_a' )

    scrub.connect(scrubbed_preprocessed, 'out_file', outputNode, 'preprocessed')
    scrub.connect(scrubbed_movement_parameters, 'out_file', outputNode, 'scrubbed_movement_parameters')

    return scrub


def get_mov_parameters(infile_a, infile_b):
    """
    Method to get the new movement parameters
    file after removing the offending time frames 
    (i.e., those exceeding FD 0.5mm/0.2mm threshold)
    
    Parameters
    ----------
    infile_a : string
        path to file containing the valid time frames
    
    infile_b : string
        path to the file containing  motion parameters 
    
    Returns
    -------
    out_file : string
        path to the file containing motion parameters
        for the valid time frames 
        
    """
    import os

    out_file = os.path.join(os.getcwd(), 'rest_mc_scrubbed.1D')

    f1= open(infile_a)
    f2=open(infile_b)
    l1=f1.readline()
    l2=f2.readlines()
    f1.close()
    f2.close()

    l1=l1.rstrip(',').split(',')

    f = open(out_file, 'a')
    for l in l1:
        data=l2[int(l.strip())]
        f.write(data)
    f.close()
    return out_file


def get_indx(in_file):
    """
    Method to get the list of time 
    frames that are to be included
    
    Parameters
    ----------
    in_file : string
        path to file containing the valid time frames
    
    Returns
    -------
    indx : list
        list of frame indexes
    
    """

    indx = []

    if(isinstance(in_file, list)):
        for file in in_file:
            f = open(file, 'r')
            line = f.readline()
            line = line.strip(',')
            indx.append(map(int, line.split(",")))
            f.close()
        print "indx ", indx
        return indx
    else:
            f = open(file, 'r')
            line = f.readline()
            line = line.strip(',')
            indx.append(map(int, line.split(",")))
            f.close()
            print "indx in else", indx
            return indx
