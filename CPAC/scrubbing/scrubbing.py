import nipype.interfaces.afni.preprocess as e_afni
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def create_scrubbing_preproc(wf_name = 'scrubbing'):
    """
    This workflow essentially takes the list of offending timepoints that are to be removed
    and removes it from the motion corrected input image. Also, it removes the information
    of discarded time points from the movement parameters file obtained during motion correction.
    
    Parameters
    ----------
    wf_name : string
        Name of the workflow
    
    Returns
    -------
    scrub : object
        Scrubbing workfow object
    
    Notes
    -----
    `Source <https://github.com/openconnectome/C-PAC/blob/master/CPAC/scrubbing/scrubbing.py>`_
    
    Workflow Inputs::
        
        inputspec.frames_in_ID : string (mat file)
            path to file containing list of time points for which FD > threshold
        inputspec.movement_parameters : string (mat file)
            path to file containing 1D file containing six movement/motion parameters
            (3 Translation, 3 Rotations) in different columns 
        inputspec.preprocessed : string (nifti file)
            preprocessed input image path
            
    Workflow Outputs::
        
        outputspec.preprocessed : string (nifti file)
            preprocessed scrubbed output image 
        outputspec.scrubbed_movement_parameters : string (mat file)
            path to 1D file containing six movement/motion parameters
            for the timepoints which are not discarded by scrubbing
        
    Order of Commands:
    
    - Remove all movement parameters for all the time frames other than those that are present
      in the frames_in_1D file
      
    - Remove the discarded timepoints from the input image. For details see `3dcalc <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_::
        
        3dcalc -a bandpassed_demeaned_filtered.nii.gz[0,1,5,6,7,8,9,10,15,16,17,18,19,20,24,25,287,288,289,290,291,292,293,294,295] 
               -expr 'a' -prefix bandpassed_demeaned_filtered_3dc.nii.gz
               
    High Level Workflow Graph:
    
    .. exec::
        from CPAC.scrubbing import create_scrubbing_preproc
        wf = create_scrubbing_preproc()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/scrubbing.dot'
        )

    .. image:: ../../images/generated/scrubbing.png
       :width: 500
    
    Detailed Workflow Graph:
    
    .. image:: ../../images/generated/scrubbing_detailed.png
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

    scrub = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['frames_in_1D',
                                                       'movement_parameters',
                                                       'preprocessed'
                                                    ]),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['preprocessed',
                                                        'scrubbed_movement_parameters']),
                         name='outputspec')

    craft_scrub_input = pe.Node(util.Function(input_names=['scrub_input', 'frames_in_1D_file'],
                                              output_names=['scrub_input_string'],
                                              function=get_indx),
                                name = 'scrubbing_craft_input_string')

    scrubbed_movement_parameters = pe.Node(util.Function(input_names=['infile_a', 'infile_b'],
                                                 output_names=['out_file'],
                                                 function=get_mov_parameters),
                                   name='scrubbed_movement_parameters')

    # THIS commented out until Nipype has an input for this interface that
    # allows for the selection of specific volumes to include

    #scrubbed_preprocessed = pe.Node(interface=e_afni.Calc(),
    #                           name='scrubbed_preprocessed')
    #scrubbed_preprocessed.inputs.expr = 'a'
    #scrubbed_preprocessed.inputs.outputtype = 'NIFTI_GZ'

    scrubbed_preprocessed = pe.Node(util.Function(input_names=['scrub_input'],
                                                  output_names=['scrubbed_image'],
                                                  function=scrub_image),
                                    name='scrubbed_preprocessed')

    scrub.connect(inputNode, 'preprocessed', craft_scrub_input, 'scrub_input')
    scrub.connect(inputNode, 'frames_in_1D', craft_scrub_input, 'frames_in_1D_file')

    scrub.connect(craft_scrub_input, 'scrub_input_string', scrubbed_preprocessed, 'scrub_input')

    scrub.connect(inputNode, 'movement_parameters', scrubbed_movement_parameters, 'infile_b')
    scrub.connect(inputNode, 'frames_in_1D', scrubbed_movement_parameters, 'infile_a' )

    scrub.connect(scrubbed_preprocessed, 'scrubbed_image', outputNode, 'preprocessed')
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
    import warnings

    out_file = os.path.join(os.getcwd(), 'rest_mc_scrubbed.1D')

    f1= open(infile_a)
    f2=open(infile_b)
    l1=f1.readline()
    l2=f2.readlines()
    f1.close()
    f2.close()

    if l1:
        l1=l1.rstrip(',').split(',')
        warnings.warn("number of timepoints remaining after scrubbing -> %d"%len(l1))
    else:
        raise Exception("No time points remaining after scrubbing.")

    f = open(out_file, 'a')
    for l in l1:
        data=l2[int(l.strip())]
        f.write(data)
    f.close()
    return out_file


def get_indx(scrub_input, frames_in_1D_file):
    """
    Method to get the list of time 
    frames that are to be included
    
    Parameters
    ----------
    in_file : string
        path to file containing the valid time frames
    
    Returns
    -------
    scrub_input_string : string
        input string for 3dCalc in scrubbing workflow,
        looks something like " 4dfile.nii.gz[0,1,2,..100] "
    
    """

    with open(frames_in_1D_file, 'r') as f:
        line = f.readline()

    line = line.strip(',')
    if line:
        indx = map(int, line.split(","))
    else:
        raise Exception("No time points remaining after scrubbing.")

    scrub_input_string = scrub_input + str(indx).replace(" ", "")

    return scrub_input_string


def scrub_image(scrub_input):
    """
    Method to run 3dcalc in order to scrub the image. This is used instead of
    the Nipype interface for 3dcalc because functionality is needed for
    specifying an input file with specifically-selected volumes. For example:
        input.nii.gz[2,3,4,..98], etc.
        
    Parameters
    ----------
    scrub_input : string
        path to 4D file to be scrubbed, plus with selected volumes to be
        included
        
    Returns
    -------
    scrubbed_image : string
        path to the scrubbed 4D file
        
    """

    import os

    os.system("3dcalc -a %s -expr 'a' -prefix scrubbed_preprocessed.nii.gz" % scrub_input)

    scrubbed_image = os.path.join(os.getcwd(), "scrubbed_preprocessed.nii.gz")

    return scrubbed_image
