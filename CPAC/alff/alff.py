import os
import sys
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from CPAC.alff.alff import *
from CPAC.alff.utils import *
from nipype.interfaces.afni import preprocess

def create_alff(wf_name = 'alff_workflow'):

    """
    Calculate Amplitude of low frequency oscillations(ALFF) and fractional ALFF maps

    Parameters
    ----------
    wf_name : string
        Workflow name

    Returns
    -------
    alff_workflow : workflow object
        ALFF workflow

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/alff/alff.py>`_

    Workflow Inputs::

        hp_input.hp : list (float) 
            high pass frequencies

        lp_input.lp : list (float) 
            low pass frequencies

        inputspec.rest_res : string (existing nifti file)
            Nuisance signal regressed functional image

        inputspec.rest_mask : string (existing nifti file)
            A mask volume(derived by dilating the motion corrected functional volume) in native space


    Workflow Outputs::

        outputspec.alff_img : string (nifti file)
            outputs image containing the sum of the amplitudes in the low frequency band

        outputspec.falff_img : string (nifti file)
            outputs image containing the sum of the amplitudes in the low frequency band divided by the amplitude of the total frequency

        outputspec.alff_Z_img : string (nifti file)
            outputs image containing Normalized ALFF Z scores across full brain in native space

        outputspec.falff_Z_img : string (nifti file)
            outputs image containing Normalized fALFF Z scores across full brain in native space


    Order of Commands:

    - Filter the input file rest file( slice-time, motion corrected and nuisance regressed) ::
        3dBandpass -prefix residual_filtered.nii.gz 
                    0.009 0.08 residual.nii.gz
                    
    - Calculate ALFF by taking the standard deviation of the filtered file ::
        3dTstat -stdev 
                -mask rest_mask.nii.gz 
                -prefix residual_filtered_3dT.nii.gz
                residual_filtered.nii.gz
                  
    - Calculate the standard deviation of the unfiltered file ::
        3dTstat -stdev 
                -mask rest_mask.nii.gz 
                -prefix residual_3dT.nii.gz
                residual.nii.gz
  
    - Calculate fALFF ::
        3dcalc -a rest_mask.nii.gz 
               -b residual_filtered_3dT.nii.gz
               -c residual_3dT.nii.gz  
               -expr '(1.0*bool(a))*((1.0*b)/(1.0*c))' -float

    - Normalize ALFF/fALFF to Z-score across full brain ::

        fslstats
        ALFF.nii.gz
        -k rest_mask.nii.gz
        -m > mean_ALFF.txt ; mean=$( cat mean_ALFF.txt )

        fslstats
        ALFF.nii.gz
        -k rest_mask.nii.gz
        -s > std_ALFF.txt ; std=$( cat std_ALFF.txt )

        fslmaths
        ALFF.nii.gz
        -sub ${mean}
        -div ${std}
        -mas rest_mask.nii.gz ALFF_Z.nii.gz

        fslstats
        fALFF.nii.gz
        -k rest_mask.nii.gz
        -m > mean_fALFF.txt ; mean=$( cat mean_fALFF.txt )

        fslstats
        fALFF.nii.gz
        -k rest_mask.nii.gz 
        -s > std_fALFF.txt
        std=$( cat std_fALFF.txt )

        fslmaths
        fALFF.nii.gz
        -sub ${mean}
        -div ${std}
        -mas rest_mask.nii.gz
        fALFF_Z.nii.gz

    High Level Workflow Graph:
    
    .. image:: ../images/alff.dot.png
        :width: 500

    Detailed Workflow Graph:

    .. image:: ../images/alff_detailed.dot.png
        :width: 500

    
    References
    ----------

    .. [1] Zou, Q.-H., Zhu, C.-Z., Yang, Y., Zuo, X.-N., Long, X.-Y., Cao, Q.-J., Wang, Y.-F., et al. (2008). An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: fractional ALFF. Journal of neuroscience methods, 172(1), 137-41. doi:10.10

    Examples
    --------

    >>> alff_w = create_alff()
    >>> alff_w.inputs.hp_input.hp = [0.01]
    >>> alff_w.inputs.lp_input.lp = [0.1]
    >>> alff_w.get_node('hp_input').iterables = ('hp',[0.01])
    >>> alff_w.get_node('lp_input').iterables = ('lp',[0.1])
    >>> alff_w.inputs.inputspec.rest_res = '/home/data/subject/func/rest_bandpassed.nii.gz'
    >>> alff_w.inputs.inputspec.rest_mask= '/home/data/subject/func/rest_mask.nii.gz' 
    >>> alff_w.run() # doctest: +SKIP


    """

    wf = pe.Workflow(name= wf_name)
    inputNode = pe.Node(util.IdentityInterface(fields=['rest_res',
                                                       'rest_mask']),
                        name='inputspec')
    
    inputnode_hp = pe.Node(util.IdentityInterface(fields=['hp']),
                             name='hp_input')

    inputnode_lp = pe.Node(util.IdentityInterface(fields=['lp']),
                             name='lp_input')
    
    outputNode = pe.Node(util.IdentityInterface(fields=[ 'alff_img',
                                                         'falff_img']),
                          name='outputspec')
    
    #filtering
    bandpass = pe.Node(interface= preprocess.Bandpass(),
                       name = 'bandpass_filtering')
    bandpass.inputs.outputtype = 'NIFTI_GZ'
    bandpass.inputs.out_file = os.path.join(os.path.curdir, 'residual_filtered.nii.gz')
    wf.connect(inputnode_hp, 'hp',
                 bandpass, 'highpass')
    wf.connect(inputnode_lp, 'lp',
                 bandpass, 'lowpass')
    wf.connect(inputNode, 'rest_res',
                 bandpass, 'in_file') 
    
    get_option_string = pe.Node(util.Function(input_names = ['mask'],
                                              output_names = ['option_string'],
                                              function = get_opt_string), 
                                name = 'get_option_string')
    wf.connect(inputNode, 'rest_mask',
                 get_option_string, 'mask')
    
    #standard deviation over frequency
    stddev_fltrd = pe.Node(interface = preprocess.TStat(),
                            name = 'stddev_fltrd')
    stddev_fltrd.inputs.outputtype = 'NIFTI_GZ'
    stddev_fltrd.inputs.out_file = os.path.join(os.path.curdir, 'residual_filtered_3dT.nii.gz')
    wf.connect(bandpass, 'out_file',
               stddev_fltrd, 'in_file')
    wf.connect(get_option_string, 'option_string',
               stddev_fltrd, 'options') 
    
    wf.connect(stddev_fltrd, 'out_file',
                 outputNode, 'alff_img')
    
    #standard deviation of the unfiltered nuisance corrected image
    stddev_unfltrd = pe.Node(interface = preprocess.TStat(),
                            name = 'stddev_unfltrd')
    stddev_unfltrd.inputs.outputtype = 'NIFTI_GZ'
    stddev_unfltrd.inputs.out_file = os.path.join(os.path.curdir, 'residual_3dT.nii.gz')
    wf.connect(inputNode, 'rest_res',
                 stddev_unfltrd, 'in_file')
    wf.connect(get_option_string, 'option_string',
                 stddev_unfltrd, 'options') 
    
    #falff calculations
    falff = pe.Node(interface = preprocess.Calc(),
                    name = 'falff')
    falff.inputs.args = '-float'
    falff.inputs.expr = '(1.0*bool(a))*((1.0*b)/(1.0*c))'
    falff.inputs.outputtype = 'NIFTI_GZ'
    wf.connect(inputNode, 'rest_mask',
               falff, 'in_file_a')
    wf.connect(stddev_fltrd, 'out_file',
               falff, 'in_file_b')
    wf.connect(stddev_unfltrd, 'out_file',
               falff, 'in_file_c')
    
    wf.connect(falff, 'out_file',
               outputNode, 'falff_img') 
    
    
    return wf


def get_zscore(wf_name = 'z_score'):
    
    """
    Workflow to calculate z-scores
    
    Parameters
    ----------
    wf_name : string
        name of the workflow
        
    Returns
    -------
    wf : workflow object
    
    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/network_centrality/z_score.py>`_
    
    
    Workflow Inputs::
        
        inputspec.input_file : string
            path to input functional derivative file for which z score has to be calculated
        inputspec.mask_file : string
            path to whole brain functional mask file required to calculate zscore
    
    Workflow Outputs::
        
        outputspec.z_score_img : string
             path to image containing Normalized Input Image Z scores across full brain.
    
    High Level Workflow Graph:
    
    .. image:: ../images/zscore.dot.png
       :width: 500
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/zscore_detailed.dot.png
       :width: 500
    
    Example
    -------
    >>> import get_zscore as z
    >>> wf = z.get_zscore()
    >>> wf.inputs.inputspec.input_file = '/home/data/graph_working_dir/calculate_centrality/degree_centrality_binarize.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/graphs/GraphGeneration/new_mask_3m.nii.gz'
    >>> wf.run()
    
    """
    
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl
    
    wflow = pe.Workflow(name = wf_name)
    
    inputNode = pe.Node(util.IdentityInterface(fields=['input_file',
                                                       'mask_file']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['z_score_img']),
                          name='outputspec')

    mean = pe.Node(interface=fsl.ImageStats(),
                   name='mean')
    mean.inputs.op_string = '-k %s -m'    
    wflow.connect(inputNode, 'input_file',
                  mean, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  mean, 'mask_file')


    standard_deviation = pe.Node(interface=fsl.ImageStats(),
                                 name='standard_deviation')
    standard_deviation.inputs.op_string = '-k %s -s'
    wflow.connect(inputNode, 'input_file',
                  standard_deviation, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  standard_deviation, 'mask_file')
    
    
    op_string = pe.Node(util.Function(input_names=['mean','std_dev'],
                                      output_names=['op_string'],
                                      function=get_operand_string),
                        name='op_string')
    wflow.connect(mean, 'out_stat',
                  op_string, 'mean')
    wflow.connect(standard_deviation, 'out_stat',
                  op_string, 'std_dev')
    
    
    z_score = pe.Node(interface=fsl.MultiImageMaths(),
                        name='z_score')
    wflow.connect(op_string, 'op_string',
                  z_score, 'op_string')
    wflow.connect(inputNode, 'input_file',
                  z_score, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  z_score, 'operand_files')
    
    wflow.connect(z_score, 'out_file',
                  outputNode, 'z_score_img')
    
    return wflow
    

