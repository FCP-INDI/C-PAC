import os
import sys
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from CPAC.alff.alff import *
from CPAC.alff.utils import *

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

    Examples
    --------


    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/alff/alff.py>`_

    Workflow Inputs: ::

        hp_input.hp : list (float) 
            high pass frequencies

        lp_input.lp : list (float) 
            low pass frequencies

        inputspec.rest_res : string (existing nifti file)
            Nuisance signal regressed functional image

        inputspec.rest_mask : string (existing nifti file)
            A mask volume(derived by dilating the motion corrected functional volume) in native space

        inputspec.tr : float
            scan TR, if input as None, TR is extracted from the nifti file header

    Workflow Outputs: ::

        outputspec.power_spectrum_distribution : string (nifti file)
            outputs image containing the spectral power density of residual functional image

        outputspec.alff_img : string (nifti file)
            outputs image containing the sum of the amplitudes in the low frequency band

        outputspec.falff_img : string (nifti file)
            outputs image containing the sum of the amplitudes in the low frequency band divided by the amplitude of the total frequency

        outputspec.alff_Z_img : string (nifti file)
            outputs image containing Normalized ALFF Z scores across full brain in native space

        outputspec.falff_Z_img : string (nifti file)
            outputs image containing Normalized fALFF Z scores across full brain in native space


    Order of Commands:

    - Compute the Power Spectrum ::

        fslpspec
        rest_res.nii.gz
        power_spectrum_distribution.nii.gz

    - Compute square root of power spectrum ::

        fslmaths
        power_spectrum_distribution.nii.gz
        -sqrt prealff_func_data_ps_sqrt.nii.gz

    - Extract the TR and Number of volumes(NVOLS) from the Nuisance signal regressed functional image

    - Calculate the High Frequency Point(n2) and Low Frequency Point(n1) ::

        n1 = (HP * NVOLS * TR) -1
        n2 = (LP * NVOLS * TR) - (HP * NVOLS * TR) + 1

    - Cut the low frequency data from the the whole frequency band ::

        fslroi
        prealff_func_data_ps_sqrt.nii.gz
        prealff_func_ps_slow.nii.gz
        n1
        n2

    - Calculate ALFF as the sum of the amplitudes in the low frequency band ::

        fslmaths 
        prealff_func_ps_slow.nii.gz
        -Tmean
        -mul n2
        ALFF.nii.gz

    - Compute amplitude of total frequency ::

        fslmaths
        prealff_func_data_ps_sqrt.nii.gz
        -Tmean
        -mul NVOLS
        -div 2
        prealff_func_pssum_amplitudes_low_frequency.nii.gz

    - Compute fALFF as ALFF/amplitude of total frequency ::

        fslmaths
        ALFF.nii.gz
        -div prealff_func_pssum_amplitudes_low_frequency.nii.gz
        fALFF.nii.gz

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

    Workflow ALFF and fractional ALFF:
    
    .. image:: ../images/alff_graph.dot.png
        :width: 500

    Workflow Detailed:

    .. image:: ../images/alff_detailed_graph.dot.png
        :width: 500

    
    References
    ----------

    .. [1] Zou, Q.-H., Zhu, C.-Z., Yang, Y., Zuo, X.-N., Long, X.-Y., Cao, Q.-J., Wang, Y.-F., et al. (2008). An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: fractional ALFF. Journal of neuroscience methods, 172(1), 137-41. doi:10.10

    Examples
    --------

    >>> alff_w = create_alff()
    >>> alff_w.inputs.hp_input.hp = [0.01]
    >>> alff_w.inputs.lp_input.lp = [0.1]
    >>> alff_w.get_node('hp_input').iterables = ('hp',
                                                [0.01])
    >>> alff_w.get_node('lp_input').iterables = ('lp',
                                                [0.1])
    >>> alff_w.inputs.inputspec.rest_res = '/home/data/subject/func/rest_bandpassed.nii.gz'
    >>> alff_w.inputs.inputspec.rest_mask= '/home/data/subject/func/rest_mask.nii.gz' 
    >>> alff_w.inputs.inputspec.tr = None
    >>> alff_w.run() # doctest: +SKIP


    """


    alff = pe.Workflow(name= wf_name)
    inputNode = pe.Node(util.IdentityInterface(fields=['rest_res',
                                                       'rest_mask',
                                                       'tr']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=[
                                            'power_spectrum_distribution',
                                            'alff_img',
                                            'falff_img',
                                            'alff_Z_img',
                                            'falff_Z_img'
                                            ]),
                          name='outputspec')



    inputnode_hp = pe.Node(util.IdentityInterface(fields=['hp']),
                             name='hp_input')

    inputnode_lp = pe.Node(util.IdentityInterface(fields=['lp']),
                             name='lp_input')

    TR = pe.Node(util.Function(input_names=['in_files', 'TRa'],
                               output_names=['TR'],
                 function=get_img_tr), name='TR')
    #TR.inputs.TRa = tr

    NVOLS = pe.Node(util.Function(input_names=['in_files'],
                                  output_names=['nvols'],
                    function=get_img_nvols),
                    name='NVOLS')

    cp = pe.Node(interface=fsl.ImageMaths(),
                    name='cp')


    delete_first_volume = pe.Node(interface=fsl.ExtractROI(),
                     name='delete_first_volume')
    delete_first_volume.inputs.t_min = 1

    concatnode = pe.Node(interface=util.Merge(2),
                            name='concatnode')

    selectnode = pe.Node(interface=util.Select(),
                            name='selectnode')

    pspec = pe.Node(interface=fsl.PowerSpectrum(),
                       name='pspec')

    ##compute sqrt_pspec of power spectrum
    sqrt_pspec = pe.Node(interface=fsl.ImageMaths(),
                      name='sqrt_pspec')
    sqrt_pspec.inputs.op_string = '-sqrt'

    calculate_low_frequency_point = pe.Node(util.Function(input_names=['nvols',
                                      'TR', 'HP'],
                                      output_names=['n1'],
                        function=get_N1),
                        name='calculate_low_frequency_point')

    calculate_high_frequency_point = pe.Node(util.Function(input_names=['nvols',
                                      'TR', 'LP', 'HP'],
                                      output_names=['n2'],
                        function=get_N2),
                        name='calculate_high_frequency_point')
    cut_low_frequency_data = pe.Node(interface=fsl.ExtractROI(),
                      name='cut_low_frequency_data')

    ## calculate ALFF as the sum_amplitudes_low_frequency of the amplitudes
    ## in the low frequency band
    sum_amplitudes_low_frequency = pe.Node(interface=fsl.ImageMaths(),
                      name='sum_amplitudes_low_frequency')

    ## 4. Calculate fALFF
    amplitude_of_total_frequency = pe.Node(interface=fsl.ImageMaths(),
                       name='amplitude_of_total_frequency')

    fALFF = pe.Node(interface=fsl.MultiImageMaths(),
                        name='fALFF')
    fALFF.inputs.op_string = '-div %s'

    ## 5. Z-normalisation across whole brain
    ALFF_mean = pe.Node(interface=fsl.ImageStats(),
                       name='ALFF_mean')
    ALFF_mean.inputs.op_string = '-k %s -m'

    ALFF_std = pe.Node(interface=fsl.ImageStats(),
                       name='ALFF_std')
    ALFF_std.inputs.op_string = '-k %s -s'

    fALFF_mean = pe.Node(interface=fsl.ImageStats(),
                        name='fALFF_mean')
    fALFF_mean.inputs.op_string = '-k %s -m'

    fALFF_std = pe.Node(interface=fsl.ImageStats(),
                        name='fALFF_std')
    fALFF_std.inputs.op_string = '-k %s -s'

    op_string = pe.Node(util.Function(input_names=['mean',
                                         'std_dev'],
                                         output_names=['op_string'],
                           function=get_operand_string),
                           name='alff_op_string')

    op_string1 = op_string.clone('op_string1')

    alff_Z = pe.Node(interface=fsl.MultiImageMaths(),
                        name='alff_Z')

    falff_Z = pe.Node(interface=fsl.MultiImageMaths(),
                         name='falff_Z')

    alff.connect(inputNode, 'tr',
                 TR, 'TRa')
    alff.connect(inputNode, 'rest_res',
                 TR, 'in_files')
    alff.connect(inputNode, 'rest_res',
                 NVOLS, 'in_files')
    alff.connect(inputNode, 'rest_res',
                 delete_first_volume, 'in_file')
    alff.connect(NVOLS, 'nvols',
                 delete_first_volume, 't_size')
    alff.connect(inputNode, 'rest_res',
                 cp, 'in_file')
    alff.connect(delete_first_volume, 'roi_file',
                 concatnode, 'in1')
    alff.connect(cp, 'out_file',
                 concatnode, 'in2')
    alff.connect(concatnode, 'out',
                 selectnode, 'inlist')
    alff.connect(NVOLS, ('nvols', takemod),
                 selectnode, 'index')
    alff.connect(selectnode, 'out',
                 pspec, 'in_file')
    alff.connect(pspec, 'out_file',
                 sqrt_pspec, 'in_file')

    alff.connect(NVOLS, 'nvols',
                 calculate_low_frequency_point, 'nvols')
    alff.connect(TR, 'TR',
                 calculate_low_frequency_point, 'TR')
    alff.connect(inputnode_hp, 'hp',
                 calculate_low_frequency_point, 'HP')

    alff.connect(NVOLS, 'nvols',
                 calculate_high_frequency_point, 'nvols')
    alff.connect(TR, 'TR',
                 calculate_high_frequency_point, 'TR')
    alff.connect(inputnode_lp, 'lp',
                 calculate_high_frequency_point, 'LP')
    alff.connect(inputnode_hp, 'hp',
                 calculate_high_frequency_point, 'HP')

    alff.connect(sqrt_pspec, 'out_file',
                 cut_low_frequency_data, 'in_file')
    alff.connect(calculate_low_frequency_point, 'n1',
                 cut_low_frequency_data, 't_min')
    alff.connect(calculate_high_frequency_point, 'n2',
                 cut_low_frequency_data, 't_size')
    alff.connect(cut_low_frequency_data, 'roi_file',
                 sum_amplitudes_low_frequency, 'in_file')
    alff.connect(calculate_high_frequency_point, ('n2', set_op_str),
                 sum_amplitudes_low_frequency, 'op_string')

    alff.connect(sqrt_pspec, 'out_file',
                 amplitude_of_total_frequency, 'in_file')
    alff.connect(NVOLS, ('nvols', set_op1_str),
                 amplitude_of_total_frequency, 'op_string')
    alff.connect(sum_amplitudes_low_frequency, 'out_file',
                 fALFF, 'in_file')
    alff.connect(amplitude_of_total_frequency, 'out_file',
                 fALFF, 'operand_files')

    alff.connect(sum_amplitudes_low_frequency, 'out_file',
                 ALFF_mean, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 ALFF_mean, 'mask_file')
    alff.connect(sum_amplitudes_low_frequency, 'out_file',
                 ALFF_std, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 ALFF_std, 'mask_file')
    alff.connect(fALFF, 'out_file',
                 fALFF_mean, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 fALFF_mean, 'mask_file')
    alff.connect(fALFF, 'out_file',
                 fALFF_std, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 fALFF_std, 'mask_file')

    alff.connect(ALFF_mean, 'out_stat',
                 op_string, 'mean')
    alff.connect(ALFF_std, 'out_stat',
                 op_string, 'std_dev')
    alff.connect(op_string, 'op_string',
                 alff_Z, 'op_string')
    alff.connect(sum_amplitudes_low_frequency, 'out_file',
                 alff_Z, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 alff_Z, 'operand_files')

    alff.connect(fALFF_mean, 'out_stat',
                 op_string1, 'mean')
    alff.connect(fALFF_std, 'out_stat',
                 op_string1, 'std_dev')
    alff.connect(op_string1, 'op_string',
                 falff_Z, 'op_string')
    alff.connect(fALFF, 'out_file',
                 falff_Z, 'in_file')
    alff.connect(inputNode, 'rest_mask',
                 falff_Z, 'operand_files')


    alff.connect(pspec, 'out_file',
                 outputNode, 'power_spectrum_distribution')
    alff.connect(sum_amplitudes_low_frequency, 'out_file',
                 outputNode, 'alff_img')
    alff.connect(fALFF, 'out_file',
                 outputNode, 'falff_img')
    alff.connect(alff_Z, 'out_file',
                 outputNode, 'alff_Z_img')
    alff.connect(falff_Z, 'out_file',
                 outputNode, 'falff_Z_img')
    return alff

