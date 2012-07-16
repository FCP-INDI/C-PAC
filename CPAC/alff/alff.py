import os
import sys
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from CPAC.utils.utils import *

def create_alff(tr):

    """
    Calculate Amplitude of low frequency oscillations(ALFF) and fractional ALFF maps

    Parameters
    ----------

    tr : Temporal Resolution of the functional alff

    Returns
    -------

    alff_workflow : workflow object

        ALFF workflow

    Examples
    --------


    Notes
    -----

    `Source <https://github.com/openconnectome/C-PAC/blob/master/CPAC/alff/alff.py>`_

    Workflow Inputs: ::

        hp_input.hp : (A list of floating point numbers)
            high pass frequencies

        lp_input.lp : (A list of floating point numbers)
            low pass frequencies

        fwhm_input.fwhm : (A list of floating point numbers)
            full width half max for spatial alff_Z_to_standard_FWHMing

        inputspec.rest_res : (an existing nifti file)
            Nuisance signal regressed functional image

        inputspec.rest_mask : (an existing nifti file)
            A mask volume(derived by dilating the motion corrected functional volume) in native space

        inputspec.rest_mask2standard : (an existing nifti file)
            A mask volume(derived from the functional volume) in standard in standard space
            Used in spatial alff_Z_to_standard_FWHMing the Z-transformed correlations in MNI space

        inputspec.premat : (an existing affine transformation .mat file)
            Specifies an affine transform that should be applied to the data prior to the non-linear warping(example_func2highres.mat).

        inputspec.standard : (an existing nifti file)
            FSL standard nifti file in user specified resolution

        inputspec.fieldcoeff_file : (an existing nifti file)
            File with warp coefficients/fields. This typically the output given by the --cout parameter of fnirt during registration step


    Workflow Outputs: ::

        outputspec.power_spectrum_distribution : (a nifti file)
            outputs image containing the spectral power density of residual functional image

        outputspec.alff_img : (a nifti file)
            outputs image containing the sum of the amplitudes in the low frequency band

        outputspec.falff_img : (a nifti file)
            outputs image containing the sum of the amplitudes in the low frequency band divided by the amplitude of the total frequency

        outputspec.alff_Z_img : (a nifti file)
            outputs image containing Normalized ALFF Z scores across full brain in native space

        outputspec.falff_Z_img : (a nifti file)
            outputs image containing Normalized fALFF Z scores across full brain in native space

        outputspec.alff_Z_2standard_img : (a nifti file)
            outputs image containing normalized  ALFF Z scores in MNI space

        outputspec.falff_Z_2standard_img : (a nifti file)
            outputs image containing normalized  fALFF Z scores in MNI space

        outputspec.alff_Z_2standard_fwhm_img : (a nifti file)
            outputs image containing normalized ALFF Z scores in MNI space with spatial alff_Z_to_standard_FWHMing applied to them

        outputspec.falff_Z_2standard_fwhm_img : (a nifti file)
            outputs image containing normalized fALFF Z scores in MNI space with spatial alff_Z_to_standard_FWHMing applied to them


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

    - Register Z-transformed ALFF to standard space ::

        applywarp
        --ref = MNI152_T1_STANDARD_RES.nii.gz
        --in = ALFF_Z.nii.gz
        --out = ALFF_Z_2standard.nii.gz
        --warp = highres2standard_warp.nii.gz
        --premat = example_func2highres.mat

    - Register Z-transformed fALFF to standard space ::

        applywarp
        --ref = MNI152_T1_STANDARD_RES.nii.gz
        --in = fALFF_Z.nii.gz
        --out = fALFF_Z_2standard.nii.gz
        --warp = highres2standard_warp.nii.gz
        --premat = example_func2highres.mat

    - Spatially Smooth the ALFF results ::

        fslmaths
        ALFF_Z_2standard.nii.gz
        -kernel gauss FWHM/ sqrt(8*ln(2))
        -fmean
        -mas rest_mask2standard.nii.gz
        ALFF_Z_2standard_FWHM.nii.gz

    - Spatially Smooth the f/ALFF results ::

        fslmaths
        fALFF_Z_2standard.nii.gz
        -kernel gauss FWHM/ sqrt(8*ln(2))
        -fmean
        -mas rest_mask2standard.nii.gz
        fALFF_Z_2standard_FWHM.nii.gz

    Workflow ALFF and fractional ALFF:
    
    .. image:: ../images/alff_graph.dot.png
        :width: 500

    
    References
    ----------

    .. [1] Zou, Q.-H., Zhu, C.-Z., Yang, Y., Zuo, X.-N., Long, X.-Y., Cao, Q.-J., Wang, Y.-F., et al. (2008). An improved approach to detection of amplitude of low-frequency fluctuation (ALFF) for resting-state fMRI: fractional ALFF. Journal of neuroscience methods, 172(1), 137-41. doi:10.10

    Examples
    --------

    >>> alff_w = create_alff(tr=2.0)
    >>> alff_w.inputs.fwhm_input.fwhm = [4.5, 6]
    >>> alff_w.get_node('fwhm_input').iterables = ('fwhm', [4.5, 6])
    >>> alff_w.inputs.inputspec.standard = '/usr/share/fsl/4.1/data/standard/MNI_152_T1_2mm_nii.gz'
    >>> alff_w.inputs.hp_input.hp = [0.01]
    >>> alff_w.inputs.lp_input.lp = [0.1]
    >>> alff_w.get_node('hp_input').iterables = ('hp',
                                                [0.01])
    >>> alff_w.get_node('lp_input').iterables = ('lp',
                                                [0.1])
    >>> alff_w.inputs.inputspec.premat = '/home/data/subject/func/example_func2highres.mat'
    >>> alff_w.inputs.inputspec.rest_res = '/home/data/subject/func/rest_bandpassed.nii.gz'
    >>> alff_w.inputs.inputspec.fieldcoeff_file = '/home/data/subject/func/highres2standard_warp.nii.gz'
    >>> alff_w.inputs.inputspec.rest_mask2standard = '/home/data/subject/func/rest_mask2standard.nii.gz' 
    >>> alff_w.run() # doctest: +SKIP


    """


    alff = pe.Workflow(name='alff_workflow')
    inputNode = pe.Node(util.IdentityInterface(fields=['rest_res',
                                                'rest_mask',
                                                'rest_mask2standard',
                                                'premat',
                                                'standard',
                                                'fieldcoeff_file',
                                                    ]),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=[
                                            'power_spectrum_distribution',
                                            'alff_img',
                                            'falff_img',
                                            'alff_Z_img',
                                            'falff_Z_img',
                                            'alff_Z_2standard_img',
                                            'falff_Z_2standard_img',
                                            'alff_Z_2standard_fwhm_img',
                                            'falff_Z_2standard_fwhm_img']),
                          name='outputspec')



    inputnode_hp = pe.Node(util.IdentityInterface(fields=['hp']),
                             name='hp_input')

    inputnode_lp = pe.Node(util.IdentityInterface(fields=['lp']),
                             name='lp_input')

    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')

    TR = pe.Node(util.Function(input_names=['in_files', 'TRa'],
                               output_names=['TR'],
                 function=getImgTR), name='TR')
    TR.inputs.TRa = tr

    NVOLS = pe.Node(util.Function(input_names=['in_files'],
                                  output_names=['nvols'],
                    function=getImgNVols),
                    name='NVOLS')



    delete_first_volume = pe.MapNode(interface=fsl.ExtractROI(),
                     name='delete_first_volume',
                     iterfield=['in_file',
                     't_size'])
    delete_first_volume.inputs.t_min = 1

    concatnode = pe.MapNode(interface=util.Merge(2),
                            name='concatnode',
                            iterfield=['in1', 'in2'])

    selectnode = pe.MapNode(interface=util.Select(),
                            name='selectnode',
                            iterfield=['inlist', 'index'])

    pspec = pe.MapNode(interface=fsl.PowerSpectrum(),
                       name='pspec',
                       iterfield=['in_file'])

    ##compute sqrt of power spectrum
    sqrt_pspec = pe.MapNode(interface=fsl.ImageMaths(),
                      name='sqrt_pspec',
                      iterfield=['in_file'])
    sqrt_pspec.inputs.op_string = '-sqrt'

    calculate_low_frequency_point = pe.MapNode(util.Function(input_names=['nvols',
                                      'TR', 'HP'],
                                      output_names=['n1'],
                        function=getN1),
                        name='calculate_low_frequency_point',
                        iterfield=['nvols', 'TR'])

    calculate_high_frequency_point = pe.MapNode(util.Function(input_names=['nvols',
                                      'TR', 'LP', 'HP'],
                                      output_names=['n2'],
                        function=getN2),
                        name='calculate_high_frequency_point',
                        iterfield=['nvols', 'TR'])
    cut_low_frequency_data = pe.MapNode(interface=fsl.ExtractROI(),
                      name='cut_low_frequency_data',
                      iterfield=['in_file',
                                 't_min', 't_size'])

    ## calculate ALFF as the sum_amplitudes_low_frequency of the amplitudes
    ## in the low frequency band
    sum_amplitudes_low_frequency = pe.MapNode(interface=fsl.ImageMaths(),
                      name='sum_amplitudes_low_frequency',
                      iterfield=['in_file',
                      'op_string'])

    ## 4. Calculate fALFF
    amplitude_of_total_frequency = pe.MapNode(interface=fsl.ImageMaths(),
                       name='amplitude_of_total_frequency',
                       iterfield=['in_file',
                        'op_string'])

    fALFF = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='fALFF',
                        iterfield=['in_file',
                        'operand_files'])
    fALFF.inputs.op_string = '-div %s'

    ## 5. Z-normalisation across whole brain
    ALFF_mean = pe.MapNode(interface=fsl.ImageStats(),
                       name='ALFF_mean',
                       iterfield=['in_file',
                        'mask_file'])
    ALFF_mean.inputs.op_string = '-k %s -m'

    ALFF_std = pe.MapNode(interface=fsl.ImageStats(),
                       name='ALFF_std',
                       iterfield=['in_file',
                       'mask_file'])
    ALFF_std.inputs.op_string = '-k %s -s'

    fALFF_mean = pe.MapNode(interface=fsl.ImageStats(),
                        name='fALFF_mean',
                        iterfield=['in_file',
                        'mask_file'])
    fALFF_mean.inputs.op_string = '-k %s -m'

    fALFF_std = pe.MapNode(interface=fsl.ImageStats(),
                        name='fALFF_std',
                        iterfield=['in_file',
                        'mask_file'])
    fALFF_std.inputs.op_string = '-k %s -s'

    op_string = pe.MapNode(util.Function(input_names=['mean',
                                         'std_dev'],
                                         output_names=['op_string'],
                           function=getOpString),
                           name='alff_op_string',
                           iterfield=['mean',
                           'std_dev'])

    op_string1 = op_string.clone('op_string1')

    alff_Z = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='alff_Z',
                        iterfield=['in_file',
                        'operand_files',
                        'op_string'])

    falff_Z = pe.MapNode(interface=fsl.MultiImageMaths(),
                         name='falff_Z',
                         iterfield=['in_file',
                            'operand_files',
                            'op_string'])

    #Registering Z-transformed ALFF to standard space
    register_alff_Z_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                           name='register_alff_Z_to_standard',
                           iterfield=['in_file',
                           'premat'])

    register_falff_Z_to_standard = pe.MapNode(interface=fsl.ApplyWarp(),
                            name='register_falff_Z_to_standard',
                            iterfield=['in_file',
                            'premat'])

    alff_Z_to_standard_FWHM = pe.MapNode(interface=fsl.MultiImageMaths(),
                        name='alff_Z_to_standard_FWHM',
                        iterfield=['in_file',
                        'operand_files'])

    falff_Z_to_standard_FWHM = pe.MapNode(interface=fsl.MultiImageMaths(),
                         name='falff_Z_to_standard_FWHM',
                         iterfield=['in_file',
                         'operand_files'])

    alff.connect(inputNode, 'rest_res',
                 TR, 'in_files')
    alff.connect(inputNode, 'rest_res',
                 NVOLS, 'in_files')
    alff.connect(inputNode, 'rest_res',
                 delete_first_volume, 'in_file')
    alff.connect(NVOLS, 'nvols',
                 delete_first_volume, 't_size')
    alff.connect(delete_first_volume, 'roi_file',
                 concatnode, 'in1')
    alff.connect(inputNode, 'rest_res',
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

    alff.connect(inputNode, 'standard',
                 register_alff_Z_to_standard, 'ref_file')
    alff.connect(alff_Z, 'out_file',
                 register_alff_Z_to_standard, 'in_file')
    alff.connect(inputNode, 'fieldcoeff_file',
                 register_alff_Z_to_standard, 'field_file')
    alff.connect(inputNode, 'premat',
                 register_alff_Z_to_standard, 'premat')

    alff.connect(inputNode, 'standard',
                 register_falff_Z_to_standard, 'ref_file')
    alff.connect(falff_Z, 'out_file',
                 register_falff_Z_to_standard, 'in_file')
    alff.connect(inputNode, 'fieldcoeff_file',
                 register_falff_Z_to_standard, 'field_file')
    alff.connect(inputNode, 'premat',
                 register_falff_Z_to_standard, 'premat')

    alff.connect(register_alff_Z_to_standard, 'out_file',
                 alff_Z_to_standard_FWHM, 'in_file')
    alff.connect(inputnode_fwhm, ('fwhm', set_gauss),
                 alff_Z_to_standard_FWHM, 'op_string')
    alff.connect(inputNode, 'rest_mask2standard',
                 alff_Z_to_standard_FWHM, 'operand_files')

    alff.connect(register_falff_Z_to_standard, 'out_file',
                 falff_Z_to_standard_FWHM, 'in_file')
    alff.connect(inputnode_fwhm, ('fwhm', set_gauss),
                 falff_Z_to_standard_FWHM, 'op_string')
    alff.connect(inputNode, 'rest_mask2standard',
                 falff_Z_to_standard_FWHM, 'operand_files')

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
    alff.connect(register_alff_Z_to_standard, 'out_file',
                 outputNode, 'alff_Z_2standard_img')
    alff.connect(register_falff_Z_to_standard, 'out_file',
                 outputNode, 'falff_Z_2standard_img')
    alff.connect(alff_Z_to_standard_FWHM, 'out_file',
                 outputNode, 'alff_Z_2standard_fwhm_img')
    alff.connect(falff_Z_to_standard_FWHM, 'out_file',
                 outputNode, 'falff_Z_2standard_fwhm_img')
    return alff

