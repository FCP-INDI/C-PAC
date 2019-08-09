# -*- coding: utf-8 -*-

from nipype.interfaces import afni
from nipype.interfaces import ants
from nipype.interfaces import fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
# from nipype.interfaces.ants import DenoiseImage
# from nipype.interfaces.ants import N4BiasFieldCorrection
from CPAC.anat_preproc.ants import init_brain_extraction_wf


from CPAC.anat_preproc.utils import create_3dskullstrip_arg_string


def create_anat_preproc(template_path=None, mask_path=None, regmask_path=None, method='afni', already_skullstripped=False,
                        non_local_means_filtering=True, n4_correction=True,
                        wf_name='anat_preproc'):
    """ 
    The main purpose of this workflow is to process T1 scans. Raw mprage file is deobliqued, reoriented
    into RPI and skullstripped. Also, a whole brain only mask is generated from the skull stripped image
    for later use in registration.

    Returns
    -------
    anat_preproc : workflow
        Anatomical Preprocessing Workflow

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/anat_preproc/anat_preproc.py>`_
    
    Workflow Inputs::
        inputspec.anat : string
            User input anatomical (T1) Image, in any of the 8 orientations
    
    Workflow Outputs::

        outputspec.refit : string
            Path to deobliqued anatomical image
    
        outputspec.reorient : string
            Path to RPI oriented anatomical image
    
        outputspec.skullstrip : string
            Path to skull stripped RPI oriented mprage file with normalized intensities.
    
        outputspec.brain : string
            Path to skull stripped RPI brain image with original intensity values and not normalized or scaled.
    
    Order of commands:
    - Deobliqing the scans. ::
        3drefit -deoblique mprage.nii.gz

    - Re-orienting the Image into Right-to-Left Posterior-to-Anterior Inferior-to-Superior  (RPI) orientation ::
        3dresample -orient RPI
                   -prefix mprage_RPI.nii.gz
                   -inset mprage.nii.gz
                   
    - Skull-Stripping the image ::
        Using AFNI ::
            3dSkullStrip -input mprage_RPI.nii.gz
                         -o_ply mprage_RPI_3dT.nii.gz
        or using BET ::
            bet mprage_RPI.nii.gz

    - The skull-stripping step modifies the intensity values. To get back the original intensity values, we do an element wise product of RPI data with step function of skull-stripped data ::
        3dcalc -a mprage_RPI.nii.gz
               -b mprage_RPI_3dT.nii.gz
               -expr 'a*step(b)'
               -prefix mprage_RPI_3dc.nii.gz
    
    High Level Workflow Graph:
    .. image:: ../images/anatpreproc_graph.dot.png
       :width: 500
    
    Detailed Workflow Graph:
    .. image:: ../images/anatpreproc_graph_detailed.dot.png
       :width: 500

    Examples
    --------
    >>> from CPAC.anat_preproc import create_anat_preproc
    >>> preproc = create_anat_preproc()
    >>> preproc.inputs.inputspec.anat = 'sub1/anat/mprage.nii.gz'
    >>> preproc.run() #doctest: +SKIP
    """

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(
        fields=['anat', 'brain_mask']), name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['refit',
                                                        'reorient',
                                                        'skullstrip',
                                                        'brain']),
                         name='outputspec')

    anat_deoblique = pe.Node(interface=afni.Refit(),
                             name='anat_deoblique')

    anat_deoblique.inputs.deoblique = True
    preproc.connect(inputnode, 'anat', anat_deoblique, 'in_file')
    preproc.connect(anat_deoblique, 'out_file', outputnode, 'refit')

    if non_local_means_filtering and n4_correction:
        denoise = pe.Node(interface = ants.DenoiseImage(), name = 'anat_denoise')
        preproc.connect(anat_deoblique, 'out_file', denoise, 'input_image')
        n4 = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, shrink_factor=2, copy_header=True),
            name='anat_n4')
        preproc.connect(denoise, 'output_image', n4, 'input_image')
    elif non_local_means_filtering and not n4_correction:
        denoise = pe.Node(interface = ants.DenoiseImage(), name = 'anat_denoise')
        preproc.connect(anat_deoblique, 'out_file', denoise, 'input_image')
    elif not non_local_means_filtering and n4_correction:
        n4 = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, shrink_factor=2, copy_header=True),
            name='anat_n4')
        preproc.connect(anat_deoblique, 'out_file', n4, 'input_image')

    # Anatomical reorientation
    anat_reorient = pe.Node(interface=afni.Resample(),
                            name='anat_reorient')

    anat_reorient.inputs.orientation = 'RPI'
    anat_reorient.inputs.outputtype = 'NIFTI_GZ'
    
    if n4_correction:
        preproc.connect(n4, 'output_image', anat_reorient, 'in_file')
    elif non_local_means_filtering and not n4_correction:
        preproc.connect(denoise, 'output_image', anat_reorient, 'in_file')
    else:
        preproc.connect(anat_deoblique, 'out_file', anat_reorient, 'in_file')

    preproc.connect(anat_reorient, 'out_file', outputnode, 'reorient')

    if already_skullstripped:

        anat_skullstrip = pe.Node(interface=util.IdentityInterface(fields=['out_file']),
                                    name='anat_skullstrip')

        preproc.connect(anat_reorient, 'out_file',
                        anat_skullstrip, 'out_file')

        preproc.connect(anat_skullstrip, 'out_file',
                        outputnode, 'skullstrip')

        preproc.connect(anat_skullstrip, 'out_file',
                        outputnode, 'brain')

    else:

        if method == 'afni':
            # Skull-stripping using AFNI 3dSkullStrip
            inputnode_afni = pe.Node(
                util.IdentityInterface(fields=['shrink_factor',
                                               'var_shrink_fac',
                                               'shrink_fac_bot_lim',
                                               'avoid_vent',
                                               'niter',
                                               'pushout',
                                               'touchup',
                                               'fill_hole',
                                               'avoid_eyes',
                                               'use_edge',
                                               'exp_frac',
                                               'smooth_final',
                                               'push_to_edge',
                                               'use_skull',
                                               'perc_int',
                                               'max_inter_iter',
                                               'blur_fwhm',
                                               'fac']),
                name='AFNI_options')

            skullstrip_args = pe.Node(util.Function(input_names=['spat_norm',
                                                                 'spat_norm_dxyz',
                                                                 'shrink_fac',
                                                                 'var_shrink_fac',
                                                                 'shrink_fac_bot_lim',
                                                                 'avoid_vent',
                                                                 'niter',
                                                                 'pushout',
                                                                 'touchup',
                                                                 'fill_hole',
                                                                 'avoid_eyes',
                                                                 'use_edge',
                                                                 'exp_frac',
                                                                 'smooth_final',
                                                                 'push_to_edge',
                                                                 'use_skull',
                                                                 'perc_int',
                                                                 'max_inter_iter',
                                                                 'blur_fwhm',
                                                                 'fac'],
                                                    output_names=['expr'],
                                                    function=create_3dskullstrip_arg_string),
                                      name='anat_skullstrip_args')

            preproc.connect([
                (inputnode_afni, skullstrip_args, [
                    ('shrink_factor', 'shrink_fac'),
                    ('var_shrink_fac', 'var_shrink_fac'),
                    ('shrink_fac_bot_lim', 'shrink_fac_bot_lim'),
                    ('avoid_vent', 'avoid_vent'),
                    ('niter', 'niter'),
                    ('pushout', 'pushout'),
                    ('touchup', 'touchup'),
                    ('fill_hole', 'fill_hole'),
                    ('avoid_eyes', 'avoid_eyes'),
                    ('use_edge', 'use_edge'),
                    ('exp_frac', 'exp_frac'),
                    ('smooth_final', 'smooth_final'),
                    ('push_to_edge', 'push_to_edge'),
                    ('use_skull', 'use_skull'),
                    ('perc_int', 'perc_int'),
                    ('max_inter_iter', 'max_inter_iter'),
                    ('blur_fwhm', 'blur_fwhm'),
                    ('fac', 'fac'),
                ])
            ])

            anat_skullstrip = pe.Node(interface=afni.SkullStrip(),
                                      name='anat_skullstrip')

            anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip, 'in_file')
            preproc.connect(skullstrip_args, 'expr',
                            anat_skullstrip, 'args')

            preproc.connect(anat_skullstrip, 'out_file',
                            outputnode, 'skullstrip')
            # Apply skull-stripping step mask to original volume
            anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                            name='anat_skullstrip_orig_vol')

            anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
            anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip_orig_vol, 'in_file_a')
            
            if method == 'mask':
                preproc.connect(inputnode, 'brain_mask',
                                anat_skullstrip_orig_vol, 'in_file_b')
            else:
                preproc.connect(anat_skullstrip, 'out_file',
                                anat_skullstrip_orig_vol, 'in_file_b')

        elif method == 'fsl':
            inputnode_bet = pe.Node(
                util.IdentityInterface(fields=['frac',
                                               'mask_boolean',
                                               'mesh_boolean',
                                               'outline',
                                               'padding',
                                               'radius',
                                               'reduce_bias',
                                               'remove_eyes',
                                               'robust',
                                               'skull',
                                               'surfaces',
                                               'threshold',
                                               'vertical_gradient']),
                name='BET_options')
            # Skull-stripping using FSL BET
            anat_skullstrip = pe.Node(
                interface=fsl.BET(), name='anat_skullstrip')

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip, 'in_file')

            preproc.connect([
                (inputnode_bet, anat_skullstrip, [
                    ('frac', 'frac'),
                    ('mask_boolean', 'mask'),
                    ('mesh_boolean', 'mesh'),
                    ('outline', 'outline'),
                    ('padding', 'padding'),
                    ('radius', 'radius'),
                    ('reduce_bias', 'reduce_bias'),
                    ('remove_eyes', 'remove_eyes'),
                    ('robust', 'robust'),
                    ('skull', 'skull'),
                    ('surfaces', 'surfaces'),
                    ('threshold', 'threshold'),
                    ('vertical_gradient', 'vertical_gradient'),
                ])
            ])

            preproc.connect(anat_skullstrip, 'out_file',
                            outputnode, 'skullstrip')
            
                # Apply skull-stripping step mask to original volume
            anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                            name='anat_skullstrip_orig_vol')

            anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
            anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip_orig_vol, 'in_file_a')  
                            
            if method == 'mask':
                preproc.connect(inputnode, 'brain_mask',
                                anat_skullstrip_orig_vol, 'in_file_b')
            else:
                preproc.connect(anat_skullstrip, 'out_file',
                                anat_skullstrip_orig_vol, 'in_file_b')

        elif method == 'antsBrainExtraction':   
            anat_skullstrip_ants = init_brain_extraction_wf(tpl_target_path=template_path,
                                                            tpl_mask_path=mask_path,
                                                            tpl_regmask_path=regmask_path,
                                                            name='antsBrainExtraction_options')
                                                            # use_float=True,
                                                            # normalization_quality='precise',
                                                            # omp_nthreads=None,
                                                            # mem_gb=3.0,
                                                            # bids_suffix='T1w',
                                                            # atropos_refine=True,
                                                            # atropos_use_random_seed=True,
                                                            # atropos_model=None,
                                                            # use_laplacian=True,
                                                            # bspline_fitting_distance=200)
            

            # if non_local_means_filtering or n4_correction:
             
            anat_skullstrip_ants.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_deoblique, 'out_file', 
                            anat_reorient, 'in_file')
            
            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip_ants, 'in_files')

            preproc.connect(anat_skullstrip_ants, 'out_file',
                            outputnode, 'skullstrip')
                            
            anat_skullstrip_orig_vol = anat_skullstrip_ants.out_file

            ### if provide mask, antsBrainExtraction how to deal with it?
        # if method == 'mask':
        #     preproc.connect(inputnode, 'brain_mask',
        #                     anat_skullstrip_orig_vol, 'in_file')
        # else:
        #     preproc.connect(anat_skullstrip, 'out_file',
        #                     anat_skullstrip_orig_vol, 'in_file')


        preproc.connect(anat_skullstrip_orig_vol, 'out_file',
                        outputnode, 'brain')

    return preproc
