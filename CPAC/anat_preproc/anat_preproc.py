# -*- coding: utf-8 -*-

from nipype.interfaces import afni
from nipype.interfaces import fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.anat_preproc.utils import create_3dskullstrip_arg_string


def create_anat_preproc(method='afni', already_skullstripped=False,
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

    # Anatomical reorientation
    anat_reorient = pe.Node(interface=afni.Resample(),
                            name='anat_reorient')

    anat_reorient.inputs.orientation = 'RPI'
    anat_reorient.inputs.outputtype = 'NIFTI_GZ'
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

            preproc.connect(inputnode_afni, 'shrink_factor',
                            skullstrip_args, 'shrink_fac')
            preproc.connect(inputnode_afni, 'var_shrink_fac',
                            skullstrip_args, 'var_shrink_fac')
            preproc.connect(inputnode_afni, 'shrink_fac_bot_lim',
                            skullstrip_args, 'shrink_fac_bot_lim')
            preproc.connect(inputnode_afni, 'avoid_vent',
                            skullstrip_args, 'avoid_vent')
            preproc.connect(inputnode_afni, 'niter',
                            skullstrip_args, 'niter')
            preproc.connect(inputnode_afni, 'pushout',
                            skullstrip_args, 'pushout')
            preproc.connect(inputnode_afni, 'touchup',
                            skullstrip_args, 'touchup')
            preproc.connect(inputnode_afni, 'fill_hole',
                            skullstrip_args, 'fill_hole')
            preproc.connect(inputnode_afni, 'avoid_eyes',
                            skullstrip_args, 'avoid_eyes')
            preproc.connect(inputnode_afni, 'use_edge',
                            skullstrip_args, 'use_edge')
            preproc.connect(inputnode_afni, 'exp_frac',
                            skullstrip_args, 'exp_frac')
            preproc.connect(inputnode_afni, 'smooth_final',
                            skullstrip_args, 'smooth_final')
            preproc.connect(inputnode_afni, 'push_to_edge',
                            skullstrip_args, 'push_to_edge')
            preproc.connect(inputnode_afni, 'use_skull',
                            skullstrip_args, 'use_skull')
            preproc.connect(inputnode_afni, 'perc_int',
                            skullstrip_args, 'perc_int')
            preproc.connect(inputnode_afni, 'max_inter_iter',
                            skullstrip_args, 'max_inter_iter')
            preproc.connect(inputnode_afni, 'blur_fwhm',
                            skullstrip_args, 'blur_fwhm')
            preproc.connect(inputnode_afni, 'fac',
                            skullstrip_args, 'fac')

            anat_skullstrip = pe.Node(interface=afni.SkullStrip(),
                                      name='anat_skullstrip')

            anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip, 'in_file')
            preproc.connect(skullstrip_args, 'expr',
                            anat_skullstrip, 'args')

            preproc.connect(anat_skullstrip, 'out_file',
                            outputnode, 'skullstrip')

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

            preproc.connect(inputnode_bet, 'frac',
                            anat_skullstrip, 'frac')
            preproc.connect(inputnode_bet, 'mask_boolean',
                            anat_skullstrip, 'mask')
            preproc.connect(inputnode_bet, 'mesh_boolean',
                            anat_skullstrip, 'mesh')
            preproc.connect(inputnode_bet, 'outline',
                            anat_skullstrip, 'outline')
            preproc.connect(inputnode_bet, 'padding',
                            anat_skullstrip, 'padding')
            preproc.connect(inputnode_bet, 'radius',
                            anat_skullstrip, 'radius')
            preproc.connect(inputnode_bet, 'reduce_bias',
                            anat_skullstrip, 'reduce_bias')
            preproc.connect(inputnode_bet, 'remove_eyes',
                            anat_skullstrip, 'remove_eyes')
            preproc.connect(inputnode_bet, 'robust',
                            anat_skullstrip, 'robust')
            preproc.connect(inputnode_bet, 'skull',
                            anat_skullstrip, 'skull')
            preproc.connect(inputnode_bet, 'surfaces',
                            anat_skullstrip, 'surfaces')
            preproc.connect(inputnode_bet, 'threshold',
                            anat_skullstrip, 'threshold')
            preproc.connect(inputnode_bet, 'vertical_gradient',
                            anat_skullstrip, 'vertical_gradient')

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

        preproc.connect(anat_skullstrip_orig_vol, 'out_file',
                        outputnode, 'brain')

    return preproc
