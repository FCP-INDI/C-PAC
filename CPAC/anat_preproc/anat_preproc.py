from nipype.interfaces.afni import preprocess
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces import fsl

def create_anat_preproc(use_AFNI, already_skullstripped=False):
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

        inputspec.anat : mprage file or a list of mprage nifti file
            User input anatomical(T1) Image, in any of the 8 orientations

    Workflow Outputs::

        outputspec.refit : nifti file
            Deobliqued anatomical data
        outputspec.reorient : nifti file
            RPI oriented anatomical data
        outputspec.skullstrip : nifti file
            Skull Stripped RPI oriented mprage file with normalized intensities.
        outputspec.brain : nifti file
            Skull Stripped RPI Brain Image with original intensity values and not normalized or scaled.

    Order of commands:

    - Deobliqing the scans.  For details see `3drefit <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3drefit.html>`_::

        3drefit -deoblique mprage.nii.gz

    - Re-orienting the Image into Right-to-Left Posterior-to-Anterior Inferior-to-Superior  (RPI) orientation.  For details see `3dresample <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dresample.html>`_::

        3dresample -orient RPI -prefix mprage_RPI.nii.gz -inset mprage.nii.gz

    - SkullStripping the image.  For details see `3dSkullStrip <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dSkullStrip.html>`_::

        3dSkullStrip -input mprage_RPI.nii.gz -o_ply mprage_RPI_3dT.nii.gz

    - The skull stripping step modifies the intensity values. To get back the original intensity values, we do an element wise product of RPI data with step function of skull Stripped data.  For details see `3dcalc <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_::

        3dcalc -a mprage_RPI.nii.gz -b mprage_RPI_3dT.nii.gz -expr 'a*step(b)' -prefix mprage_RPI_3dc.nii.gz

    High Level Workflow Graph:

    .. image:: ../images/anatpreproc_graph.dot.png
       :width: 500


    Detailed Workflow Graph:

    .. image:: ../images/anatpreproc_graph_detailed.dot.png
       :width: 500

    Examples
    --------

    >>> import anat
    >>> preproc = create_anat_preproc()
    >>> preproc.inputs.inputspec.anat='sub1/anat/mprage.nii.gz'
    >>> preproc.run() #doctest: +SKIP

    """

    preproc = pe.Workflow(name='anat_preproc')
    inputNode = pe.Node(util.IdentityInterface(fields=['anat']),name='inputspec')
    inputNode_AFNI = pe.Node(util.IdentityInterface(fields =['shrink_factor','var_shrink_fac','shrink_factor_bot_lim','avoid_vent','niter','pushout','touchup','fill_hole','avoid_eyes','use_edges','exp_frac','smooth_final','push_to_edge','use_skull','perc_init','max_inter_init','blur_fwhm','fac']),name ='AFNI_options')
    inputNode_BET = pe.Node(util.IdentityInterface(fields=['center','mask_boolean','mesh_boolean','outline','padding','radius','reduce_bias','remove_eyes','robust','skull','surfaces','threshold','vertical_gradient']),name = 'BET_options')
    
    outputNode = pe.Node(util.IdentityInterface(fields=['refit',
                                                        'reorient',
                                                        'skullstrip',
                                                        'brain']),
                         name='outputspec')
    

    try:
        from nipype.interfaces.afni import utils as afni_utils
        anat_deoblique = pe.Node(interface=afni_utils.Refit(),
                                 name='anat_deoblique')
    except ImportError:
        anat_deoblique = pe.Node(interface=preprocess.Refit(),
                                 name='anat_deoblique')

    anat_deoblique.inputs.deoblique = True

    try:
        anat_reorient = pe.Node(interface=afni_utils.Resample(),
                                name='anat_reorient')
    except UnboundLocalError:
        anat_reorient = pe.Node(interface=preprocess.Resample(),
                                name='anat_reorient')

    anat_reorient.inputs.orientation = 'RPI'
    anat_reorient.inputs.outputtype = 'NIFTI_GZ'

    if not already_skullstripped:
        if use_AFNI == True:
    
            anat_skullstrip = pe.Node(interface=preprocess.SkullStrip(),name='anat_skullstrip')
            anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'
            anat_skullstrip.inputs.args = '%(shrink_fac)s %(var_shrink_fac)s %(shrink_fac_bot_lim)s %(avoid_vent)s %(niter)s %(pushout)s %(touchup)s %(fill_hole)s %(NN_smooth)s %(avoid_eyes)s %(use_edge)s %(exp_frac)s %(smooth_final)s %(push_to_edge)s %(use_skull)s %(perc_init)s %(max_inter_iter)s %(blur_fwhm)s %(fac)s' % ({'shrink_fac': "inputNode_AFNI[shrink_factor]", 'var_shrink_fac': "inputNode_AFNI[var_shrink_fac]", 'shrink_fac_bot_lim': "inputNode_AFNI[shrink_factor_bot_lim]", 'avoid_vent' : "inputNode_AFNI[avoid_vent]", 'niter' : "inputNode_AFNI[niter]",'pushout': "inputNode_AFNI[pushout]", 'touchup': "inputNode_AFNI[touchup]", 'fill_hole' : "inputNode_AFNI[fill_hole]", 'avoid_eyes' : "inputNode_AFNI[avoid_eyes]", 'use_edge': "inputNode_AFNI[use_edges]", 'NN_smooth': "inputNode_AFNI[NN_smooth]",'exp_frac': "inputNode_AFNI[exp_frac]", 'smooth_final':"inputNode_AFNI[smooth_final]", 'push_to_edge' : "inputNode_AFNI[push_to_edge]",'use_skull':"inputNode_AFNI[use_skull]", 'perc_init':"inputNode_AFNI[perc_init]", 'max_inter_iter': "inputNode_AFNI[max_inter_init]", 'blur_fwhm': "inputNode_AFNI[blur_fwhm]", 'fac' : "inputNode_AFNI[fac]"})
        else:
            anat_skullstrip = pe.Node(interface=fsl.BET(),name = 'anat_skullstrip')

    try:
        anat_skullstrip_orig_vol = pe.Node(interface=afni_utils.Calc(),
                                           name='anat_skullstrip_orig_vol')
    except UnboundLocalError:
        anat_skullstrip_orig_vol = pe.Node(interface=preprocess.Calc(),
                                           name='anat_skullstrip_orig_vol')

    anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
    anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

    preproc.connect(inputNode, 'anat',anat_deoblique, 'in_file')
    preproc.connect(anat_deoblique,'out_file',anat_reorient,'in_file')
    if not already_skullstripped:
        if use_AFNI == True:
            
            preproc.connect(anat_reorient, 'out_file',anat_skullstrip, 'in_file')
            preproc.connect(anat_skullstrip, 'out_file',anat_skullstrip_orig_vol, 'in_file_b')
        
        
        else:
            preproc.connect(anat_reorient,'out_file',anat_skullstrip,'in_file')
            preproc.connect(inputNode_BET,'center',anat_skullstrip,'center')
            preproc.connect(inputNode_AFNI,'shrink_factor', anat_skullstrip, 'frac')
            preproc.connect(inputNode_BET, 'mask_boolean',anat_skullstrip,'mask')
            preproc.connect(inputNode_BET,'mesh_boolean',anat_skullstrip,'mesh')
            preproc.connect(inputNode_BET,'outline',anat_skullstrip,'outline')
            preproc.connect(inputNode_BET,'padding',anat_skullstrip,'padding')
            preproc.connect(inputNode_BET,'radius',anat_skullstrip,'radius')
            preproc.connect(inputNode_BET,'reduce_bias',anat_skullstrip,'reduce_bias')
            preproc.connect(inputNode_BET,'remove_eyes',anat_skullstrip,'remove_eyes')
            preproc.connect(inputNode_BET,'robust',anat_skullstrip,'robust')
            preproc.connect(inputNode_BET,'skull',anat_skullstrip,'skull')
            preproc.connect(inputNode_BET,'surfaces',anat_skullstrip,'surfaces')
            preproc.connect(inputNode_BET,'threshold',anat_skullstrip,'threshold')
            preproc.connect(inputNode_BET,'vertical_gradient',anat_skullstrip,'vertical_gradient')
            preproc.connect(anat_skullstrip,'out_file',anat_skullstrip_orig_vol,'in_file_b')

    else:
        preproc.connect(anat_reorient, 'out_file',
                        anat_skullstrip_orig_vol, 'in_file_b')
    preproc.connect(anat_reorient, 'out_file',
                    anat_skullstrip_orig_vol, 'in_file_a')

    preproc.connect(anat_deoblique, 'out_file',
                    outputNode, 'refit')
    preproc.connect(anat_reorient, 'out_file',
                    outputNode, 'reorient')
    if not already_skullstripped:
        preproc.connect(anat_skullstrip, 'out_file',
                        outputNode, 'skullstrip')
    preproc.connect(anat_skullstrip_orig_vol, 'out_file',
                    outputNode, 'brain')

    return preproc
