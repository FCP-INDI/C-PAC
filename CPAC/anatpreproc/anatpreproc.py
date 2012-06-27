#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import sys
#import e_afni
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
#from utils import *


def create_anat_preproc():
    """ 
    
    The main purpose of this workflow is to process T1 scans. Raw mprage file is deobliqued, reoriented 
    into RPI and skullstripped. Also, a whole brain only mask is generated from the skull stripped image 
    for later use in registration.
    
    Returns 
    -------
    anat_preproc : workflow
        Anatomical preprocessing workflow
    
    Notes
    -----
    
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
    
    Source code: `anat_preproc <https://github.com/ssikka/NKI_NYU_Nipype/blob/development/base.py#L250>`_
    
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
    >>> prproc = create_anat_preproc()
    >>> preproc.inputs.inputspec.anat='sub1/anat/mprage.nii.gz'
    >>> preporc.run() #doctest: +SKIP
            
    """
    preproc = pe.Workflow(name='anatpreproc')
    
    inputNode = pe.Node(util.IdentityInterface(fields=['anat']),
                        name='inputspec')
    
    outputNode = pe.Node(util.IdentityInterface(fields=['refit',
                                                    'reorient',
                                                    'skullstrip',
                                                    'brain']),
                         name='outputspec')
    
    anat_deoblique = pe.Node(interface=e_afni.Threedrefit(),
                         name='anat_deoblique')
    anat_deoblique.inputs.deoblique = True
    
    anat_reorient = pe.Node(interface=afni.Resample(),
                            name='anat_reorient')
    anat_reorient.inputs.orientation = 'RPI'
    
    anat_skullstrip = pe.Node(interface=afni.SkullStrip(),
                              name='anat_skullstrip')
    anat_skullstrip.inputs.options = '-o_ply'
    
    anat_brain_only = pe.Node(interface=e_afni.Threedcalc(),
                        name='anat_brain_only')
    anat_brain_only.inputs.expr = '\'a*step(b)\''
    
    preproc.connect(inputNode, 'anat',
                    anat_deoblique, 'in_file')
    preproc.connect(anat_deoblique, 'out_file',
                    anat_reorient, 'in_file')
    preproc.connect(anat_reorient, 'out_file',
                    anat_skullstrip, 'in_file')
    preproc.connect(anat_skullstrip, 'out_file',
                    anat_brain_only, 'infile_b')
    preproc.connect(anat_reorient, 'out_file',
                    anat_brain_only, 'infile_a')
    
    preproc.connect(anat_deoblique, 'out_file',
                    outputNode, 'refit')
    preproc.connect(anat_reorient, 'out_file',
                    outputNode, 'reorient')
    preproc.connect(anat_skullstrip, 'out_file',
                    outputNode, 'skullstrip')
    preproc.connect(anat_brain_only, 'out_file',
                    outputNode, 'brain')
    
    return preproc
