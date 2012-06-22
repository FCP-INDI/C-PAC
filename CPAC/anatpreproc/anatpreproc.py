#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import sys
import e_afni
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from utils import *


def create_anat_preproc():
    """ 
        
    Anatomical Preprocessing
    ===============================
    
    The main purpose of this workflow is to process T1 scans. Raw mprage file is deobliqued, reoriented 
    into RPI and skullstripped. Also, a whole brain only mask is generated from the skull stripped image 
    for later use in registration.
    
    Source code: `anat_preproc <https://github.com/ssikka/NKI_NYU_Nipype/blob/development/base.py#L250>`_
    
    **Inputs**:
    -----------
    
    - *inputspec.anat* : (mprage file or a list of mprage nifti file) 
        User input anatomical(T1) Image, in any of the 8 orientations
    
    **Outputs**:
    ------------
    
    - *outputspec.refit* : (a nifti file) 
        Deobliqued anatomical data 
    - *outputspec.reorient* : (a nifti file) 
        RPI oriented anatomical data 
    - *outputspec.skullstrip* : (a nifti file)
        Skull Stripped RPI oriented mprage file with normalized intensities.
    - *outputspec.brain* : (a nifti file)
        Skull Stripped RPI Brain Image with original intensity values and not normalized or scaled.
    
    Example
    -------
    
    >>> import anat
    >>> prproc = create_anat_preproc()
    >>> preproc.inputs.inputspec.anat='sub1/anat/mprage.nii.gz'
    >>> preporc.run() #doctest: +SKIP
    
    
    **Commands in Order of Execution**
    ----------------------------------
    
    - Deobliqing the scans
    .. code-block:: python
    
        3drefit -deoblique mprage.nii.gz
        
    For information on Command and the options used, please refer to : `3drefit <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3drefit.html>`_
    
    - Re-orienting the Image into Right-to-Left Posterior-to-Anterior Inferior-to-Superior  (RPI) orientation
    .. code-block:: python
    
        3dresample -orient RPI -prefix mprage_RPI.nii.gz -inset mprage.nii.gz 
    
    For information on Command and the options used, please refer to : `3dresample <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dresample.html>`_ 
    
    - SkullStripping the Image
    .. code-block:: python
        
        3dSkullStrip -input mprage_RPI.nii.gz -o_ply mprage_RPI_3dT.nii.gz
    
    For information on Command and the options used, please refer to : `3dSkullStrip <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dSkullStrip.html>`_
    
    - The skull strippping step modifies the intensity values. To get back the original intensity values,
      we do an element wise product of  RPI data with step function of skull Stripped data.
    .. code-block:: python
        
        3dcalc -a mprage_RPI.nii.gz -b mprage_RPI_3dT.nii.gz -expr 'a*step(b)' -prefix mprage_RPI_3dc.nii.gz
     
    For information on Command and the options used, please refer to : `3dcalc <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_
    
    High Level Workflow Graph:
    --------------------------
    
    .. image:: anatpreproc_graph.dot.png
       :width: 500
       :target: ../images/anatpreproc_graph.dot.png
    
    
    Detailed Workflow Graph:
    ------------------------
    
    .. image:: anatpreproc_graph_detailed.dot.png
       :width: 500
       :target: ../images/anatpreproc_graph_detailed.dot.png
            
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
