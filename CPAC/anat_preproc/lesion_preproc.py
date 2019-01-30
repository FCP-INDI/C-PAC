# -*- coding: utf-8 -*-

from nipype.interfaces import afni
from nipype.interfaces import fsl
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.anat_preproc.utils import create_3dskullstrip_arg_string


def create_lesion_preproc(wf_name='lesion_preproc'):
    """
    The main purpose of this workflow is to process lesions masks.
    Lesion mask file is deobliqued and reoriented in the same way as the T1 in
    the anat_preproc function.

    Returns
    -------
    lesion_preproc : workflow
        Lesion Preprocessing Workflow

    Workflow Inputs::
        inputspec.lesion : string
            User input lesion mask, in any of the 8 orientations

    Workflow Outputs::

        outputspec.refit : string
            Path to deobliqued anatomical image

        outputspec.reorient : string
            Path to RPI oriented anatomical image

    Order of commands:
    - Deobliqing the scans. ::
        3drefit -deoblique mprage.nii.gz

    - Re-orienting the Image into Right-to-Left Posterior-to-Anterior Inferior-to-Superior  (RPI) orientation ::
        3dresample -orient RPI
                   -prefix mprage_RPI.nii.gz
                   -inset mprage.nii.gz

    Examples
    --------
    >>> from CPAC.anat_preproc.lesion_preproc import create_lesion_preproc
    >>> preproc = create_lesion_preproc()
    >>> preproc.inputs.inputspec.lesion = 'sub1/anat/lesion-mask.nii.gz'
    >>> preproc.run() #doctest: +SKIP
    """

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(
        fields=['lesion']), name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['refit',
                                                        'reorient']),
                         name='outputspec')

    lesion_deoblique = pe.Node(interface=afni.Refit(),
                               name='lesion_deoblique')

    lesion_deoblique.inputs.deoblique = True
    preproc.connect(
        inputnode, 'lesion', lesion_deoblique, 'in_file')
    preproc.connect(
        lesion_deoblique, 'out_file', outputnode, 'refit')

    # Anatomical reorientation
    lesion_reorient = pe.Node(interface=afni.Resample(),
                              name='lesion_reorient')

    lesion_reorient.inputs.orientation = 'RPI'
    lesion_reorient.inputs.outputtype = 'NIFTI_GZ'
    preproc.connect(
        lesion_deoblique, 'out_file', lesion_reorient,
        'in_file')
    preproc.connect(
        lesion_reorient, 'out_file', outputnode, 'reorient')

    return preproc
