import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, \
    BaseInterfaceInputSpec, traits, File, TraitedSpec
from nipype.interfaces import fsl

from nilearn import input_data, masking, image, datasets
from nilearn.image import resample_to_img, concat_imgs
from nilearn.input_data import NiftiMasker, NiftiLabelsMasker

from CPAC.utils.function import Function

import os
import copy
import numpy as np
import nibabel as nb



def create_randomise(name='randomise', working_dir=None, crash_dir=None):
    """
    
    Parameters
    ----------
        
    Returns
    -------
    workflow : nipype.pipeline.engine.Workflow
        Randomise workflow.
        
    Notes
    -----
    
    Workflow Inputs::
    
        
    Workflow Outputs::

    
    References
    ----------
    
    """

    if not working_dir:
        working_dir = os.path.join(os.getcwd(), 'Randomise_work_dir')
    if not crash_dir:
        crash_dir = os.path.join(os.getcwd(), 'Randomise_crash_dir')

    wf = pe.Workflow(name=name)
    wf.base_dir = working_dir
    wf.config['execution'] = {'hash_method': 'timestamp',
                              'crashdump_dir': os.path.abspath(crash_dir)}

    inputspec = pe.Node(
        util.IdentityInterface(fields=[
            'subjects',
            'design_matrix_file',
            'constrast_file',
            'permutations',
        ]),
        name='inputspec'
    )

    outputspec = pe.Node(
        util.IdentityInterface(fields=[
            'correlations',
            'significance'
        ]),
        name='outputspec'
    )

    merge = pe.Node(interface=fsl.Merge(), name='fsl_merge')
    merge.inputs.dimension = 't'
    merge.inputs.merged_file = "randomise_merged.nii.gz"

    wf.connect(inputspec, 'subjects', merge, 'in_files')

    mask = pe.Node(interface=fsl.maths.MathsCommand(), name='fsl_maths')
    mask.inputs.args = '-abs -Tmin -bin'
    mask.inputs.out_file = "randomise_mask.nii.gz"

    wf.connect(merge, 'merged_file', mask, 'in_file')

    randomise = pe.Node(interface=fsl.Randomise(), name='randomise')
    wf.connect(mask, 'out_file', randomise, 'mask')
    randomise.inputs.base_name = "randomise"
    randomise.inputs.demean = True
    randomise.inputs.tfce = True
    wf.connect([
        (merge, randomise, [('merged_file', 'in_file')]),
        (inputspec, randomise, [
            ('design_matrix_file', 'design_mat'),
            ('constrast_file', 'tcon'),
            ('permutations', 'num_perm'),
        ]),
    ])

    select_t_corrected = pe.Node(niu.Function(input_names=["input_list"],
                                                output_names=['out_file'],
                                                function=select),
                                    name='select_t_cor{0}'.format(current_contrast))

    wf.connect(randomise, "t_corrected_p_files", select_t_corrected, "input_list")