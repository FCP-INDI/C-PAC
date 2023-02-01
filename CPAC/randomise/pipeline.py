from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, traits, File, TraitedSpec
from nipype.interfaces import fsl

from nilearn import input_data, masking, image, datasets
from nilearn.image import resample_to_img, concat_imgs
from nilearn.input_data import NiftiMasker, NiftiLabelsMasker

from CPAC.utils.interfaces.function import Function

import os
import copy
import numpy as np
import nibabel as nb



def create_randomise(name='randomise',working_dir=None,crash_dir=None):

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

    inputspec = pe.Node(util.IdentityInterface(fields=['subjects_list','pipeline_output_folder','permutations','mask_boolean','demean','c_thresh']),name='inputspec')
    
    outputspec = pe.Node(util.IdentityInterface(fields=['tstat_files' ,'t_corrected_p_files','index_file','threshold_file','localmax_txt_file','localmax_vol_file','max_file','mean_file','pval_file','size_file']), name='outputspec')
    
    
    #merge = pe.Node(interface=fslMerge(), name='fsl_merge')
    #merge.inputs.dimension = 't'
    #merge.inputs.merged_file = "randomise_merged.nii.gz"

    #wf.connect(inputspec, 'subjects', merge, 'in_files')

    #mask = pe.Node(interface=fsl.maths.MathsCommand(), name='fsl_maths')
    #mask.inputs.args = '-abs -Tmin -bin'
    #mask.inputs.out_file = "randomise_mask.nii.gz"
    #wf.connect(inputspec, 'subjects', mask, 'in_file')

    randomise = pe.Node(interface=fsl.Randomise(), name='randomise')
    randomise.inputs.base_name = "randomise"
    randomise.inputs.demean = True
    randomise.inputs.tfce = True
    wf.connect([(inputspec, randomise, [('subjects', 'in_file'),
            ('design_matrix_file', 'design_mat'),
            ('constrast_file', 'tcon'),
            ('permutations', 'num_perm'),
            ])])
    wf.connect(randomise,'tstat_files',outputspec,'tstat_files')
    wf.connect(randomise,'t_corrected_p_files',outputspec,'t_corrected_p_files')
#------------- issue here arises while using tfce. By not using tfce, you don't get t_corrected_p files. R V in a conundrum? --------------------#

    
    select_tcorrp_files = pe.Node(Function(input_names=['input_list'],output_names=['out_file'],function=select),name='select_t_corrp')

    wf.connect(randomise, 't_corrected_p_files',select_tcorrp_files, 'input_list')
    wf.connect(select_tcorrp_files,'out_file',outputspec,'out_tcorr_corrected')

    select_tstat_files = pe.Node(Function(input_names=['input_list'],output_names=['out_file'],function=select),name='select_t_stat')

    wf.connect(randomise, 'tstat_files',select_tstat_files, 'input_list')
    wf.connect(select_tstat_files,'out_file',outputspec,'out_tstat_corrected')

    thresh = pe.Node(interface=fsl.Threshold(),name='fsl_threshold_contrast')
    thresh.inputs.thresh = 0.95
    thresh.inputs.out_file = 'rando_pipe_thresh_tstat.nii.gz'
    wf.connect(select_tstat_files, 'out_file', thresh, 'in_file')
    wf.connect(thresh,'out_file',outputspec,'rando_pipe_thresh_tstat.nii.gz')

    thresh_bin = pe.Node(interface=fsl.UnaryMaths(),name='fsl_threshold_bin_contrast')
    thresh_bin.inputs.operation = 'bin'
    wf.connect(thresh, 'out_file', thresh_bin, 'in_file')
    wf.connect(thresh_bin,'out_file',outputspec,'thresh_bin_out')

    apply_mask = pe.Node(interface=fsl.ApplyMask(),name='fsl_applymask_contrast')
    wf.connect(select_tstat_files, 'out_file', apply_mask, 'in_file')
    wf.connect(thresh_bin, 'out_file', apply_mask, 'mask_file')


    cluster = pe.Node(interface=fsl.Cluster(),name='cluster_contrast')
    cluster.inputs.threshold = 0.0001
    cluster.inputs.out_index_file = "index_file"
    cluster.inputs.out_localmax_txt_file = "lmax_contrast.txt"
    cluster.inputs.out_size_file = "cluster_size_contrast"
    cluster.inputs.out_threshold_file = True
    cluster.inputs.out_max_file = True
    cluster.inputs.out_mean_file = True
    cluster.inputs.out_pval_file = True
    cluster.inputs.out_size_file = True
    

    wf.connect(apply_mask, 'out_file', cluster, 'in_file')
    
    wf.connect(cluster,'index_file',outputspec,'index_file')
    wf.connect(cluster,'threshold_file',outputspec,'threshold_file')
    wf.connect(cluster,'localmax_txt_file',outputspec,'localmax_txt_file')
    wf.connect(cluster,'localmax_vol_file',outputspec,'localmax_vol_file')
    wf.connect(cluster,'max_file',outputspec,'max_file')
    wf.connect(cluster,'mean_file',outputspec,'meal_file')
    wf.connect(cluster,'pval_file',outputspec,'pval_file')
    wf.connect(cluster,'size_file',outputspec,'size_file')
    
  

    return wf
