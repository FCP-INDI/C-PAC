import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe
import re
import os
import sys
import glob
from CPAC.utils import Configuration
from CPAC.pipeline.cpac_group_runner import load_config_yml


def load_subject_file(group_config_path):
    group_config_obj = load_config_yml(group_config_path)
    pipeline_output_folder = group_config_obj.pipeline_dir
    
    if not group_config_obj.participant_list == None:
        s_paths = group_config_obj.participant_list
    else:
        s_paths = x for x in os.listdir(pipeline_output_folder) if os.path.isdir(x) 
    return s_paths

def randomise_merged_file(s_paths):
    
    merge = pe.Node(interface=fsl.Merge(), name='fsl_merge')
    merge.inputs.dimension = 't'
    merge.inputs.merged_file = "randomise_merged.nii.gz"
    merge.inputs.in_files = s_paths   

    return merged_file

def randomise_merged_mask(s_paths):

    mask = pe.Node(interface=fsl.maths.MathsCommand(), name='fsl_maths')
    mask.inputs.args = '-abs -Tmin -bin'
    mask.inputs.out_file = "randomise_mask.nii.gz"
    mask.inputs.in_file = s_paths

    return out_file

def prep_randomise_workflow(c, subject_infos):
    print 'Preparing Randomise workflow'
    #p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    #print 'Subjects', s_ids

    wf = pe.Workflow(name='randomise_workflow')
    wf.base_dir = c.workingDirectory

    from CPAC.randomise import create_randomise
    import numpy as np
    
    rw = create_randomise()

    inputspec = pe.Node(util.IdentityInterface(fields=['permutations','demean','c_thresh','contrast_file','design_matrix_file']),name='inputspec')
    
    outputspec = pe.Node(util.IdentityInterface(fields=['tstat_files' ,'t_corrected_p_files','index_file','threshold_file','localmax_txt_file','localmax_vol_file','max_file','mean_file','pval_file','size_file']), name='outputspec')

    randomise = pe.Node(interface=fsl.Randomise(), name='randomise')
    randomise.inputs.base_name = "randomise"
    randomise.inputs.subjects = merged_file 
    randomise.inputs.mask = out_file
    randomise.inputs.demean = True
    randomise.inputs.tfce = True
    
    wf.connect(inputspec,'design_matrix_file',randomsie,'design_mat')
    wf.connect(inputspec,'constrast_file',randomise,'tcon')
    wf.connect(inputspec,'permutations',randomise,'num_perm')
    wf.connect(randomise,'tstat_files',outputspec,'tstat_files')
    wf.connect(randomise,'t_corrected_p_files',outputspec,'t_corrected_p_files')

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
