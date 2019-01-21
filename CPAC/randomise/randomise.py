
import nipype.pipeline.engine as pe
from CPAC.pipeline.cpac_group_runner import load_config_yml


def select(input_list):
    import nibabel as nb
    for i in input_list:
        img = nb.load(i)
        hdr = img.header
        if hdr['cal_min'] == 0 and hdr['cal_max'] == 0:
            print ("Warning! {} is an empty image because of no positive "
                   "values in the unpermuted statistic image, and it could "
                   "not be processed with tfce.".format('i'))
        if not hdr['cal_max'] == 0 and hdr['cal_min'] == 0:
            selected_file = i

    return i


def prep_randomise_workflow(c, merged_file, mask_file, f_test, mat_file,
                            con_file, grp_file, output_dir, working_dir,
                            log_dir, model_name, fts_file=None):

    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl

    wf = pe.Workflow(name='randomise_workflow')
    wf.base_dir = c.working_dir

    outputspec = pe.Node(util.IdentityInterface(fields=['tstat_files',
                                                        't_corrected_p_files',
                                                        'index_file',
                                                        'threshold_file',
                                                        'localmax_txt_file',
                                                        'localmax_vol_file',
                                                        'max_file',
                                                        'mean_file',
                                                        'pval_file',
                                                        'size_file']),
                         name='outputspec')

    randomise = pe.Node(interface=fsl.Randomise(),
                        name='fsl-randomise_{0}'.format(model_name))
    randomise.inputs.base_name = model_name
    randomise.inputs.in_file = merged_file
    randomise.inputs.mask = mask_file
    randomise.inputs.num_perm = c.permutations
    randomise.inputs.demean = c.demean
    randomise.inputs.c_thresh = c.randomise_thresh
    randomise.inputs.tfce = c.tfce

    randomise.inputs.design_mat = mat_file
    randomise.inputs.tcon = con_file

    if fts_file:
        randomise.inputs.fcon = fts_file

    wf.connect(randomise,'tstat_files',outputspec,'tstat_files')
    wf.connect(randomise,'t_corrected_p_files',outputspec,'t_corrected_p_files')

    select_tcorrp_files = pe.Node(util.Function(input_names=['input_list'],
                                                output_names=['out_file'],
                                                function=select),
                                  name='select_t_corrp')

    wf.connect(randomise, 't_corrected_p_files', select_tcorrp_files, 'input_list')
    wf.connect(select_tcorrp_files,'out_file', outputspec,'out_tcorr_corrected')

    select_tstat_files = pe.Node(util.Function(input_names=['input_list'],
                                               output_names=['out_file'],
                                               function=select),
                                 name='select_t_stat')

    wf.connect(randomise, 'tstat_files',select_tstat_files, 'input_list')
    wf.connect(select_tstat_files,'out_file',outputspec,'out_tstat_corrected')

    thresh = pe.Node(interface=fsl.Threshold(),
                     name='fsl_threshold_contrast')
    thresh.inputs.thresh = 0.95
    thresh.inputs.out_file = 'randomise_pipe_thresh_tstat.nii.gz'
    wf.connect(select_tstat_files, 'out_file', thresh, 'in_file')
    wf.connect(thresh,'out_file',outputspec,'randomise_pipe_thresh_tstat.nii.gz')

    thresh_bin = pe.Node(interface=fsl.UnaryMaths(),
                         name='fsl_threshold_bin_contrast')
    thresh_bin.inputs.operation = 'bin'
    wf.connect(thresh, 'out_file', thresh_bin, 'in_file')
    wf.connect(thresh_bin,'out_file',outputspec,'thresh_bin_out')

    apply_mask = pe.Node(interface=fsl.ApplyMask(),
                         name='fsl_applymask_contrast')
    wf.connect(select_tstat_files, 'out_file', apply_mask, 'in_file')
    wf.connect(thresh_bin, 'out_file', apply_mask, 'mask_file')

    cluster = pe.Node(interface=fsl.Cluster(),
                      name='cluster_contrast')
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


def run(group_config_path):
    import re
    import commands
    commands.getoutput('source ~/.bashrc')
    import os
    import sys
    import pickle
    import yaml

    group_config_obj = load_config_yml(group_config_path)
    pipeline_output_folder = group_config_obj.pipeline_dir
    
    if not group_config_obj.participant_list == None:
        s_paths = group_config_obj.participant_list
    else:
        s_paths = [x for x in os.listdir(pipeline_output_folder) if os.path.isdir(x)]
    
    merged_file = randomise_merged_file(s_paths)

    out_file = randomise_merged_mask(s_paths)

    prep_randomise_workflow(group_config_obj, merged_file=merged_file,mask_file=out_file,working_dir=None,output_dir=None,crash_dir=None)
