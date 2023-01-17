from logging import raiseExceptions
import os
import nipype.interfaces.utility as util
from CPAC.utils.interfaces.function import Function
from CPAC.pipeline import nipype_pipeline_engine as pe



def run_surface(post_freesurfer_folder,
                freesurfer_folder,
                subject,
                t1w_restore_image,
                atlas_space_t1w_image,
                atlas_transform, 
                inverse_atlas_transform,
                atlas_space_bold,
                scout_bold,
                surf_atlas_dir,
                gray_ordinates_dir,
                gray_ordinates_res,
                high_res_mesh,
                low_res_mesh,
                subcortical_gray_labels,
                freesurfer_labels,
                fmri_res,
                smooth_fwhm):
    """
    Returns
    -------
    dtseries : str
        Path to the dtseries file.

    desikan_killiany : str
        Path to the Desikan-Killiany parcellation file.

    destrieux : str
        Path to the Destrieux parcellation file.
    """

    import os
    import subprocess

    freesurfer_folder = os.path.join(freesurfer_folder, 'recon_all')

    # DCAN-HCP PostFreeSurfer
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/PostFreeSurferPipeline.sh
    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/run.sh', '--post_freesurfer_folder', post_freesurfer_folder, \
        '--freesurfer_folder', freesurfer_folder, '--subject', subject, \
        '--t1w_restore', t1w_restore_image, '--atlas_t1w', atlas_space_t1w_image, \
        '--atlas_transform', atlas_transform, '--inverse_atlas_transform', inverse_atlas_transform, \
        '--surfatlasdir', surf_atlas_dir, '--grayordinatesdir', gray_ordinates_dir, '--grayordinatesres', gray_ordinates_res, \
        '--hiresmesh', high_res_mesh, '--lowresmesh', low_res_mesh, \
        '--subcortgraylabels', subcortical_gray_labels, '--freesurferlabels', freesurfer_labels]

    subprocess.check_output(cmd)

    # DCAN-HCP fMRISurface
    # https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh
    cmd = ['bash', '/code/CPAC/surface/fMRISurface/run.sh',
           '--post_freesurfer_folder', post_freesurfer_folder,
           '--subject', subject, '--fmri', atlas_space_bold, '--scout',
           scout_bold, '--lowresmesh', low_res_mesh, '--grayordinatesres',
           gray_ordinates_res, '--fmrires', fmri_res, '--smoothingFWHM',
           smooth_fwhm]
    subprocess.check_output(cmd)

    dtseries = os.path.join(post_freesurfer_folder,
                            'MNINonLinear/Results/task-rest01/'
                            'task-rest01_Atlas.dtseries.nii')
    aparc = {'desikan_killiany': {
        164: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                          f'{subject}.aparc.164k_fs_LR.dlabel.nii'),
        32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                         'fsaverage_LR32k',
                         f'{subject}.aparc.32k_fs_LR.dlabel.nii')},
             'destrieux': {
        164: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                          f'{subject}.aparc.a2009s.164k_fs_LR.dlabel.nii'),
        32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                         'fsaverage_LR32k',
                         f'{subject}.aparc.a2009s.32k_fs_LR.dlabel.nii')}}

    return (dtseries, aparc['desikan_killiany'][164], aparc['destrieux'][164],
            aparc['desikan_killiany'][32], aparc['destrieux'][32])


def surface_connector(wf, cfg, strat_pool, pipe_num, opt):
    
    surf = pe.Node(util.Function(input_names=['post_freesurfer_folder',
                                              'freesurfer_folder',
                                              'subject',
                                              't1w_restore_image',
                                              'atlas_space_t1w_image',
                                              'atlas_transform',
                                              'inverse_atlas_transform',
                                              'atlas_space_bold',
                                              'scout_bold',
                                              'surf_atlas_dir',
                                              'gray_ordinates_dir',
                                              'gray_ordinates_res',
                                              'high_res_mesh',
                                              'low_res_mesh',
                                              'subcortical_gray_labels',
                                              'freesurfer_labels',
                                              'fmri_res',
                                              'smooth_fwhm'],
                                 output_names=['dtseries',
                                               'desikan_killiany_164',
                                               'destrieux_164',
                                               'desikan_killiany_32',
                                               'destrieux_32',
                                              ],
                                 function=run_surface),
                   name=f'post_freesurfer_{pipe_num}')

    surf.inputs.subject = cfg['subject_id']

    surf.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
        'cpac_'+cfg['subject_id'],
        f'post_freesurfer_{pipe_num}')

    surf.inputs.surf_atlas_dir = cfg.surface_analysis['post_freesurfer']['surf_atlas_dir']
    surf.inputs.gray_ordinates_dir = cfg.surface_analysis['post_freesurfer']['gray_ordinates_dir']
    surf.inputs.subcortical_gray_labels = cfg.surface_analysis['post_freesurfer']['subcortical_gray_labels']
    surf.inputs.freesurfer_labels = cfg.surface_analysis['post_freesurfer']['freesurfer_labels']

    # convert integers to strings as subprocess requires string inputs
    surf.inputs.gray_ordinates_res = str(cfg.surface_analysis['post_freesurfer']['gray_ordinates_res'])
    surf.inputs.high_res_mesh = str(cfg.surface_analysis['post_freesurfer']['high_res_mesh'])
    surf.inputs.low_res_mesh = str(cfg.surface_analysis['post_freesurfer']['low_res_mesh'])
    surf.inputs.fmri_res = str(cfg.surface_analysis['post_freesurfer']['fmri_res'])
    surf.inputs.smooth_fwhm = str(cfg.surface_analysis['post_freesurfer']['smooth_fwhm'])

    restore = ["desc-restore_T1w", "desc-preproc_T1w", "desc-reorient_T1w", "T1w",
                  "space-longitudinal_desc-reorient_T1w"]
    space_temp = ["space-template_desc-head_T1w", "space-template_desc-brain_T1w", "space-template_desc-T1w_mask",]
    atlas_xfm = ["from-T1w_to-template_mode-image_xfm", "from-T1w_to-template_mode-image_desc-linear_xfm"]
    atlas_xfm_inv = ["from-template_to-T1w_mode-image_xfm", "from-template_to-T1w_mode-image_desc-linear_xfm"]
    atlas_space_bold = ["space-template_desc-brain_bold", "space-template_desc-preproc_bold"]
    scout_bold = ["space-template_desc-scout_bold", "space-template_desc-cleaned_bold", "space-template_desc-brain_bold",
                  "space-template_desc-preproc_bold", "space-template_desc-motion_bold", "space-template_bold"]


    node, out = strat_pool.get_data('freesurfer-subject-dir')
    wf.connect(node, out, surf, 'freesurfer_folder')

    node, out = strat_pool.get_data(restore) 
    wf.connect(node, out, surf, 't1w_restore_image')

    node, out = strat_pool.get_data(space_temp) 
    wf.connect(node, out, surf, 'atlas_space_t1w_image')

    node, out = strat_pool.get_data(atlas_xfm) 
    wf.connect(node, out, surf, 'atlas_transform')

    node, out = strat_pool.get_data(atlas_xfm_inv) 
    wf.connect(node, out, surf, 'inverse_atlas_transform')

    node, out = strat_pool.get_data(atlas_space_bold) 
    wf.connect(node, out, surf, 'atlas_space_bold')

    node, out = strat_pool.get_data(scout_bold)
    wf.connect(node, out, surf, 'scout_bold')
    
    falff = pe.Node(util.Function(input_names=['dtseries'], 
                                output_names=['falff'],
                                function=run_surf_falff),
                   name=f'surf_falff_{pipe_num}')
    #falff.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
    #    'cpac_'+cfg['subject_id'],
    #    f'post_freesurfer_{pipe_num}')
   
    wf.connect(surf, 'dtseries', falff, 'dtseries')
    
    alff = pe.Node(util.Function(input_names=['dtseries'], 
                                output_names=['alff'],
                                function=run_surf_alff),
                  name=f'surf_alff_{pipe_num}')

  
    wf.connect(surf, 'dtseries', alff, 'dtseries')


    outputs = {
        'atlas-DesikanKilliany_space-fsLR_den-32k_dlabel': (surf,
                                                            'desikan_'
                                                            'killiany_32'),
        'atlas-Destrieux_space-fsLR_den-32k_dlabel': (surf, 'destrieux_32'),
        'atlas-DesikanKilliany_space-fsLR_den-164k_dlabel': (surf,
                                                             'desikan_'
                                                             'killiany_164'),
        'atlas-Destrieux_space-fsLR_den-164k_dlabel': (surf, 'destrieux_164'),
        'space-fsLR_den-32k_bold-dtseries': (surf, 'dtseries'),
        'falff-surf_dscalar': (falff, 'falff'),
        'alff-surf_dscalar': (alff, 'alff')


    }

    
    
  
    # L_cortex_file = pe.Node(util.Function(input_names=['dtseries', 'structure', 'post_freesurfer_folder', 'cortex_filename'], 
    #                             output_names=['L_cortex_file'],
    #                             function=run_get_cortex),
    #               name=f'L_surf_cortex_{pipe_num}')

    # L_cortex_file.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
    #     'cpac_'+cfg['subject_id'],
    #     f'post_freesurfer_{pipe_num}')
    # L_cortex_file.inputs.structure = "LEFT"
    # L_cortex_file.inputs.cortex_filename = "L_cortex.func.gii"
    # wf.connect(surf, 'dtseries', L_cortex_file, 'dtseries')

    # R_cortex_file = pe.Node(util.Function(input_names=['dtseries', 'structure', 'post_freesurfer_folder', 'cortex_filename'], 
    #                             output_names=['R_cortex_file'],
    #                             function=run_get_cortex),
    #               name=f'R_surf_cortex_{pipe_num}')

    # R_cortex_file.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
    #     'cpac_'+cfg['subject_id'],
    #     f'post_freesurfer_{pipe_num}')
    # R_cortex_file.inputs.structure = "RIGHT"
    # R_cortex_file.inputs.cortex_filename = "R_cortex.func.gii"
    # wf.connect(surf, 'dtseries', R_cortex_file, 'dtseries')


    # mean_timeseries = pe.Node(util.Function(input_names=['post_freesurfer_folder', 'dtseries'], 
    #                             output_names=['mean_timeseries'],
    #                             function=run_mean_timeseries),
    #               name=f'mean_timeseries_{pipe_num}')

    # mean_timeseries.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
    #     'cpac_'+cfg['subject_id'],
    #     f'post_freesurfer_{pipe_num}')
    # wf.connect(surf, 'dtseries', mean_timeseries, 'dtseries')

    
    # L_reho = pe.Node(util.Function(input_names=['dtseries', 'mask' , 'cortex_file', 'surface_file', 'mean_timeseries', 
    #                                         'post_freesurfer_folder', 'reho_filename'], 
    #                             output_names=['L_reho'],
    #                             function=run_surf_reho),
    #                 name=f'surf_reho{pipe_num}')

    # wf.connect(get_L_cortex_file, 'L_cortex_file', L_reho, 'cortex_file')
    # wf.connect(surf, 'L_surface_file', L_reho, 'surface_file')
    # wf.connect(surf, 'L_mask', L_reho, 'mask')
    # wf.connect(mean_timeseries, 'mean_timeseries', L_reho, 'mean_timeseries') 
    # L_reho.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
    #     'cpac_'+cfg['subject_id'],
    #     f'post_freesurfer_{pipe_num}')
    # L_reho.inputs.reho_filename = L_surf_reho.dscalar.nii

    # R_reho = pe.Node(util.Function(input_names=['dtseries', 'mask' , 'cortex_file', 'surface_file', 'mean_timeseries', 
    #                                         'post_freesurfer_folder', 'reho_filename'], 
    #                             output_names=['R_reho'],
    #                             function=run_surf_reho),
    #                 name=f'surf_reho{pipe_num}')

    # wf.connect(get_R_cortex_file, 'R_cortex_file', R_reho, 'cortex_file')
    # wf.connect(surf, 'R_surface_file', R_reho, 'surface_file')
    # wf.connect(surf, 'R_mask', R_reho, 'mask')
    # wf.connect(mean_timeseries, 'mean_timeseries', R_reho, 'mean_timeseries')  
    # R_reho.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
    #     'cpac_'+cfg['subject_id'],
    #     f'post_freesurfer_{pipe_num}')
    # R_reho.inputs.reho_filename = R_surf_reho.dscalar.nii


    
    # connectivity_parcellation = pe.Node(util.Function(input_names=['dtseries', 'surf_atlaslabel',
    #                                             'post_freesurfer_folder'], 
    #                             output_names=['parcellation_file'],
    #                             function=run_ciftiparcellate),
    #                         name=f'connectivity_parcellation{pipe_num}')
                
    
    # wf.connect(surf, 'dtseries', connectivity, 'dtseries') 
    # connectivity_parcellation.inputs.surf_atlaslabel = ## path to the label file

    # correlation_matrix = pe.Node(util.Function(input_names=['ptseries','post_freesurfer_folder'], 
    #                             output_names=['correlation_matrix'],
    #                             function=run_cifticorrelation),
    #                         name=f'correlation_matrix{pipe_num}')
                
    
    # wf.connect(connectivity_parcellation, 'parcellation_file', correlation_matrix 'ptseries') 

    return wf, outputs

def surface_postproc(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "surface_preproc",
     "config": ["surface_analysis", "post_freesurfer"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["freesurfer-subject-dir",
                ["desc-restore_T1w", "desc-preproc_T1w", "desc-reorient_T1w", "T1w", 
                "space-longitudinal_desc-reorient_T1w"],
                ["space-template_desc-head_T1w", "space-template_desc-brain_T1w", "space-template_desc-T1w_mask"],
                ["from-T1w_to-template_mode-image_xfm", "from-T1w_to-template_mode-image_desc-linear_xfm"],
                ["from-template_to-T1w_mode-image_xfm", "from-template_to-T1w_mode-image_desc-linear_xfm"],
                ["space-template_desc-brain_bold", "space-template_desc-preproc_bold"],
                ["space-template_desc-scout_bold", "space-template_desc-cleaned_bold", "space-template_desc-brain_bold", 
                "space-template_desc-preproc_bold", "space-template_desc-motion_bold", "space-template_bold"]],
     "outputs": ["atlas-DesikanKilliany_space-fsLR_den-32k_dlabel",
                 "atlas-Destrieux_space-fsLR_den-32k_dlabel",
                 "atlas-DesikanKilliany_space-fsLR_den-164k_dlabel",
                 "atlas-Destrieux_space-fsLR_den-164k_dlabel",
                 "space-fsLR_den-32k_bold-dtseries",
                 "falff-surf_dscalar",
                 "alff-surf_dscalar"]}
    '''
    wf, outputs = surface_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


def run_surf_falff(dtseries):
    import os
    import subprocess
    falff = os.path.join(os.getcwd(), 'falff_surf.dscalar.nii.gz')
    cmd = ['ciftify_falff', dtseries, falff, '--min-low-freq', '0.01', '--max-low-freq' , '0.1']
    subprocess.check_output(cmd)
    return falff


def run_surf_alff(dtseries):
    import os
    import subprocess
    alff = os.path.join(os.getcwd(), 'alff_surf.dscalar.nii.gz')
    cmd = ['ciftify_falff', dtseries, alff, '--min-low-freq', '0.01', '--max-low-freq' , '0.1' , '--calc-alff']
    subprocess.check_output(cmd)
    return alff
    
#cmd = ['ciftify_falff', dtseries , 'alff_surf.dscalar.nii', '--min-low-freq', '0.01', '--max-low-freq' , '0.1' , '--calc-alff']
# def run_get_cortex(dtseries, structure, post_freesurfer_folder, cortex_filename):
#     import os
#     import subprocess
#     cmd = ['wb_command', '-cifti-separate', dtseries , 'COLUMN', '-label', structure, cortex_filename]
#     subprocess.check_output(cmd)
#     cortex_file = os.path.join(post_freesurfer_folder, cortex_filename)
#     return cortex_file

# def run_mean_timeseries(dtseries,post_freesurfer_folder):

#     import os
#     import subprocess
#     cmd = ['wb_command', '-cifti-reduce', dtseries, 'MEAN', 'mean.dscalar.nii']
#     subprocess.check_output(cmd)
#     mean_timeseries = os.path.join(post_freesurfer_folder, mean.dscalar.nii)
#     return mean_timeseries
    
# def run_surf_reho(dtseries, mask, cortex_file, surface_file,mean_timeseries,post_freesurfer_folder, reho_filename):

#     import os
#     import subprocess
#     cmd = ['python', '/code/CPAC/surface/PostFreeSurfer/surf_reho.py', dtseries, mask, cortex_file, surface_file, mean_timeseries, reho_filename]
#     subprocess.check_output(cmd)
#     surf_reho = os.path.join(post_freesurfer_folder, reho_filename)
#     return surf_reho
    
# def run_ciftiparcellate(dtseries, surf_atlaslabel):
#     import os  
#     import subprocess
#     cmd = ['wb_command', '-cifti-parcellate', dtseries , surf_atlaslabel, 'COLUMN', 'parcellation.ptseries.nii']
#     subprocess.check_output(cmd)
#     parcellation_file = os.path.join(post_freesurfer_folder, 'parcellation.ptseries.nii')
#     return parcellation_file

# def run_cifticorrelation(ptseries):
#     import os  
#     import subprocess
#     cmd = ['wb_command', '-cifti-correlation ', ptseries , 'cifti_corr.pconn.nii']
#     subprocess.check_output(cmd)
#     correlation_matrix = os.path.join(post_freesurfer_folder, 'cifti_corr.pconn.nii')
#     return correlation_matrix
