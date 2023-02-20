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

    #freesurfer_folder = os.path.join(freesurfer_folder, 'recon_all')
    print('run_surface 1')

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

    print('run_surface 12')
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
                                               'destrieux_32'],
                                 function=run_surface),
                   name=f'post_freesurfer_{pipe_num}')

    
    surf.inputs.subject = cfg['subject_id']
    print('surf_postproc1')
    surf.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
        'cpac_'+cfg['subject_id'],
        f'post_freesurfer_{pipe_num}')
    print('surf_postproc2')
    surf.inputs.surf_atlas_dir = cfg.surface_analysis['post_freesurfer']['surf_atlas_dir']
    surf.inputs.gray_ordinates_dir = cfg.surface_analysis['post_freesurfer']['gray_ordinates_dir']
    surf.inputs.subcortical_gray_labels = cfg.surface_analysis['post_freesurfer']['subcortical_gray_labels']
    surf.inputs.freesurfer_labels = cfg.surface_analysis['post_freesurfer']['freesurfer_labels']
    print('surf_postproc3')
    # convert integers to strings as subprocess requires string inputs
    surf.inputs.gray_ordinates_res = str(cfg.surface_analysis['post_freesurfer']['gray_ordinates_res'])
    surf.inputs.high_res_mesh = str(cfg.surface_analysis['post_freesurfer']['high_res_mesh'])
    surf.inputs.low_res_mesh = str(cfg.surface_analysis['post_freesurfer']['low_res_mesh'])
    surf.inputs.fmri_res = str(cfg.surface_analysis['post_freesurfer']['fmri_res'])
    surf.inputs.smooth_fwhm = str(cfg.surface_analysis['post_freesurfer']['smooth_fwhm'])

    restore = ["pipeline-fs_desc-restore_T1w", "desc-preproc_T1w", "desc-reorient_T1w", "T1w",
                  "space-longitudinal_desc-reorient_T1w"]
    space_temp = ["space-template_desc-head_T1w", "space-template_desc-brain_T1w", "space-template_desc-T1w_mask"]
    atlas_xfm = ["from-T1w_to-template_mode-image_xfm", "from-T1w_to-template_mode-image_desc-linear_xfm"]
    atlas_xfm_inv = ["from-template_to-T1w_mode-image_xfm", "from-template_to-T1w_mode-image_desc-linear_xfm"]
    atlas_space_bold = ["space-template_desc-brain_bold", "space-template_desc-preproc_bold"]
    scout_bold = ["space-template_desc-scout_bold", "space-template_desc-cleaned_bold", "space-template_desc-brain_bold",
                  "space-template_desc-preproc_bold", "space-template_desc-motion_bold", "space-template_bold"]

    print('surf_postproc4')
    node, out = strat_pool.get_data('freesurfer-subject-dir')
    wf.connect(node, out, surf, 'freesurfer_folder')
    
    print('surf_postproc5')
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
    print('surf_postproc6')
    outputs = {
        'atlas-DesikanKilliany_space-fsLR_den-32k_dlabel': (surf,
                                                            'desikan_'
                                                            'killiany_32'),
        'atlas-Destrieux_space-fsLR_den-32k_dlabel': (surf, 'destrieux_32'),
        'atlas-DesikanKilliany_space-fsLR_den-164k_dlabel': (surf,
                                                             'desikan_'
                                                             'killiany_164'),
        'atlas-Destrieux_space-fsLR_den-164k_dlabel': (surf, 'destrieux_164'),
        'space-fsLR_den-32k_bold-dtseries': (surf, 'dtseries')
    }

    return wf, outputs

def surface_postproc(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "surface_postproc",
     "config": ["surface_analysis", "post_freesurfer"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["freesurfer-subject-dir",
                (["pipeline-fs_desc-restore_T1w", "desc-preproc_T1w", "desc-reorient_T1w", "T1w", 
                "space-longitudinal_desc-reorient_T1w"],
                ["space-template_desc-head_T1w", "space-template_desc-brain_T1w", "space-template_desc-T1w_mask"],
                ["from-T1w_to-template_mode-image_xfm", "from-T1w_to-template_mode-image_desc-linear_xfm"],
                ["from-template_to-T1w_mode-image_xfm", "from-template_to-T1w_mode-image_desc-linear_xfm"],
                ["space-template_desc-brain_bold", "space-template_desc-preproc_bold"],
                ["space-template_desc-scout_bold", "space-template_desc-cleaned_bold", "space-template_desc-brain_bold", 
                "space-template_desc-preproc_bold", "space-template_desc-motion_bold", "space-template_bold"])],
     "outputs": ["atlas-DesikanKilliany_space-fsLR_den-32k_dlabel",
                 "atlas-Destrieux_space-fsLR_den-32k_dlabel",
                 "atlas-DesikanKilliany_space-fsLR_den-164k_dlabel",
                 "atlas-Destrieux_space-fsLR_den-164k_dlabel",
                 "space-fsLR_den-32k_bold-dtseries"]}
    '''
    wf, outputs = surface_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)

