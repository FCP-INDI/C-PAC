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
    cmd = ['bash', '/code/CPAC/surface/fMRISurface/run.sh', '--post_freesurfer_folder', post_freesurfer_folder,\
        '--subject', subject, '--fmri', atlas_space_bold, '--scout', scout_bold,
        '--lowresmesh', low_res_mesh, '--grayordinatesres', gray_ordinates_res,
        '--fmrires', fmri_res, '--smoothingFWHM', smooth_fwhm]
    subprocess.check_output(cmd)

    out_file = os.path.join(post_freesurfer_folder, 'MNINonLinear/Results/task-rest01/task-rest01_Atlas.dtseries.nii')

    return out_file


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
                                output_names=['out_file'],
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

    node, out = strat_pool.get_data('freesurfer-subject-dir')
    wf.connect(node, out, surf, 'freesurfer_folder')

    node, out = strat_pool.get_data('desc-restore_T1w')
    wf.connect(node, out, surf, 't1w_restore_image')

    node, out = strat_pool.get_data('space-template_desc-head_T1w')
    wf.connect(node, out, surf, 'atlas_space_t1w_image')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, surf, 'atlas_transform')

    node, out = strat_pool.get_data('from-template_to-T1w_mode-image_xfm')
    wf.connect(node, out, surf, 'inverse_atlas_transform')

    node, out = strat_pool.get_data('space-template_desc-brain_bold')
    wf.connect(node, out, surf, 'atlas_space_bold')

    node, out = strat_pool.get_data('space-template_desc-scout_bold')
    wf.connect(node, out, surf, 'scout_bold')

    outputs = {
        'space-fsLR_den-32k_bold.dtseries': (surf, 'out_file')
    }

    return wf, outputs


def surface_preproc(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "surface_preproc",
     "config": ["surface_analysis", "post_freesurfer"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["freesurfer-subject-dir",
                "desc-restore_T1w",
                "space-template_desc-head_T1w",
                "from-T1w_to-template_mode-image_xfm",
                "from-template_to-T1w_mode-image_xfm",
                "space-template_desc-brain_bold",
                "space-template_desc-scout_bold"],
     "outputs": ["space-fsLR_den-32k_bold.dtseries"]}
    '''

    wf, outputs = surface_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)