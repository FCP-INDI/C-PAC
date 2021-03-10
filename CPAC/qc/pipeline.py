import os
import pkg_resources as p
from CPAC.pipeline import nipype_pipeline_engine as pe
from nipype.interfaces import afni
from nipype.interfaces import fsl
from CPAC.utils.interfaces.function import Function

from CPAC.qc.qc import (
    create_montage,
    create_montage_gm_wm_csf,
    qa_montages,
    create_qc_snr,
    create_qc_motion,
    create_qc_fd,
    create_qc_skullstrip,
    create_qc_carpet,
    afni_Edge3
)

from CPAC.qc.utils import (
    register_pallete,
    generate_qc_pages,
)

# register color palettes
palletes = ['red', 'green', 'blue', 'red_to_blue', 'cyan_to_yellow']
for pallete in palletes:
    register_pallete(
        p.resource_filename('CPAC', 'qc/colors/%s.txt' % pallete),
        pallete
    )


def qc_snr_plot(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_snr_plot",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [(["desc-brain_bold", "desc-motion_bold", "desc-preproc_bold"],
                 "space-bold_desc-brain_mask"),
                "from-bold_to-T1w_mode-image_desc-linear_xfm",
                "desc-brain_T1w",
                "space-T1w_desc-mean_bold"],
     "outputs": ["bold-snr-axial-qc",
                 "bold-snr-sagittal-qc",
                 "bold-snr-hist-qc",
                 "bold-snr-qc"]}
    '''

    # make SNR plot
    qc_workflow = create_qc_snr(f'qc_snr_{pipe_num}')

    node, out = strat_pool.get_data(["desc-brain_bold",
                                     "desc-motion_bold",
                                     "desc-preproc_bold"])
    wf.connect(node, out, qc_workflow, 'inputspec.functional_preprocessed')

    node, out = strat_pool.get_data("space-bold_desc-brain_mask")
    wf.connect(node, out, qc_workflow, 'inputspec.functional_brain_mask')

    node, out = \
        strat_pool.get_data('from-bold_to-T1w_mode-image_desc-linear_xfm')
    wf.connect(node, out,
               qc_workflow, 'inputspec.functional_to_anat_linear_xfm')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, qc_workflow, 'inputspec.anatomical_brain')

    node, out = strat_pool.get_data('space-T1w_desc-mean_bold')
    wf.connect(node, out, qc_workflow, 'inputspec.mean_functional_in_anat')

    outputs = {
        'bold-snr-axial-qc': (qc_workflow, 'outputspec.snr_axial_image'),
        'bold-snr-sagittal-qc':
            (qc_workflow, 'outputspec.snr_sagittal_image'),
        'bold-snr-hist-qc': (
            qc_workflow, 'outputspec.snr_histogram_image'),
        'bold-snr-qc': (qc_workflow, 'outputspec.snr_mean')
    }

    return (wf, outputs)


def qc_motion_plot(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_motion_plot",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["movement-parameters"],
     "outputs": ["movement-parameters-trans-qc",
                 "movement-parameters-rot-qc"]}
    '''

    # make motion parameters plot
    qc_workflow = create_qc_motion(f'qc_motion_{pipe_num}')

    node, out = strat_pool.get_data("movement-parameters")
    wf.connect(node, out, qc_workflow, 'inputspec.motion_parameters')

    outputs = {
        'movement-parameters-trans-qc': (
            qc_workflow, 'outputspec.motion_translation_plot'),
        'movement-parameters-rot-qc': (
            qc_workflow, 'outputspec.motion_rotation_plot')
    }

    return (wf, outputs)


def qc_fd_plot(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_fd_plot",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["framewise-displacement-jenkinson"],
     "outputs": ["framewise-displacement-jenkinson-plot-qc"]}
    '''

    qc_workflow = create_qc_fd(f'qc_fd_{pipe_num}')

    node, out = strat_pool.get_data('framewise-displacement-jenkinson')
    wf.connect(node, out, qc_workflow, 'inputspec.fd')

    outputs = {
        'framewise-displacement-jenkinson-plot-qc':
            (qc_workflow, 'outputspec.fd_histogram_plot')
    }

    return (wf, outputs)


def qc_brain_extraction(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_brain_extraction",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["desc-brain_T1w",
                "desc-reorient_T1w"],
     "outputs": ["desc-brain_T1w-axial-qc",
                 "desc-brain_T1w-sagittal-qc"]}
    '''

    # make QC montages for Skull Stripping Visualization
    qc_workflow = create_qc_skullstrip(
        f'qc_skullstrip_{pipe_num}'
    )

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, qc_workflow, 'inputspec.anatomical_brain')

    node, out = strat_pool.get_data('desc-reorient_T1w')
    wf.connect(node, out, qc_workflow, 'inputspec.anatomical_reorient')

    outputs = {
        'desc-brain_T1w-axial-qc': (qc_workflow, 'outputspec.axial_image'),
        'desc-brain_T1w-sagittal-qc':
            (qc_workflow, 'outputspec.sagittal_image')
    }

    return (wf, outputs)


def qc_T1w_standard(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_brain_extraction",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["space-template_desc-brain_T1w",
                "T1w_brain_template"],
     "outputs": ["space-template_desc-brain_T1w-axial-qc",
                 "space-template_desc-brain_T1w-sagittal-qc"]}
    '''

    # make QC montages for mni normalized anatomical image
    montage_mni_anat = create_montage(f'montage_mni_anat_{pipe_num}',
                                      'red', 'mni_anat',
                                      mapnode=False)

    node, out = strat_pool.get_data('space-template_desc-brain_T1w')
    wf.connect(node, out, montage_mni_anat, 'inputspec.underlay')

    anat_template_edge = pe.Node(Function(input_names=['in_file'],
                                          output_names=['out_file'],
                                          function=afni_Edge3,
                                          as_module=True),
                                 name=f'anat_template_edge_{pipe_num}')

    node, out = strat_pool.get_data('T1w_brain_template')
    wf.connect(node, out, anat_template_edge, 'in_file')

    wf.connect(anat_template_edge, 'out_file',
               montage_mni_anat, 'inputspec.overlay')

    outputs = {
        'space-template_desc-brain_T1w-axial-qc':
            (montage_mni_anat, 'outputspec.axial_png'),
        'space-template_desc-brain_T1w-sagittal-qc':
            (montage_mni_anat, 'outputspec.sagittal_png')
    }

    return (wf, outputs)


def qc_segmentation(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_segmentation",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [("desc-brain_T1w",
                 ["label-CSF_desc-preproc_mask",
                  "label-CSF_mask"],
                 ["label-WM_desc-preproc_mask",
                  "label-WM_mask"],
                 ["label-GM_desc-preproc_mask",
                  "label-GM_mask"])],
     "outputs": ["dseg-axial-qc",
                 "dseg-sagittal-qc"]}
    '''

    # make QC montages for CSF WM GM
    montage_csf_gm_wm = create_montage_gm_wm_csf(
        f'montage_csf_gm_wm_{pipe_num}', 'montage_csf_gm_wm')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, montage_csf_gm_wm, 'inputspec.underlay')

    node, out = strat_pool.get_data(['label-CSF_desc-preproc_mask',
                                     'label-CSF_mask'])
    wf.connect(node, out, montage_csf_gm_wm, 'inputspec.overlay_csf')

    node, out = strat_pool.get_data(['label-WM_desc-preproc_mask',
                                     'label-WM_mask'])
    wf.connect(node, out, montage_csf_gm_wm, 'inputspec.overlay_wm')

    node, out = strat_pool.get_data(['label-GM_desc-preproc_mask',
                                     'label-GM_mask'])
    wf.connect(node, out, montage_csf_gm_wm, 'inputspec.overlay_gm')

    outputs = {
        'dseg-axial-qc': (montage_csf_gm_wm, 'outputspec.axial_png'),
        'dseg-sagittal-qc': (montage_csf_gm_wm, 'outputspec.sagittal_png')
    }

    return (wf, outputs)
    
    
def qc_epi_segmentation(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_epi_segmentation",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [("desc-brain_bold",
                 ["space-bold_label-CSF_desc-preproc_mask",
                  "space-bold_label-CSF_mask"],
                 ["space-bold_label-WM_desc-preproc_mask",
                  "space-bold_label-WM_mask"],
                 ["space-bold_label-GM_desc-preproc_mask",
                  "space-bold_label-GM_mask"])],
     "outputs": ["epi-dseg-axial-qc",
                 "epi-dseg-sagittal-qc"]}
    '''

    # make QC montages for CSF WM GM
    montage_csf_gm_wm = create_montage_gm_wm_csf(
        f'montage_csf_gm_wm_{pipe_num}', 'montage_csf_gm_wm')

    node, out = strat_pool.get_data('desc-brain_bold')
    wf.connect(node, out, montage_csf_gm_wm, 'inputspec.underlay')

    node, out = strat_pool.get_data(['space-bold_label-CSF_desc-preproc_mask',
                                     'space-bold_label-CSF_mask'])
    wf.connect(node, out, montage_csf_gm_wm, 'inputspec.overlay_csf')

    node, out = strat_pool.get_data(['space-bold_label-WM_desc-preproc_mask',
                                     'space-bold_label-WM_mask'])
    wf.connect(node, out, montage_csf_gm_wm, 'inputspec.overlay_wm')

    node, out = strat_pool.get_data(['space-bold_label-GM_desc-preproc_mask',
                                     'space-bold_label-GM_mask'])
    wf.connect(node, out, montage_csf_gm_wm, 'inputspec.overlay_gm')

    outputs = {
        'epi-dseg-axial-qc': (montage_csf_gm_wm, 'outputspec.axial_png'),
        'epi-dseg-sagittal-qc': (montage_csf_gm_wm, 'outputspec.sagittal_png')
    }

    return (wf, outputs)


def qc_carpet_plot(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_carpet_plot",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [(["space-template_desc-cleaned_bold",
                  "space-template_desc-brain_bold",
                  "space-template_desc-preproc_bold",
                  "space-template_bold"],
                 "space-template_desc-mean_bold"),
                "GM_path",
                "WM_path",
                "CSF_path"],
     "outputs": ["space-template_desc-cleaned_bold-carpet-qc",
                 "space-template_desc-brain_bold-carpet-qc",
                 "space-template_desc-preproc_bold-carpet-qc",
                 "space-template_bold-carpet-qc"]}
    '''

    # make QC Carpet plot
    carpet_seg = create_qc_carpet(f'carpet_seg_{pipe_num}', 'carpet_seg')

    connection, resource = \
        strat_pool.get_data(["space-template_desc-cleaned_bold",
                             "space-template_desc-brain_bold",
                             "space-template_desc-preproc_bold",
                             "space-template_bold"],
                            report_fetched=True)
    node, out = connection
    wf.connect(node, out, carpet_seg, 'inputspec.functional_to_standard')

    node, out = strat_pool.get_data("space-template_desc-mean_bold")
    wf.connect(node, out, carpet_seg, 'inputspec.mean_functional_to_standard')

    node, out = strat_pool.get_data("GM_path")
    wf.connect(node, out, carpet_seg, 'inputspec.anatomical_gm_mask')

    node, out = strat_pool.get_data("WM_path")
    wf.connect(node, out, carpet_seg, 'inputspec.anatomical_wm_mask')

    node, out = strat_pool.get_data("CSF_path")
    wf.connect(node, out, carpet_seg, 'inputspec.anatomical_csf_mask')

    outputs = {
        f'{resource}-carpet-qc': (carpet_seg, 'outputspec.carpet_plot'),
    }

    return (wf, outputs)


def qc_coregistration(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_coregistration",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [("desc-brain_T1w",
                 "space-T1w_desc-mean_bold")],
     "outputs": ["space-T1w_desc-mean_bold-axial-qc",
                 "space-T1w_desc-mean_bold-sagittal-qc"]}
    '''

    # make QC montage for Mean Functional in T1 with T1 edge
    anat_edge = pe.Node(Function(input_names=['in_file'],
                                 output_names=['out_file'],
                                 function=afni_Edge3,
                                 as_module=True),
                        name=f'anat_edge_{pipe_num}')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, anat_edge, 'in_file')

    montage_anat = create_montage(f'montage_anat_{pipe_num}', 'red',
                                  't1_edge_on_mean_func_in_t1', 
                                  mapnode=False)

    wf.connect(anat_edge, 'out_file', montage_anat, 'inputspec.overlay')

    node, out = strat_pool.get_data('space-T1w_desc-mean_bold')
    wf.connect(node, out, montage_anat, 'inputspec.underlay')

    outputs = {
        'space-T1w_desc-mean_bold-axial-qc':
            (montage_anat, 'outputspec.axial_png'),
        'space-T1w_desc-mean_bold-sagittal-qc':
            (montage_anat, 'outputspec.sagittal_png')
    }

    return (wf, outputs)


def qc_bold_registration(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "qc_bold_registration",
     "config": ["pipeline_setup", "output_directory"],
     "switch": ["generate_quality_control_images"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["space-template_desc-mean_bold",
                "T1w_brain_template_funcreg"],
     "outputs": ["space-template_desc-mean_bold-axial-qc",
                 "space-template_desc-mean_bold-sagittal-qc"]}
    '''

    # make QC montage for Mean Functional in MNI with MNI edge
    montage_mfi = create_montage(f'montage_mfi_{pipe_num}', 'red',
                                 'MNI_edge_on_mean_func_mni',
                                 mapnode=False)

    node, out = strat_pool.get_data('space-template_desc-mean_bold')
    wf.connect(node, out,  montage_mfi, 'inputspec.underlay')

    func_template_edge = pe.Node(Function(input_names=['in_file'],
                                          output_names=['out_file'],
                                          function=afni_Edge3,
                                          as_module=True),
                                 name=f'func_template_edge_{pipe_num}')

    node, out = strat_pool.get_data("T1w_brain_template_funcreg")
    wf.connect(node, out, func_template_edge, 'in_file')

    wf.connect(func_template_edge, 'out_file',
               montage_mfi, 'inputspec.overlay')

    outputs = {
        'space-template_desc-mean_bold-axial-qc':
            (montage_mfi, 'outputspec.axial_png'),
        'space-template_desc-mean_bold-sagittal-qc':
            (montage_mfi, 'outputspec.sagittal_png')
    }

    return (wf, outputs)


def create_qc_workflow(cfg):
    qc_montage_id_a = {}
    qc_montage_id_s = {}
    qc_plot_id = {}
    qc_hist_id = {}

    qc_stack = []

    if cfg.functional_preproc['run'] and cfg.registration_workflows[
        'functional_registration']['coregistration']['run']:

        qc_stack += [qc_snr_plot, qc_coregistration]

        if not 0 in qc_montage_id_a:
            qc_montage_id_a[0] = 'snr_a'
            qc_montage_id_s[0] = 'snr_s'
            qc_hist_id[0] = 'snr_hist'

        if not 9 in qc_montage_id_a:
            qc_montage_id_a[9] = 'mean_func_with_t1_edge_a'
            qc_montage_id_s[9] = 'mean_func_with_t1_edge_s'

    if cfg.functional_preproc['run']:

        qc_stack += [qc_motion_plot, qc_fd_plot]

        if not 1 in qc_plot_id:
            qc_plot_id[1] = 'movement_trans_plot'
        if not 2 in qc_plot_id:
            qc_plot_id[2] = 'movement_rot_plot'
        if not 3 in qc_plot_id:
            qc_plot_id[3] = 'fd_plot'

    if cfg.anatomical_preproc['run']:

        qc_stack.append(qc_brain_extraction)

        if not 4 in qc_montage_id_a:
            qc_montage_id_a[4] = 'skullstrip_vis_a'
            qc_montage_id_s[4] = 'skullstrip_vis_s'

    if cfg.registration_workflows['anatomical_registration']['run']:

        qc_stack.append(qc_T1w_standard)

        if not 5 in qc_montage_id_a:
            qc_montage_id_a[5] = 'mni_normalized_anatomical_a'
            qc_montage_id_s[5] = 'mni_normalized_anatomical_s'

    if cfg.anatomical_preproc['run'] and cfg.segmentation['run']:

        if 'T1_Template' in cfg.segmentation['tissue_segmentation'][
            'Template_Based']['template_for_segmentation']:
            qc_stack.append(qc_segmentation)

        if 'EPI_Template' in cfg.segmentation['tissue_segmentation'][
            'Template_Based']['template_for_segmentation']:
            qc_stack.append(qc_epi_segmentation)

        if not 7 in qc_montage_id_a:
            qc_montage_id_a[7] = 'csf_gm_wm_a'
            qc_montage_id_s[7] = 'csf_gm_wm_s'

        if cfg.registration_workflows['functional_registration'][
            'func_registration_to_template']['run']:

            qc_stack.append(qc_carpet_plot)

            if not 8 in qc_plot_id:
                qc_plot_id[8] = 'carpet'

    if cfg.registration_workflows['functional_registration'][
        'func_registration_to_template']['run']:

        qc_stack.append(qc_bold_registration)

        if not 10 in qc_montage_id_a:
            qc_montage_id_a[10] = 'mean_func_with_mni_edge_a'
            qc_montage_id_s[10] = 'mean_func_with_mni_edge_s'

    return qc_stack, qc_montage_id_a, qc_montage_id_s, qc_hist_id, qc_plot_id
