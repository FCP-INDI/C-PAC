import os
import pkg_resources as p
import nipype.pipeline.engine as pe
from nipype.interfaces import afni

from CPAC.utils.function import Function

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


def create_qc_workflow(workflow, c, strategies, qc_outputs):

    qc_montage_id_a = {}
    qc_montage_id_s = {}
    qc_plot_id = {}
    qc_hist_id = {}

    for num_strat, strat in enumerate(strategies):

        nodes = strat.get_nodes_names()

        if 'functional_preprocessed' in strat:

            preproc, out_file = strat['functional_preprocessed']
            brain_mask, mask_file = strat['functional_brain_mask']
            func_to_anat_xfm, xfm_file = strat['functional_to_anat_linear_xfm']
            anat_ref, ref_file = strat['anatomical_brain']
            mfa, mfa_file = strat['mean_functional_in_anat']

            # make SNR plot
            qc_workflow = create_qc_snr('qc_snr_{0}'.format(num_strat))
            workflow.connect(preproc, out_file, qc_workflow, 'inputspec.functional_preprocessed')
            workflow.connect(brain_mask, mask_file, qc_workflow, 'inputspec.functional_brain_mask')
            workflow.connect(func_to_anat_xfm, xfm_file, qc_workflow, 'inputspec.functional_to_anat_linear_xfm')
            workflow.connect(anat_ref, ref_file, qc_workflow, 'inputspec.anatomical_brain')
            workflow.connect(mfa, mfa_file, qc_workflow, 'inputspec.mean_functional_in_anat')

            strat.update_resource_pool({
                'qc___snr_a': (qc_workflow, 'outputspec.snr_axial_image'),
                'qc___snr_s': (qc_workflow, 'outputspec.snr_sagittal_image'),
                'qc___snr_hist': (qc_workflow, 'outputspec.snr_histogram_image'),
                'qc___snr_val': (qc_workflow, 'outputspec.snr_mean')
            })

            if not 3 in qc_montage_id_a:
                qc_montage_id_a[3] = 'snr_a'
                qc_montage_id_s[3] = 'snr_s'
                qc_hist_id[3] = 'snr_hist'

            # make motion parameters plot
            mov_param, out_file = strat['movement_parameters']

            qc_workflow = create_qc_motion('qc_motion_{0}'.format(num_strat))
            workflow.connect(mov_param, out_file, qc_workflow,
                                'inputspec.motion_parameters')

            strat.update_resource_pool({
                'qc___movement_trans_plot': (qc_workflow, 'outputspec.motion_translation_plot'),
                'qc___movement_rot_plot': (qc_workflow, 'outputspec.motion_rotation_plot')
            })

            if not 6 in qc_plot_id:
                qc_plot_id[6] = 'movement_trans_plot'

            if not 7 in qc_plot_id:
                qc_plot_id[7] = 'movement_rot_plot'

            # make FD plot and volumes removed
            if 'gen_motion_stats' in nodes and 1 in c.runNuisance:
                if c.fdCalc == 'Power':
                    fd, out_file = strat['frame_wise_displacement_power']
                else:
                    fd, out_file = strat['frame_wise_displacement_jenkinson']

                qc_workflow = create_qc_fd('qc_fd_{0}'.format(num_strat))

                workflow.connect(fd, out_file, qc_workflow, 'inputspec.fd')

                if "De-Spiking" in c.runMotionSpike:
                    excluded, out_file_ex = strat['despiking_frames_excluded']
                    workflow.connect(excluded, out_file_ex,
                                        qc_workflow, 'inputspec.excluded_volumes')

                elif "Scrubbing" in c.runMotionSpike:
                    excluded, out_file_ex = strat['scrubbing_frames_excluded']
                    workflow.connect(excluded, out_file_ex,
                                        qc_workflow, 'inputspec.excluded_volumes')

                strat.update_resource_pool({
                    'qc___fd_plot': (qc_workflow, 'outputspec.fd_histogram_plot')
                })

                if not 8 in qc_plot_id:
                    qc_plot_id[8] = 'fd_plot'

        # make QC montages for Skull Stripping Visualization
        anat_underlay, out_file = strat['anatomical_brain']
        skull, out_file_s = strat['anatomical_reorient']

        qc_workflow = create_qc_skullstrip(
            'qc_skullstrip_{0}'.format(num_strat)
        )
        workflow.connect(anat_underlay, out_file,
                            qc_workflow, 'inputspec.anatomical_brain')
        workflow.connect(skull, out_file_s,
                            qc_workflow, 'inputspec.anatomical_reorient')

        strat.update_resource_pool({
            'qc___skullstrip_vis_a': (qc_workflow, 'outputspec.axial_image'),
            'qc___skullstrip_vis_s': (qc_workflow, 'outputspec.sagittal_image')
        })

        if not 1 in qc_montage_id_a:
            qc_montage_id_a[1] = 'skullstrip_vis_a'
            qc_montage_id_s[1] = 'skullstrip_vis_s'

        if 'functional_preprocessed' in strat:

            # make QC montages for mni normalized anatomical image
            mni_anat_underlay, out_file = strat['mean_functional_in_anat']

            montage_mni_anat = create_montage('montage_mni_anat_{0}'.format(num_strat),
                                                'red',
                                                'mni_anat')

            workflow.connect(mni_anat_underlay, out_file,
                                montage_mni_anat, 'inputspec.underlay')
            montage_mni_anat.inputs.inputspec.overlay = os.path.abspath(
                p.resource_filename(
                    'CPAC', 'resources/templates/MNI152_Edge_AllTissues.nii.gz'
                )
            )

            strat.update_resource_pool({'qc___mni_normalized_anatomical_a': (montage_mni_anat, 'outputspec.axial_png'),
                                        'qc___mni_normalized_anatomical_s': (montage_mni_anat, 'outputspec.sagittal_png')})

            if not 6 in qc_montage_id_a:
                qc_montage_id_a[6] = 'mni_normalized_anatomical_a'
                qc_montage_id_s[6] = 'mni_normalized_anatomical_s'

        # make QC montages for CSF WM GM
        if 'seg_preproc' in nodes:

            anat_underlay, out_file_anat = strat['anatomical_brain']
            csf_overlay, out_file_csf = strat['anatomical_csf_mask']
            wm_overlay, out_file_wm = strat['anatomical_wm_mask']
            gm_overlay, out_file_gm = strat['anatomical_gm_mask']

            montage_csf_gm_wm = create_montage_gm_wm_csf('montage_csf_gm_wm_%d' % num_strat,
                                                            'montage_csf_gm_wm')

            workflow.connect(anat_underlay, out_file_anat,
                                montage_csf_gm_wm, 'inputspec.underlay')
            workflow.connect(csf_overlay, out_file_csf,
                                montage_csf_gm_wm, 'inputspec.overlay_csf')
            workflow.connect(wm_overlay, out_file_wm,
                                montage_csf_gm_wm, 'inputspec.overlay_wm')
            workflow.connect(gm_overlay, out_file_gm,
                                montage_csf_gm_wm, 'inputspec.overlay_gm')

            strat.update_resource_pool({
                'qc___csf_gm_wm_a': (montage_csf_gm_wm, 'outputspec.axial_png'),
                'qc___csf_gm_wm_s': (montage_csf_gm_wm, 'outputspec.sagittal_png')
            })

            if not 2 in qc_montage_id_a:
                qc_montage_id_a[2] = 'csf_gm_wm_a'
                qc_montage_id_s[2] = 'csf_gm_wm_s'

            if 'functional_preprocessed' in strat:
                preproc, out_file_preproc = strat['functional_to_standard']
                mean_preproc, out_file_mean_preproc = strat['mean_functional_to_standard']

                # make QC Carpet plot
                carpet_seg = create_qc_carpet('carpet_seg_%d' % num_strat,
                                              'carpet_seg')

                workflow.connect(
                    preproc, out_file_preproc,
                    carpet_seg, 'inputspec.functional_to_standard'
                )
                workflow.connect(
                    mean_preproc, out_file_mean_preproc,
                    carpet_seg, 'inputspec.mean_functional_to_standard'
                )


                workflow.connect(
                    c.PRIORS_GRAY, 'local_path',
                    carpet_seg, 'inputspec.anatomical_gm_mask'
                )
                workflow.connect(
                    c.PRIORS_WHITE, 'local_path',
                    carpet_seg, 'inputspec.anatomical_wm_mask'
                )
                workflow.connect(
                    c.PRIORS_CSF, 'local_path',
                    carpet_seg, 'inputspec.anatomical_csf_mask'
                )

                strat.update_resource_pool({
                    'qc___carpet': (carpet_seg, 'outputspec.carpet_plot'),
                })

                if not 9 in qc_plot_id:
                    qc_plot_id[9] = 'carpet'

        if 'functional_preprocessed' in strat:
                
            # make QC montage for Mean Functional in T1 with T1 edge
            anat, out_file = strat['anatomical_brain']
            m_f_a, out_file_mfa = strat['mean_functional_in_anat']

            anat_edge = pe.Node(Function(input_names=['in_file'],
                                        output_names=['out_file'],
                                        function=afni_Edge3,
                                        as_module=True),
                                name='anat_edge_%d' % num_strat)

            montage_anat = create_montage(
                'montage_anat_%d' % num_strat, 'red', 't1_edge_on_mean_func_in_t1')

            workflow.connect(anat, out_file, anat_edge, 'in_file')
            workflow.connect(anat_edge, 'out_file',
                                montage_anat, 'inputspec.overlay')
            workflow.connect(m_f_a, out_file_mfa,
                                montage_anat, 'inputspec.underlay')

            strat.update_resource_pool({
                'qc___mean_func_with_t1_edge_a': (montage_anat, 'outputspec.axial_png'),
                'qc___mean_func_with_t1_edge_s': (montage_anat, 'outputspec.sagittal_png')
            })

            if not 4 in qc_montage_id_a:
                qc_montage_id_a[4] = 'mean_func_with_t1_edge_a'
                qc_montage_id_s[4] = 'mean_func_with_t1_edge_s'

            # make QC montage for Mean Functional in MNI with MNI edge
            m_f_i, out_file = strat['mean_functional_to_standard']

            montage_mfi = create_montage(
                'montage_mfi_%d' % num_strat, 'red', 'MNI_edge_on_mean_func_mni')
            workflow.connect(m_f_i, out_file, montage_mfi,
                                'inputspec.underlay')
            montage_mfi.inputs.inputspec.overlay =  os.path.abspath(
                p.resource_filename(
                    'CPAC', 'resources/templates/MNI152_Edge_AllTissues.nii.gz'
                )
            )

            strat.update_resource_pool({'qc___mean_func_with_mni_edge_a': (montage_mfi, 'outputspec.axial_png'),
                                        'qc___mean_func_with_mni_edge_s': (montage_mfi, 'outputspec.sagittal_png')})

            if not 5 in qc_montage_id_a:
                qc_montage_id_a[5] = 'mean_func_with_mni_edge_a'
                qc_montage_id_s[5] = 'mean_func_with_mni_edge_s'

        # Link all the derivatives to the QC pages
        idx = 7
        rp = strat.get_resource_pool()
        for key in sorted(rp.keys()):
            # qc_outputs is from the outputs CSV
            if key in qc_outputs:
                qa_montages(workflow, c, strat, num_strat,
                            qc_montage_id_a, qc_montage_id_s, qc_hist_id,
                            key, idx)
                idx += 1

    return qc_montage_id_a, qc_montage_id_s, qc_hist_id, qc_plot_id
