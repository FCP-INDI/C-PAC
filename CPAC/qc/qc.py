import matplotlib
matplotlib.use('Agg', warn=False, force=True)

from nipype.interfaces import afni 
from CPAC.utils.interfaces.function import Function
from CPAC.qc.utils import (
    resample_1mm, montage_axial, montage_sagittal,
    montage_gm_wm_csf_axial, montage_gm_wm_csf_sagittal,
    cal_snr_val, gen_histogram, drop_percent, gen_motion_plt,
    gen_plot_png,
    gen_carpet_plt
)

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from nipype.interfaces import afni
import nipype.interfaces.fsl as fsl


def create_montage(wf_name, cbar_name, png_name):

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['underlay',
                                                       'overlay']),
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['axial_png',
                                                        'sagittal_png',
                                                        'resampled_underlay',
                                                        'resampled_overlay']),
                         name='outputspec')

    # node for resampling create_montage images to 1mm for QC pages
    resample_u = pe.Node(Function(input_names=['file_'],
                                  output_names=['new_fname'],
                                  function=resample_1mm,
                                  as_module=True),
                         name='resample_u')
    
    wf.connect(inputnode, 'underlay', resample_u, 'file_')
    wf.connect(resample_u, 'new_fname', outputnode,'resampled_underlay')

    # same for overlays (resampling to 1mm)
    resample_o = pe.Node(Function(input_names=['file_'],
                                  output_names=['new_fname'],
                                  function=resample_1mm,
                                  as_module=True),
                         name='resample_o')

    wf.connect(inputnode, 'overlay', resample_o, 'file_')
    wf.connect(resample_o, 'new_fname', outputnode,'resampled_overlay')

    # node for axial montages
    montage_a = pe.MapNode(Function(input_names=['overlay',
                                                 'underlay',
                                                 'png_name',
                                                 'cbar_name'],
                                    output_names=['png_name'],
                                    function=montage_axial,
                                    as_module=True),
                           name='montage_a',
                           iterfield=['overlay'])
    montage_a.inputs.cbar_name = cbar_name
    montage_a.inputs.png_name = png_name + '_a.png'

    wf.connect(resample_u, 'new_fname', montage_a, 'underlay')

    wf.connect(resample_o, 'new_fname', montage_a, 'overlay')

    # node for sagittal montages
    montage_s = pe.MapNode(Function(input_names=['overlay',
                                                 'underlay',
                                                 'png_name',
                                                 'cbar_name'],
                                    output_names=['png_name'],
                                    function=montage_sagittal,
                                    as_module=True),
                           name='montage_s',
                           iterfield=['overlay'])
    montage_s.inputs.cbar_name = cbar_name
    montage_s.inputs.png_name = png_name + '_s.png'

    wf.connect(resample_u, 'new_fname', montage_s, 'underlay')
    wf.connect(resample_o, 'new_fname', montage_s, 'overlay')

    wf.connect(montage_a, 'png_name', outputnode, 'axial_png')
    wf.connect(montage_s, 'png_name', outputnode, 'sagittal_png')

    return wf


def create_montage_gm_wm_csf(wf_name, png_name):

    wf = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['underlay',
                                                       'overlay_csf',
                                                       'overlay_wm',
                                                       'overlay_gm']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['axial_png',
                                                        'sagittal_png',
                                                        'resampled_underlay',
                                                        'resampled_overlay_csf',
                                                        'resampled_overlay_wm',
                                                        'resampled_overlay_gm']),
                         name='outputspec')

    resample_u = pe.Node(Function(input_names=['file_'],
                                  output_names=['new_fname'],
                                  function=resample_1mm,
                                  as_module=True),
                         name='resample_u')

    resample_o_csf = resample_u.clone('resample_o_csf')
    resample_o_wm = resample_u.clone('resample_o_wm')
    resample_o_gm = resample_u.clone('resample_o_gm')

    wf.connect(inputNode, 'underlay', resample_u, 'file_')
    wf.connect(inputNode, 'overlay_csf', resample_o_csf, 'file_')
    wf.connect(inputNode, 'overlay_gm', resample_o_gm, 'file_')
    wf.connect(inputNode, 'overlay_wm', resample_o_wm, 'file_')

    montage_a = pe.Node(Function(input_names=['overlay_csf',
                                              'overlay_wm',
                                              'overlay_gm',
                                              'underlay',
                                              'png_name'],
                                 output_names=['png_name'],
                                 function=montage_gm_wm_csf_axial,
                                 as_module=True),
                        name='montage_a')


    wf.connect(resample_u, 'new_fname', montage_a, 'underlay')
    wf.connect(resample_o_csf, 'new_fname', montage_a, 'overlay_csf')
    wf.connect(resample_o_gm, 'new_fname', montage_a, 'overlay_gm')
    wf.connect(resample_o_wm, 'new_fname', montage_a, 'overlay_wm')
    montage_a.inputs.png_name = png_name + '_a.png'

    montage_s = pe.Node(Function(input_names=['overlay_csf',
                                              'overlay_wm',
                                              'overlay_gm',
                                              'underlay',
                                              'png_name'],
                                 output_names=['png_name'],
                                 function=montage_gm_wm_csf_sagittal,
                                 as_module=True),
                        name='montage_s')
    montage_s.inputs.png_name = png_name + '_s.png'

    wf.connect(resample_u, 'new_fname', montage_s, 'underlay')
    wf.connect(resample_o_csf, 'new_fname', montage_s, 'overlay_csf')
    wf.connect(resample_o_gm, 'new_fname', montage_s, 'overlay_gm')
    wf.connect(resample_o_wm, 'new_fname', montage_s, 'overlay_wm')

    wf.connect(resample_u, 'new_fname', outputNode, 'resampled_underlay')
    wf.connect(resample_o_csf, 'new_fname', outputNode, 'resampled_overlay_csf')
    wf.connect(resample_o_wm, 'new_fname', outputNode, 'resampled_overlay_wm')
    wf.connect(resample_o_gm, 'new_fname', outputNode, 'resampled_overlay_gm')
    wf.connect(montage_a, 'png_name', outputNode, 'axial_png')
    wf.connect(montage_s, 'png_name', outputNode, 'sagittal_png')

    return wf


def qa_montages(workflow, c, strat, num_strat,
                qc_montage_id_a, qc_montage_id_s, qc_hist_id,
                measure, idx):

    try:
        overlay, out_file = strat[measure]

        overlay_drop_percent = pe.MapNode(
            Function(input_names=['measure_file',
                                  'percent'],
                     output_names=['modified_measure_file'],
                     function=drop_percent,
                     as_module=True),
            name='dp_%s_%d' % (
                measure, num_strat),
            iterfield=['measure_file'])

        overlay_drop_percent.inputs.percent = 99.999

        workflow.connect(overlay, out_file,
                         overlay_drop_percent, 'measure_file')

        montage = create_montage('montage_%s_%d' % (measure, num_strat),
                                 'cyan_to_yellow', measure)
        
        node, out_file_template = strat['template_brain_for_func_derivative']
        workflow.connect(node, out_file_template, 
                        montage, 'inputspec.underlay')

        workflow.connect(overlay_drop_percent, 'modified_measure_file',
                         montage, 'inputspec.overlay')

        if 'centrality' in measure:
            histogram = pe.MapNode(
                Function(input_names=['measure_file',
                                      'measure'],
                         output_names=['hist_path'],
                         function=gen_histogram,
                         as_module=True),
                name='hist_{0}_{1}'.format(measure, num_strat),
                iterfield=['measure_file'])
        else:
            histogram = pe.Node(
                Function(input_names=['measure_file',
                                      'measure'],
                         output_names=['hist_path'],
                         function=gen_histogram,
                         as_module=True),
                name='hist_{0}_{1}'.format(measure, num_strat))

        histogram.inputs.measure = measure

        workflow.connect(overlay, out_file,
                         histogram, 'measure_file')

        strat.update_resource_pool({'qc___%s_a' % measure: (montage, 'outputspec.axial_png'),
                                    'qc___%s_s' % measure: (montage, 'outputspec.sagittal_png'),
                                    'qc___%s_hist' % measure: (histogram, 'hist_path')})

        if not idx in qc_montage_id_a:
            qc_montage_id_a[idx] = '%s_a' % measure
            qc_montage_id_s[idx] = '%s_s' % measure
            qc_hist_id[idx] = '%s_hist' % measure

    except Exception as e:
        print("[!] Connection of QA montages workflow for %s " \
                "has failed.\n" % measure)
        print("Error: %s" % e)
        pass


def create_qc_snr(wf_name='qc_snr'):

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['functional_preprocessed',
                                                        'functional_brain_mask',
                                                        'functional_to_anat_linear_xfm',
                                                        'anatomical_brain',
                                                        'mean_functional_in_anat']),
                        name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['snr_axial_image',
                                                         'snr_sagittal_image',
                                                         'snr_histogram_image',
                                                         'snr_mean']),
                          name='outputspec')

    std_dev = pe.Node(afni.TStat(args='-stdev'),
                      name='std_dev')

    std_dev.inputs.outputtype = 'NIFTI_GZ'
    wf.connect(input_node, 'functional_preprocessed', std_dev, 'in_file')
    wf.connect(input_node, 'functional_brain_mask', std_dev, 'mask')


    std_dev_anat = pe.Node(fsl.ApplyWarp(interp='trilinear'), 
                           name='std_dev_anat')
    wf.connect(input_node, 'functional_to_anat_linear_xfm', std_dev_anat, 'premat')
    wf.connect(std_dev, 'out_file', std_dev_anat, 'in_file')
    wf.connect(input_node, 'anatomical_brain', std_dev_anat, 'ref_file')

    snr = pe.Node(afni.Calc(expr='b/a'), name='snr')
    snr.inputs.outputtype = 'NIFTI_GZ'
    wf.connect(input_node, 'mean_functional_in_anat', snr, 'in_file_b')
    wf.connect(std_dev_anat, 'out_file', snr, 'in_file_a')

    snr_val = pe.Node(Function(input_names=['measure_file'],
                                output_names=['snr_storefl'],
                                function=cal_snr_val,
                                as_module=True),
                        name='snr_val')

    wf.connect(snr, 'out_file', snr_val, 'measure_file')

    hist_snr = pe.Node(Function(input_names=['measure_file', 'measure'],
                                output_names=['hist_path'],
                                function=gen_histogram,
                                as_module=True),
                        name='hist_snr')

    hist_snr.inputs.measure = 'snr'

    wf.connect(snr, 'out_file', hist_snr, 'measure_file')

    snr_drop_percent = pe.Node(
        Function(input_names=['measure_file',
                              'percent'],
                 output_names=['modified_measure_file'],
                 function=drop_percent,
                 as_module=True),
        name='dp_snr'
    )

    snr_drop_percent.inputs.percent = 99

    wf.connect(snr, 'out_file', snr_drop_percent, 'measure_file')

    montage_snr = create_montage('montage_snr',
                                 'red_to_blue',
                                 'snr')

    wf.connect(snr_drop_percent, 'modified_measure_file', montage_snr, 'inputspec.overlay')
    wf.connect(input_node, 'anatomical_brain', montage_snr, 'inputspec.underlay')

    wf.connect(montage_snr, 'outputspec.axial_png', output_node, 'snr_axial_image')
    wf.connect(montage_snr, 'outputspec.sagittal_png', output_node, 'snr_sagittal_image')
    wf.connect(hist_snr, 'hist_path', output_node, 'snr_histogram_image')
    wf.connect(snr_val, 'snr_storefl', output_node, 'snr_mean')

    return wf


def create_qc_motion(wf_name='qc_motion'):

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['motion_parameters']),
                        name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['motion_translation_plot',
                                                         'motion_rotation_plot']),
                          name='outputspec')

    mov_plot = pe.Node(Function(input_names=['motion_parameters'],
                                output_names=['translation_plot',
                                              'rotation_plot'],
                                function=gen_motion_plt,
                                as_module=True),
                       name='motion_plot')

    wf.connect(input_node, 'motion_parameters', mov_plot, 'motion_parameters')
    wf.connect(mov_plot, 'translation_plot', output_node, 'motion_translation_plot')
    wf.connect(mov_plot, 'rotation_plot', output_node, 'motion_rotation_plot')

    return wf


def create_qc_fd(wf_name='qc_fd'):

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['fd', 'excluded_volumes']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['fd_histogram_plot']),
                          name='outputspec')

    fd_plot = pe.Node(Function(input_names=['arr',
                                            'measure',
                                            'ex_vol'],
                               output_names=['hist_path'],
                               function=gen_plot_png,
                               as_module=True),
                      name='fd_plot')

    fd_plot.inputs.measure = 'FD'

    wf.connect(input_node, 'fd', fd_plot, 'arr')
    wf.connect(input_node, 'excluded_volumes', fd_plot, 'ex_vol')
    wf.connect(fd_plot, 'hist_path', output_node, 'fd_histogram_plot')

    return wf

def create_qc_carpet(wf_name='qc_carpet', output_image='qc_carpet'):

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['functional_to_standard',
                                                        'mean_functional_to_standard',
                                                        'anatomical_gm_mask',
                                                        'anatomical_wm_mask',
                                                        'anatomical_csf_mask']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['carpet_plot']),
                          name='outputspec')


    gm_resample = pe.Node(afni.Resample(), name='gm_resample')
    gm_resample.inputs.outputtype = 'NIFTI'
    wf.connect(input_node, 'anatomical_gm_mask', gm_resample, 'in_file')
    wf.connect(input_node, 'mean_functional_to_standard', gm_resample, 'master')

    gm_mask = pe.Node(afni.Calc(), name="gm_mask")
    gm_mask.inputs.expr = 'astep(a, 0.5)'
    gm_mask.inputs.outputtype = 'NIFTI'
    wf.connect(gm_resample, 'out_file', gm_mask, 'in_file_a')



    wm_resample = pe.Node(afni.Resample(), name='wm_resample')
    wm_resample.inputs.outputtype = 'NIFTI'
    wf.connect(input_node, 'anatomical_wm_mask', wm_resample, 'in_file')
    wf.connect(input_node, 'mean_functional_to_standard', wm_resample, 'master')

    wm_mask = pe.Node(afni.Calc(), name="wm_mask")
    wm_mask.inputs.expr = 'astep(a, 0.5)'
    wm_mask.inputs.outputtype = 'NIFTI'
    wf.connect(wm_resample, 'out_file', wm_mask, 'in_file_a')



    csf_resample = pe.Node(afni.Resample(), name='csf_resample')
    csf_resample.inputs.outputtype = 'NIFTI'
    wf.connect(input_node, 'anatomical_csf_mask', csf_resample, 'in_file')
    wf.connect(input_node, 'mean_functional_to_standard', csf_resample, 'master')

    csf_mask = pe.Node(afni.Calc(), name="csf_mask")
    csf_mask.inputs.expr = 'astep(a, 0.5)'
    csf_mask.inputs.outputtype = 'NIFTI'
    wf.connect(csf_resample, 'out_file', csf_mask, 'in_file_a')



    carpet_plot = pe.Node(Function(input_names=['gm_mask',
                                                'wm_mask',
                                                'csf_mask',
                                                'functional_to_standard',
                                                'output'],
                                   output_names=['carpet_plot'],
                                   function=gen_carpet_plt,
                                   as_module=True),
                          name='carpet_plot')

    carpet_plot.inputs.output = output_image
    wf.connect(gm_mask, 'out_file', carpet_plot, 'gm_mask')
    wf.connect(wm_mask, 'out_file', carpet_plot, 'wm_mask')
    wf.connect(csf_mask, 'out_file', carpet_plot, 'csf_mask')
    wf.connect(input_node, 'functional_to_standard', carpet_plot, 'functional_to_standard')
    wf.connect(carpet_plot, 'carpet_plot', output_node, 'carpet_plot')

    return wf


def afni_Edge3(in_file):
    import os
    import subprocess

    comm = "3dedge3 -input %s -prefix skull_edge.nii.gz" % in_file
    subprocess.getoutput(comm)

    return os.path.abspath('./skull_edge.nii.gz')


def create_qc_skullstrip(wf_name='qc_skullstrip'):

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['anatomical_brain',
                                                        'anatomical_reorient']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['axial_image',
                                                         'sagittal_image']),
                          name='outputspec')

    skull_edge = pe.Node(Function(input_names=['in_file'],
                                  output_names=['out_file'],
                                  function=afni_Edge3,
                                  as_module=True),
                         name='skull_edge')

    montage_skull = create_montage('montage_skull', 'red', 'skull_vis')

    wf.connect(input_node, 'anatomical_reorient', skull_edge, 'in_file')
    wf.connect(input_node, 'anatomical_brain', montage_skull, 'inputspec.underlay')
    wf.connect(skull_edge, 'out_file', montage_skull, 'inputspec.overlay')
    
    wf.connect(montage_skull, 'outputspec.axial_png', output_node, 'axial_image')
    wf.connect(montage_skull, 'outputspec.sagittal_png', output_node, 'sagittal_image')

    return wf
