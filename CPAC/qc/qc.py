
from CPAC.utils.function import Function
from CPAC.qc.utils import resample_1mm, montage_axial, montage_sagittal, \
                          montage_gm_wm_csf_axial, montage_gm_wm_csf_sagittal

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def create_montage(wf_name, cbar_name, png_name):

    wf = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['underlay',
                                                       'overlay']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['axial_png',
                                                        'sagittal_png',
                                                        'resampled_underlay',
                                                        'resampled_overlay']),
                         name='outputspec')

    # node for resampling images to 1mm for QC pages
    resample_u = pe.Node(Function(input_names=['file_'],
                                  output_names=['new_fname'],
                                  function=resample_1mm,
                                  as_module=True),
                         name='resample_u')
    wf.connect(inputNode, 'underlay', resample_u, 'file_')

    # same for overlays (resampling to 1mm)
    resample_o = resample_u.clone('resample_o')
    wf.connect(inputNode, 'overlay', resample_o, 'file_')

    # node for axial montages
    montage_a = pe.Node(Function(input_names=['overlay',
                                              'underlay',
                                              'png_name',
                                              'cbar_name'],
                                 output_names=['png_name'],
                                 function=montage_axial,
                                 as_module=True),
                        name='montage_a')
    montage_a.inputs.cbar_name = cbar_name
    montage_a.inputs.png_name = png_name + '_a.png'

    wf.connect(resample_u, 'new_fname', montage_a, 'underlay')
    wf.connect(resample_o, 'new_fname', montage_a, 'overlay')

    # node for sagittal montages
    montage_s = pe.Node(Function(input_names=['overlay',
                                              'underlay',
                                              'png_name',
                                              'cbar_name'],
                                 output_names=['png_name'],
                                 function=montage_sagittal,
                                 as_module=True),
                        name='montage_s')
    montage_s.inputs.cbar_name = cbar_name
    montage_s.inputs.png_name = png_name + '_s.png'

    wf.connect(resample_u, 'new_fname', montage_s, 'underlay')
    wf.connect(resample_o, 'new_fname', montage_s, 'overlay')

    wf.connect(resample_u, 'new_fname', outputNode, 'resampled_underlay')
    wf.connect(resample_o, 'new_fname', outputNode, 'resampled_overlay')
    wf.connect(montage_a, 'png_name', outputNode, 'axial_png')
    wf.connect(montage_s, 'png_name', outputNode, 'sagittal_png')

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
    montage_a.inputs.png_name = png_name + '_a.png'

    wf.connect(resample_u, 'new_fname', montage_a, 'underlay')
    wf.connect(resample_o_csf, 'new_fname', montage_a, 'overlay_csf')
    wf.connect(resample_o_gm, 'new_fname', montage_a, 'overlay_gm')
    wf.connect(resample_o_wm, 'new_fname', montage_a, 'overlay_wm')

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
    wf.connect(resample_o_csf, 'new_fname',
               outputNode, 'resampled_overlay_csf')
    wf.connect(resample_o_wm, 'new_fname', outputNode, 'resampled_overlay_wm')
    wf.connect(resample_o_gm, 'new_fname', outputNode, 'resampled_overlay_gm')
    wf.connect(montage_a, 'png_name', outputNode, 'axial_png')
    wf.connect(montage_s, 'png_name', outputNode, 'sagittal_png')

    return wf


def QA_montages(measure, idx):
    try:
        overlay, out_file = strat.get_node_from_resource_pool(measure)

        overlay_drop_percent = pe.MapNode(function.Function(input_names=['measure_file',
                                                                            'percent'],
                                                            output_names=[
                                                                'modified_measure_file'],
                                                            function=drop_percent,
                                                            as_module=True),
                                            name='dp_%s_%d' % (
                                                measure, num_strat),
                                            iterfield=['measure_file'])
        overlay_drop_percent.inputs.percent = 99.999

        workflow.connect(overlay, out_file,
                            overlay_drop_percent, 'measure_file')

        montage = create_montage('montage_%s_%d' % (measure, num_strat), 'cyan_to_yellow', measure)
        montage.inputs.inputspec.underlay = c.template_brain_only_for_func

        workflow.connect(overlay_drop_percent, 'modified_measure_file',
                            montage, 'inputspec.overlay')

        if 'centrality' in measure:
            histogram = pe.MapNode(
                util.Function(input_names=['measure_file',
                                            'measure'],
                                output_names=['hist_path'],
                                function=gen_histogram),
                name='hist_{0}_{1}'.format(measure, num_strat),
                iterfield=['measure_file'])
        else:
            histogram = pe.Node(
                util.Function(input_names=['measure_file',
                                            'measure'],
                                output_names=['hist_path'],
                                function=gen_histogram),
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
        print "[!] Connection of QA montages workflow for %s " \
                "has failed.\n" % measure
        print "Error: %s" % e
        pass