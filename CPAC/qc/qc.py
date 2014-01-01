import os
import sys
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from CPAC.qc.qc import *
from CPAC.qc.utils import *


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


    resample_u = pe.Node(util.Function(input_names=['file_'],
                                       output_names=['new_fname'],
                                       function=resample_1mm),
                           name='resample_u')

    resample_o = resample_u.clone('resample_o')



    montage_a = pe.Node(util.Function(input_names=['overlay',
                                                  'underlay',
                                                  'png_name',
                                                  'cbar_name'],
                                      output_names=['png_name'],
                                      function=montage_axial),
                         name='montage_a')
    montage_a.inputs.cbar_name = cbar_name
    montage_a.inputs.png_name = png_name + '_a.png'


    montage_s = pe.Node(util.Function(input_names=['overlay',
                                                   'underlay',
                                                   'png_name',
                                                   'cbar_name'],
                                      output_names=['png_name'],
                                      function=montage_sagittal),
                                   name='montage_s')
    montage_s.inputs.cbar_name = cbar_name
    montage_s.inputs.png_name = png_name + '_s.png'

    wf.connect(inputNode, 'underlay',
                resample_u, 'file_')

    wf.connect(inputNode, 'overlay',
                resample_o, 'file_')

    wf.connect(resample_u, 'new_fname',
                montage_a, 'underlay')

    wf.connect(resample_o, 'new_fname',
                montage_a, 'overlay')

    wf.connect(resample_u, 'new_fname',
                montage_s, 'underlay')

    wf.connect(resample_o, 'new_fname',
                montage_s, 'overlay')

    wf.connect(resample_u, 'new_fname',
                outputNode, 'resampled_underlay')

    wf.connect(resample_o, 'new_fname',
                outputNode, 'resampled_overlay')

    wf.connect(montage_a, 'png_name',
                outputNode, 'axial_png')

    wf.connect(montage_s, 'png_name',
                outputNode, 'sagittal_png')


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


    resample_u = pe.Node(util.Function(input_names=['file_'],
                                       output_names=['new_fname'],
                                       function=resample_1mm),
                           name='resample_u')

    resample_o_csf = resample_u.clone('resample_o_csf')
    resample_o_wm = resample_u.clone('resample_o_wm')
    resample_o_gm = resample_u.clone('resample_o_gm')



    montage_a = pe.Node(util.Function(input_names=['overlay_csf',
                                                   'overlay_wm',
                                                   'overlay_gm',
                                                  'underlay',
                                                  'png_name'],
                                      output_names=['png_name'],
                                      function=montage_gm_wm_csf_axial),
                         name='montage_a')
    montage_a.inputs.png_name = png_name + '_a.png'


    montage_s = pe.Node(util.Function(input_names=['overlay_csf',
                                                   'overlay_wm',
                                                   'overlay_gm',
                                                   'underlay',
                                                   'png_name'],
                                      output_names=['png_name'],
                                      function=montage_gm_wm_csf_sagittal),
                                   name='montage_s')
    montage_s.inputs.png_name = png_name + '_s.png'

    wf.connect(inputNode, 'underlay',
                resample_u, 'file_')

    wf.connect(inputNode, 'overlay_csf',
                resample_o_csf, 'file_')

    wf.connect(inputNode, 'overlay_gm',
                resample_o_gm, 'file_')

    wf.connect(inputNode, 'overlay_wm',
                resample_o_wm, 'file_')

    wf.connect(resample_u, 'new_fname',
                montage_a, 'underlay')

    wf.connect(resample_o_csf, 'new_fname',
                montage_a, 'overlay_csf')

    wf.connect(resample_o_gm, 'new_fname',
                montage_a, 'overlay_gm')

    wf.connect(resample_o_wm, 'new_fname',
                montage_a, 'overlay_wm')

    wf.connect(resample_u, 'new_fname',
                montage_s, 'underlay')

    wf.connect(resample_o_csf, 'new_fname',
                montage_s, 'overlay_csf')

    wf.connect(resample_o_gm, 'new_fname',
                montage_s, 'overlay_gm')

    wf.connect(resample_o_wm, 'new_fname',
                montage_s, 'overlay_wm')

    wf.connect(resample_u, 'new_fname',
                outputNode, 'resampled_underlay')

    wf.connect(resample_o_csf, 'new_fname',
                outputNode, 'resampled_overlay_csf')

    wf.connect(resample_o_wm, 'new_fname',
                outputNode, 'resampled_overlay_wm')

    wf.connect(resample_o_gm, 'new_fname',
                outputNode, 'resampled_overlay_gm')

    wf.connect(montage_a, 'png_name',
                outputNode, 'axial_png')

    wf.connect(montage_s, 'png_name',
                outputNode, 'sagittal_png')


    return wf
