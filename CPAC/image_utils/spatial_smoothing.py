# Copyright (C) 2018-2022  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
from nipype.interfaces import fsl, utility as util
from nipype.interfaces.afni import preprocess as afni
from CPAC.pipeline import nipype_pipeline_engine as pe


def set_gauss(fwhm):
    """
    Compute the sigma value, given Full Width Half Max.
    Returns an operand string

    Parameters
    ----------

    fwhm : float

    Returns
    -------

    op_string : string

    """

    sigma = float(fwhm) / 2.3548

    op_string = "-kernel gauss %f -fmean -mas " % sigma + "%s"

    return op_string


def spatial_smoothing(wf_name, fwhm, input_image_type='func_derivative',
                      opt=None):

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['in_file',
                                                       'mask']),
                        name='inputspec')

    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input')
    inputnode_fwhm.iterables = ("fwhm", fwhm)

    image_types = ['func_derivative', 'func_derivative_multi', 'func_4d',
                   'func_mask']

    if input_image_type not in image_types:
        raise ValueError('Input image type {0} should be one of '
                         '{1}'.format(input_image_type,
                                      ', '.join(image_types)))

    if opt == 'FSL':
        output_smooth_mem_gb = 4.0
        if input_image_type == 'func_derivative_multi':
            output_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                       name='smooth_multi',
                                       iterfield=['in_file'],
                                       mem_gb=output_smooth_mem_gb)
        else:
            output_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                                    name='smooth',
                                    mem_gb=output_smooth_mem_gb)

    elif opt == 'AFNI':
        if input_image_type == 'func_derivative_multi':
            output_smooth = pe.MapNode(interface=afni.BlurToFWHM(),
                                       name='smooth_multi',
                                       iterfield=['in_file'])
        else:
            output_smooth = pe.Node(interface=afni.BlurToFWHM(),
                                    name='smooth',
                                    iterfield=['in_file'])
        output_smooth.inputs.outputtype = 'NIFTI_GZ'

    if opt == 'FSL':
        # wire in the resource to be smoothed
        wf.connect(inputnode, 'in_file', output_smooth, 'in_file')
        # get the parameters for fwhm
        wf.connect(inputnode_fwhm, ('fwhm', set_gauss),
                         output_smooth, 'op_string')
        wf.connect(inputnode, 'mask', output_smooth, 'operand_files')
    elif opt =='AFNI':
        wf.connect(inputnode, 'in_file', output_smooth, 'in_file')
        wf.connect(inputnode_fwhm, 'fwhm', output_smooth, 'fwhm')
        wf.connect(inputnode, 'mask', output_smooth, 'mask')

    outputnode = pe.Node(util.IdentityInterface(fields=['out_file',
                                                        'fwhm']),
                         name='outputspec')

    wf.connect(output_smooth, 'out_file', outputnode, 'out_file')
    wf.connect(inputnode_fwhm, 'fwhm', outputnode, 'fwhm')

    return wf
