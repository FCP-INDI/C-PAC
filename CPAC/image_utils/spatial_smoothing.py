import nipype.interfaces.fsl as fsl
from nipype.interfaces.afni import preprocess as afni
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from CPAC.utils import Outputs


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


def spatial_smooth_outputs(workflow, func_key, strat, num_strat,
                           pipeline_config_object):

    # set the mask for the smoothing operation
    if "centrality" in func_key:
        smoothing_mask_key = (pipeline_config_object.templateSpecificationFile, 'local_path')
    elif func_key in Outputs.functional_timeseries or "sca_tempreg" in func_key:
        smoothing_mask_key = "functional_brain_mask_to_standard"
    elif func_key in Outputs.template_nonsmooth + Outputs.template_nonsmooth_mult:
        smoothing_mask_key = "functional_brain_mask_to_standard_derivative"
    else:
        smoothing_mask_key = "functional_brain_mask"

    # set the image type for the smoothing operation
    if func_key in Outputs.template_nonsmooth_mult + Outputs.native_nonsmooth_mult:
        input_image_type = 'func_derivative_multi'
    else:
        input_image_type = 'func_derivative'

    output_name = '{0}_smooth'.format(func_key)

    # insert the smoothing workflow
    if output_name not in strat:
        workflow, strat = spatial_smooth(workflow, func_key,
                                         smoothing_mask_key, output_name,
                                         strat, num_strat,
                                         pipeline_config_object,
                                         input_image_type=input_image_type)
    else:
        print('{0} already exists in the resource pool'.format(output_name))

    return workflow, strat


def spatial_smooth(workflow, func_key, mask_key, output_name, strat,
                   num_strat, pipeline_config_object,
                   input_image_type='func_derivative'):

    method = pipeline_config_object.post_processing['spatial_smoothing']['smoothing_method']
    inputnode_fwhm = pe.Node(util.IdentityInterface(fields=['fwhm']),
                             name='fwhm_input_{0}_{1}'.format(output_name,
                                                              num_strat))
    inputnode_fwhm.iterables = ("fwhm", pipeline_config_object.post_processing['spatial_smoothing']['fwhm'])

    image_types = ['func_derivative', 'func_derivative_multi', 'func_4d',
                   'func_mask']

    if input_image_type not in image_types:
        raise ValueError('Input image type {0} should be one of {1}'.format(input_image_type,
            ', '.join(image_types)))

    if method == 'FSL':
        if input_image_type == 'func_derivative_multi':
            output_smooth = pe.MapNode(interface=fsl.MultiImageMaths(),
                                       name='{0}_multi_{1}'.format(output_name, num_strat),
                                       iterfield=['in_file'])
        else:
            output_smooth = pe.Node(interface=fsl.MultiImageMaths(),
                                    name='{0}_{1}'.format(output_name,num_strat))

    elif method == 'AFNI':
        if input_image_type == 'func_derivative_multi':
            output_smooth = pe.MapNode(interface= afni.BlurToFWHM(),
                                       name='{0}_multi_{1}'.format(output_name, num_strat),
                                       iterfield=['in_file'])
        else:
            output_smooth = pe.Node(interface= afni.BlurToFWHM(),
                                    name='{0}_multi_{1}'.format(output_name, num_strat),
                                    iterfield=['in_file'])
        output_smooth.inputs.outputtype = 'NIFTI_GZ'

    if isinstance(func_key, str):
        if func_key == 'leaf':
            func_node, func_file = strat.get_leaf_properties()
        else:
            try:
                func_node, func_file = strat[func_key]
            except KeyError as e:
                print('Could not find func_key {0} in resource pool'.format(func_key))

    elif isinstance(func_key, tuple):
        func_node, func_file = func_key

    if isinstance(mask_key, str):
        mask_node, mask_file = strat[mask_key]
    elif isinstance(mask_key, tuple):
        mask_node, mask_file = mask_key
    else:
        raise ValueError('mask {0} ({1}) could not be deciphered'.format(mask_key, type(mask_key)))

    if method =='FSL':
        # TODO review connetion to config, is the node really necessary?
         # wire in the resource to be smoothed
        workflow.connect(func_node, func_file, output_smooth, 'in_file')
        # get the parameters for fwhm
        workflow.connect(inputnode_fwhm, ('fwhm', set_gauss),
                         output_smooth, 'op_string')
        workflow.connect(mask_node, mask_file,
                         output_smooth, 'operand_files')
    elif method =='AFNI':
        workflow.connect(func_node, func_file, output_smooth, 'in_file')
        workflow.connect(inputnode_fwhm, 'fwhm', output_smooth, 'fwhm')
        workflow.connect(mask_node, mask_file, output_smooth, 'mask')

    # output the resource
    strat.append_name(output_smooth.name)
    strat.update_resource_pool({output_name: (output_smooth, 'out_file')})

    return workflow, strat

