import os

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl

def fsl_apply_transform_func_to_mni(workflow, output_name, func_key, ref_key, num_strat, strat, interpolation_method, map_node=False):

    strat_nodes = strat.get_nodes_names()

    if map_node == True:
        # func_mni_warp
        func_mni_warp = pe.MapNode(interface=fsl.ApplyWarp(),
                name='func_mni_fsl_warp_{0}_{1:d}'.format(output_name, num_strat),
                iterfield=['in_file'],
                mem_gb=1.5)
    else:
        # func_mni_warp
        func_mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                name='func_mni_fsl_warp_{0}_{1:d}'.format(output_name, num_strat))

        
    func_mni_warp.inputs.interp = interpolation_method

    if func_key == 'leaf':
        func_node, func_file = strat.get_leaf_properties()
    else:
        func_node, func_file = strat[func_key]

    workflow.connect(func_node, func_file,
                     func_mni_warp, 'in_file')

    ref_node, ref_out_file = strat['template_brain_for_func_preproc']

    ref_node, ref_out_file = strat[ref_key]
    workflow.connect(ref_node, ref_out_file,
                     func_mni_warp, 'ref_file')

    if 'anat_mni_fnirt_register' in strat_nodes:

        node, out_file = strat['functional_to_anat_linear_xfm']
        workflow.connect(node, out_file,
                         func_mni_warp, 'premat')

        node, out_file = strat['anatomical_to_mni_nonlinear_xfm']
        workflow.connect(node, out_file,
                         func_mni_warp, 'field_file')

    elif 'anat_mni_flirt_register' in strat_nodes:

        if 'functional_to_mni_linear_xfm' not in strat:

            combine_transforms = pe.Node(interface=fsl.ConvertXFM(),
                name='combine_fsl_xforms_func_mni_{0}_{1:d}'.format(output_name, num_strat))

            combine_transforms.inputs.concat_xfm = True

            node, out_file = strat['anatomical_to_mni_linear_xfm']
            workflow.connect(node, out_file,
                             combine_transforms, 'in_file2')

            node, out_file = strat['functional_to_anat_linear_xfm']
            workflow.connect(node, out_file,
                             combine_transforms, 'in_file')

            strat.update_resource_pool({ 'functional_to_mni_linear_xfm': (combine_transforms, 'out_file')}) 
            strat.append_name(combine_transforms.name)

        combine_transforms, outfile = strat['functional_to_mni_linear_xfm']

        workflow.connect(combine_transforms, outfile,
                         func_mni_warp, 'premat')

        strat.update_resource_pool({ output_name: (func_mni_warp, 'out_file')}) 
        strat.append_name(func_mni_warp.name)

    else:
        raise ValueError('Could not find anat_mni_flirt_register or anat_mni_fnirt register in nodes')
    
    return workflow
