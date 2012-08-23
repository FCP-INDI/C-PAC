import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def create_func_datasource(rest_dict, wf_name = 'func_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)
    

    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'scan'],
                                mandatory_inputs=True),
                        name='inputnode')
    inputnode.iterables = [('scan',  rest_dict.keys())]
    
    selectrest = pe.Node(util.Function(input_names=['scan','rest_dict'],
                                       output_names=['rest'],
                        function = get_rest ),
                         name  = 'selectrest')
    selectrest.inputs.rest_dict = rest_dict
    
    outputnode = pe.Node(util.IdentityInterface(fields=['subject',
                                                     'rest',
                                                     'scan' ]),
                         name='outputspec')
    
    wf.connect(inputnode, 'scan', selectrest, 'scan')
    
    wf.connect(inputnode,  'subject', outputnode, 'subject')
    wf.connect(selectrest, 'rest', outputnode, 'rest')
    wf.connect(inputnode, 'scan', outputnode, 'scan')

    return wf



def get_rest(scan, rest_dict):
    return rest_dict[scan]
  
