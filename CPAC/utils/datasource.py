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



def create_anat_datasource(wf_name='anat_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'anat'],
                                mandatory_inputs=True),
                        name='inputnode')

    outputnode = pe.Node(util.IdentityInterface(fields=['subject',
                                                     'anat' ]),
                         name='outputspec')

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(inputnode, 'anat', outputnode, 'anat')

    return wf



def create_mask_dataflow(dirPath, wf_name = 'datasource_mask'):

    import nipype.interfaces.io as nio
    import os

    masklist = [os.path.splitext(
                    os.path.splitext(f)[0])[0]
                    for f in os.listdir(dirPath)]

    datasource = pe.Node(interface=nio.DataGrabber(
                                        infields=['mask'],
                                        outfields=['out_file']),
                         name = wf_name)
    datasource.inputs.base_directory = dirPath
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(out_file='%s.nii.gz')
    datasource.inputs.template_args = dict(out_file=[['mask']])
    datasource.iterables = ('mask', masklist)

    return datasource


def create_roi_dataflow(dirPath, wf_name = 'datasource_roi'):

    import nipype.interfaces.io as nio
    import os

    unitlist = [os.path.splitext(
                    os.path.splitext(f)[0])[0]
                    for f in os.listdir(dirPath)]

    datasource = pe.Node(interface = nio.DataGrabber(
                                        infields=['roi'],
                                        outfields=['out_file']),
                         name = wf_name)
    datasource.inputs.base_directory = dirPath
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(out_file='%s.nii.gz')
    datasource.inputs.template_args = dict(out_file=[['roi']])
    datasource.iterables = ('roi', unitlist)

    return datasource

def create_gpa_dataflow(model_dict, ftest, wf_name = 'gp_dataflow'):
        """
        Dataflow to iterate over each model and 
        pick the model files for group analysis
        """
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util

        wf = pe.Workflow(name=wf_name) 
        
        inputnode = pe.Node(util.IdentityInterface(
                                fields=['model'],
                                mandatory_inputs=True),
                        name='inputspec')
        
        inputnode.iterables = [('model', model_dict.keys())]
        
        outputnode = pe.Node(util.Function(input_names=['model',
                                                        'model_dict', 
                                                        'ftest'],
                                           output_names=['fts', 
                                                         'con', 
                                                         'grp', 
                                                         'mat'],
                                           function = get_model),
                             name  = 'outputspec')
        outputnode.inputs.model_dict = model_dict
        outputnode.inputs.ftest = ftest
        
        wf.connect(inputnode, 'model', 
                   outputnode, 'model')
        
        return wf
    
    
    
def get_model(model, model_dict, ftest):
    
    try:
        files = model_dict[model]
        
        if ftest:
            fts_file = [file for file in files if file.endswith('.fts')][0]
        else:
            fts_file =''
    
        mat_file = [file for file in files if file.endswith('.mat')][0]
        grp_file = [file for file in files if file.endswith('.grp')][0]
        con_file = [file for file in files if file.endswith('.con')][0]
    
    except Exception:
        print "All the model files are not present. Please check the model folder"
        raise
    
    return fts_file, con_file, grp_file, mat_file

    
