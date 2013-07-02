import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def create_func_datasource(rest_dict, wf_name='func_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)


    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'scan'],
                                mandatory_inputs=True),
                        name='inputnode')
    inputnode.iterables = [('scan', rest_dict.keys())]

    selectrest = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['rest'],
                        function=get_rest),
                         name='selectrest')
    selectrest.inputs.rest_dict = rest_dict

    outputnode = pe.Node(util.IdentityInterface(fields=['subject',
                                                     'rest',
                                                     'scan' ]),
                         name='outputspec')

    wf.connect(inputnode, 'scan', selectrest, 'scan')

    wf.connect(inputnode, 'subject', outputnode, 'subject')
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



def create_mask_dataflow(dirPath, wf_name='datasource_mask'):

    import nipype.interfaces.io as nio
    import os

    wf = pe.Workflow(name=wf_name)
    masks = open(dirPath, 'r').readlines()

    mask_dict = {}
    for mask_file in masks:

        mask_file = mask_file.rstrip('\r\n')

        if mask_file.strip() == '' or mask_file.startswith('#'):
            continue

        base_file = os.path.basename(mask_file)
        base_name = ''
        if base_file.endswith('.nii'):
            base_name = os.path.splitext(base_file)[0]
        elif(base_file.endswith('.nii.gz')):
            base_name = os.path.splitext(os.path.splitext(base_file)[0])[0]
        else:
            raise("File extension not in  .nii and .nii.gz File: %s" % mask_file)

        if not (base_name in mask_dict):
            mask_dict[base_name] = mask_file
        else:
            raise ValueError('Files with same name not allowed %s %s' % (mask_file, mask_dict[base_name]))

    inputnode = pe.Node(util.IdentityInterface(
                            fields=['mask'],
                            mandatory_inputs=True),
                    name='inputspec')

    inputnode.iterables = [('mask', mask_dict.keys())]

    selectmask = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['out_file'],
                                       function=get_rest),
                         name='select_mask')
    selectmask.inputs.rest_dict = mask_dict

    outputnode = pe.Node(util.IdentityInterface(fields=['out_file']),
                         name='outputspec')

    wf.connect(inputnode, 'mask',
               selectmask, 'scan')

    wf.connect(selectmask, 'out_file',
               outputnode, 'out_file')
    return wf


def create_roi_dataflow(dirPath, wf_name='datasource_roi'):

    import nipype.interfaces.io as nio
    import os


    wf = pe.Workflow(name=wf_name)
    rois = open(dirPath, 'r').readlines()

    roi_dict = {}
    for roi_file in rois:

        roi_file = roi_file.rstrip('\r\n')

        if roi_file.strip() == '' or roi_file.startswith('#'):
            continue

        base_file = os.path.basename(roi_file)
        base_name = ''
        if base_file.endswith('.nii'):
            base_name = os.path.splitext(base_file)[0]
        elif(base_file.endswith('.nii.gz')):
            base_name = os.path.splitext(os.path.splitext(base_file)[0])[0]
        else:
            raise("File extension not in  .nii and .nii.gz File: %s" % roi_file)

        if not (base_name in roi_dict):
            roi_dict[base_name] = roi_file
        else:
            raise ValueError('Files with same name not allowed %s %s' % (roi_file, roi_dict[base_name]))

    inputnode = pe.Node(util.IdentityInterface(
                            fields=['roi'],
                            mandatory_inputs=True),
                    name='inputspec')

    inputnode.iterables = [('roi', roi_dict.keys())]

    selectroi = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['out_file'],
                                       function=get_rest),
                         name='select_roi')
    selectroi.inputs.rest_dict = roi_dict

    wf.connect(inputnode, 'roi',
               selectroi, 'scan')
    return wf


def create_spatial_map_dataflow(dirPath, wf_name='datasource_maps'):

    import nipype.interfaces.io as nio
    import os


    wf = pe.Workflow(name=wf_name)
    spatial_maps = open(dirPath, 'r').readlines()

    spatial_map_dict = {}
    for spatial_map_file in spatial_maps:

        spatial_map_file = spatial_map_file.rstrip('\r\n')
        base_file = os.path.basename(spatial_map_file)
        base_name = ''
        try:
            if base_file.endswith('.nii'):
                base_name = os.path.splitext(base_file)[0]
            elif(base_file.endswith('.nii.gz')):
                base_name = os.path.splitext(os.path.splitext(base_file)[0])[0]
            else:
                raise Exception("File extension not in  .nii and .nii.gz File: %s" % spatial_map_file)
        except Exception, e:
            print('error in spatial_map_dataflow: ', e)

        if not (base_name in spatial_map_dict):
            spatial_map_dict[base_name] = spatial_map_file
        else:
            raise ValueError('Files with same name not allowed %s %s' % (spatial_map_file, spatial_map_dict[base_name]))

    inputnode = pe.Node(util.IdentityInterface(
                            fields=['spatial_map'],
                            mandatory_inputs=True),
                    name='inputspec')

    inputnode.iterables = [('spatial_map', spatial_map_dict.keys())]

    select_spatial_map = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['out_file'],
                                       function=get_rest),
                         name='select_spatial_map')
    select_spatial_map.inputs.rest_dict = spatial_map_dict

    wf.connect(inputnode, 'spatial_map',
               select_spatial_map, 'scan')
    return wf


def create_grp_analysis_dataflow(wf_name='gp_dataflow'):

        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util
        from CPAC.utils import  select_model_files

        wf = pe.Workflow(name=wf_name)

        inputnode = pe.Node(util.IdentityInterface(fields=['ftest',
                                                           'grp_model'],
                                                   mandatory_inputs=True),
                            name='inputspec')

        selectmodel = pe.Node(util.Function(input_names=['model',
                                                         'ftest'],
                                           output_names=['fts_file',
                                                         'con_file',
                                                         'grp_file',
                                                         'mat_file'],
                                           function=select_model_files),
                             name='selectnode')

        wf.connect(inputnode, 'ftest',
                   selectmodel, 'ftest')
        wf.connect(inputnode, 'grp_model',
                   selectmodel, 'model')



        outputnode = pe.Node(util.IdentityInterface(fields=['fts',
                                                            'grp',
                                                            'mat',
                                                            'con'],
                                mandatory_inputs=True),
                    name='outputspec')

        wf.connect(selectmodel, 'mat_file',
                   outputnode, 'mat')
        wf.connect(selectmodel, 'grp_file',
                   outputnode, 'grp') 
        wf.connect(selectmodel, 'fts_file',
                   outputnode, 'fts')
        wf.connect(selectmodel, 'con_file',
                   outputnode, 'con')


        return wf

def create_gpa_dataflow(wf_name='gp_dataflow'):

        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util
        from CPAC.utils import modify_model, select_model_files

        wf = pe.Workflow(name=wf_name)

        inputnode = pe.Node(util.IdentityInterface(fields=['ftest',
                                                           'grp_model',
                                                           'input_sublist',
                                                           'output_sublist'],
                                                   mandatory_inputs=True),
                            name='inputspec')

        selectmodel = pe.Node(util.Function(input_names=['model',
                                                         'ftest'],
                                           output_names=['fts_file',
                                                         'con_file',
                                                         'grp_file',
                                                         'mat_file'],
                                           function=select_model_files),
                             name='selectnode')

        wf.connect(inputnode, 'ftest',
                   selectmodel, 'ftest')
        wf.connect(inputnode, 'grp_model',
                   selectmodel, 'model')

        modifymodel = pe.Node(util.Function(input_names=['input_sublist',
                                                            'output_sublist',
                                                            'mat_file',
                                                            'grp_file'],
                                             output_names=['grp_file',
                                                             'mat_file',
                                                             'sub_file'],
                                             function=modify_model),
                              name='modifymodel')

        wf.connect(selectmodel, 'mat_file',
                   modifymodel, 'mat_file')
        wf.connect(selectmodel, 'grp_file',
                   modifymodel, 'grp_file')
        wf.connect(inputnode, 'input_sublist',
                   modifymodel, 'input_sublist')
        wf.connect(inputnode, 'output_sublist',
                   modifymodel, 'output_sublist')

        outputnode = pe.Node(util.IdentityInterface(fields=['fts',
                                                            'grp',
                                                            'mat',
                                                            'con',
                                                            'sublist'],
                                mandatory_inputs=True),
                    name='outputspec')

        wf.connect(modifymodel, 'mat_file',
                   outputnode, 'mat')
        wf.connect(modifymodel, 'grp_file',
                   outputnode, 'grp')
        wf.connect(modifymodel, 'sub_file',
                   outputnode, 'sublist')
        wf.connect(selectmodel, 'fts_file',
                   outputnode, 'fts')
        wf.connect(selectmodel, 'con_file',
                   outputnode, 'con')


        return wf

