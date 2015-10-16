# CPAC/utils/datasource.py
#
#

'''
This module contains classes and functions used to interface with data
access
'''

# Import packages
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def create_func_datasource(rest_dict, wf_name='func_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)


    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'scan', 'creds_path'],
                                mandatory_inputs=True),
                        name='inputnode')
    inputnode.iterables = [('scan', rest_dict.keys())]

    selectrest = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['rest'],
                        function=get_rest),
                         name='selectrest')
    selectrest.inputs.rest_dict = rest_dict

    check_s3_node = pe.Node(util.Function(input_names=['file_path', 'creds_path'],
                                          output_names=['local_path'],
                                          function=check_for_s3),
                            name='check_for_s3')
    wf.connect(selectrest, 'rest', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')

    outputnode = pe.Node(util.IdentityInterface(fields=['subject',
                                                     'rest',
                                                     'scan' ]),
                         name='outputspec')

    wf.connect(inputnode, 'scan', selectrest, 'scan')

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(check_s3_node, 'local_path', outputnode, 'rest')
    wf.connect(inputnode, 'scan', outputnode, 'scan')

    return wf


def get_rest(scan, rest_dict):
    return rest_dict[scan]


# Check if passed in file is on S3
def check_for_s3(file_path, creds_path):
    '''
    '''

    # Import packages
    import os
    from CPAC.AWS import fetch_creds

    # Init variables
    s3_str = 's3://'
    local_download_dir = '/tmp'

    # Check for s3 string in filepaths
    if s3_str in file_path:
        # Get bucket name and bucket object
        bucket_name = file_path.replace(s3_str, '').split('/')[0]
        bucket = fetch_creds.return_bucket(creds_path, bucket_name)

        # Extract relative key path from bucket and local path
        s3_prefix = os.path.join(s3_str, bucket_name)
        rel_path = file_path.replace(s3_prefix, '').lstrip('/')
        local_path = file_path.replace(s3_prefix, local_download_dir)

        # Download file
        bucket.download_file(Key=rel_path, Filename=local_path)
    # Otherwise just return what was passed in
    else:
        local_path = file_path

    # Return the local path
    return local_path


# Anatomical datasource
def create_anat_datasource(wf_name='anat_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'anat', 'creds_path'],
                                mandatory_inputs=True),
                        name='inputnode')

    check_s3_node = pe.Node(util.Function(input_names=['file_path', 'creds_path'],
                                          output_names=['local_path'],
                                          function=check_for_s3),
                            name='check_for_s3')
    wf.connect(inputnode, 'anat', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')

    outputnode = pe.Node(util.IdentityInterface(fields=['subject',
                                                     'anat' ]),
                         name='outputspec')

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(check_s3_node, 'local_path', outputnode, 'anat')

    # Return the workflow
    return wf


# ROI mask dataflow
def create_roi_mask_dataflow(dir_path, mask_type, wf_name='datasource_roi_mask'):

    import nipype.interfaces.io as nio
    import os

    wf = pe.Workflow(name=wf_name)

    if mask_type == 'roi':
        tab = 'ROI Average TSE'
    elif mask_type == 'voxel':
        tab = 'ROI Voxelwise TSE'
    elif mask_type == 'centrality':
        tab = 'Network Centrality'


    if '.nii' in dir_path:

        masks = []
        masks.append(dir_path)

    elif '.txt' in dir_path:
        
        masks = open(dir_path, 'r').readlines()

    else:

        print '\n\n[!] CPAC says: Your ROI/mask specification file (under ' \
              '%s options) either needs to be a NIFTI file (.nii or ' \
              '.nii.gz) of an ROI/mask or a text file (.txt) containing a ' \
              'list of NIFTI files of ROI/mask files.\nPlease change this ' \
              'in your pipeline configuration file and try again.\n\n' % tab
        raise Exception


    mask_dict = {}

    for mask_file in masks:

        mask_file = mask_file.rstrip('\r\n')

        if not os.path.exists(mask_file):
            err = '\n\n[!] CPAC says: One of your ROI/mask specification ' \
                  'files (under %s options) does not have a correct path ' \
                  'or does not exist.\nTip: If all the paths are okay, ' \
                  'then ensure there are no whitespaces or blank lines in ' \
                  'your ROI specification file.\n\n' % mask_type
            raise Exception(err)

        if mask_file.strip() == '' or mask_file.startswith('#'):
            continue

        base_file = os.path.basename(mask_file)
        base_name = ''
        if base_file.endswith('.nii'):
            base_name = os.path.splitext(base_file)[0]
        elif(base_file.endswith('.nii.gz')):
            base_name = os.path.splitext(os.path.splitext(base_file)[0])[0]
        else:
            err = "\n\n[!] CPAC says: One of your ROI/mask specification " \
                  "files (under %s options) does not have '.nii' or " \
                  "'.nii.gz' as an extension.\n\nMask file: %s\n\n" \
                  % (tab, mask_file)
            raise Exception(err)

        if not (base_name in mask_dict):
            mask_dict[base_name] = mask_file
        else:
            err = "\n\n[!] CPAC says: You have two or more ROI/mask files " \
            "with the same name - please make sure these files are named " \
            "differently.\n\nDuplicate name: %s\n\nNote: This can be " \
            "changed in the ROI/mask file you specified under the %s " \
            "options.\n\n" % (mask_file, tab)
            raise Exception(err)


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




def create_spatial_map_dataflow(dirPath, wf_name='datasource_maps'):

    import nipype.interfaces.io as nio
    import os


    wf = pe.Workflow(name=wf_name)
    spatial_maps = open(dirPath, 'r').readlines()

    spatial_map_dict = {}
    
    for spatial_map_file in spatial_maps:

        spatial_map_file = spatial_map_file.rstrip('\r\n')

        if not os.path.exists(spatial_map_file):
            print "\n\n" + "ERROR: One of your spatial map files (under Spatial" + \
            " Regression options) does not have a correct path or does not exist." + \
            "\n" + "Tip: If all the paths are okay, then ensure there are no" + \
            " whitespaces or blank lines in your spatial map specification file." + \
            "\n\n" + "Error name: datasource_0001" + "\n\n"
            raise Exception

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
        from CPAC.utils import select_model_files

        wf = pe.Workflow(name=wf_name)

        inputnode = pe.Node(util.IdentityInterface(fields=['ftest',
                                                           'grp_model',
                                                           'model_name'],
                                                   mandatory_inputs=True),
                            name='inputspec')

        selectmodel = pe.Node(util.Function(input_names=['model',
                                                         'ftest',
                                                         'model_name'],
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
        wf.connect(inputnode, 'model_name', selectmodel, 'model_name')



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
