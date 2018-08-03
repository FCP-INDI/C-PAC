import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def get_rest(scan, rest_dict, resource="scan"):
    """Return the file path of the chosen resource stored in the functional
    file dictionary, if it exists.

    scan: the scan/series name or label
    rest_dict: the dictionary read in from the data configuration YAML file
               (sublist) nested under 'func:'
    resource: the dictionary key
                  scan - the functional timeseries
                  scan_parameters - path to the scan parameters JSON file, or
                                    a dictionary containing scan parameters
                                    information (to be phased out in the
                                    future)
    """
    try:
        file_path = rest_dict[scan][resource]
    except KeyError:
        file_path = None
    return file_path


def extract_scan_params_dct(scan_params_dct):
    return scan_params_dct


def get_map(map, map_dct):
    # return the spatial map required
    return map_dct[map]


def check_func_scan(func_scan_dct, scan):
    """Run some checks on the functional timeseries-related files for a given
    series/scan name or label."""

    scan_resources = func_scan_dct[scan]

    try:
        keys = scan_resources.keys()
    except AttributeError:
        err = "\n[!] The data configuration file you provided is " \
              "missing a level under the 'func:' key. CPAC versions " \
              "1.2 and later use data configurations with an " \
              "additional level of nesting.\n\nExample\nfunc:\n  " \
              "rest01:\n    scan: /path/to/rest01_func.nii.gz\n" \
              "    scan parameters: /path/to/scan_params.json\n\n" \
              "See the User Guide for more information.\n\n"
        raise Exception(err)

    # actual 4D time series file
    if "scan" not in scan_resources.keys():
        err = "\n\n[!] The {0} scan is missing its actual time-series " \
              "scan file, which should be a filepath labeled with the " \
              "'scan' key.\n\n".format(scan)
        raise Exception(err)

    # Nipype restriction (may have changed)
    if '.' in scan or '+' in scan or '*' in scan:
        raise Exception('\n\n[!] Scan names cannot contain any special '
                        'characters (., +, *, etc.). Please update this '
                        'and try again.\n\nScan: {0}'
                        '\n\n'.format(scan))


def create_func_datasource(rest_dict, wf_name='func_datasource'):
    """Return the functional timeseries-related file paths for each
    series/scan, from the dictionary of functional files described in the data
    configuration (sublist) YAML file.

    Scan input (from inputnode) is an iterable.
    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'scan', 'creds_path',
                                        'dl_dir'],
                                mandatory_inputs=True),
                        name='inputnode')

    outputnode = pe.Node(util.IdentityInterface(fields=['subject', 'rest',
                                                        'scan', 'scan_params',
                                                        'phase_diff',
                                                        'magnitude']),
                         name='outputspec')

    # have this here for now because of the big change in the data
    # configuration format
    check_scan = pe.Node(util.Function(input_names=['func_scan_dct',
                                                    'scan'],
                                       output_names=[],
                                       function=check_func_scan),
                         name='check_func_scan')

    check_scan.inputs.func_scan_dct = rest_dict
    wf.connect(inputnode, 'scan', check_scan, 'scan')

    # get the functional scan itself
    selectrest = pe.Node(util.Function(input_names=['scan',
                                                    'rest_dict',
                                                    'resource'],
                                       output_names=['file_path'],
                                       function=get_rest),
                         name='selectrest')
    selectrest.inputs.rest_dict = rest_dict
    selectrest.inputs.resource = "scan"
    wf.connect(inputnode, 'scan', selectrest, 'scan')

    # check to see if it's on an Amazon AWS S3 bucket, and download it, if it
    # is - otherwise, just return the local file path
    check_s3_node = pe.Node(util.Function(input_names=['file_path',
                                                       'creds_path',
                                                       'dl_dir',
                                                       'img_type'],
                                          output_names=['local_path'],
                                          function=check_for_s3),
                            name='check_for_s3')

    wf.connect(selectrest, 'file_path', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')
    wf.connect(inputnode, 'dl_dir', check_s3_node, 'dl_dir')
    check_s3_node.inputs.img_type = 'func'

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(check_s3_node, 'local_path', outputnode, 'rest')
    wf.connect(inputnode, 'scan', outputnode, 'scan')

    # scan parameters CSV
    select_scan_params = pe.Node(util.Function(input_names=['scan',
                                                            'rest_dict',
                                                            'resource'],
                                               output_names=['file_path'],
                                               function=get_rest),
                                 name='select_scan_params')
    select_scan_params.inputs.rest_dict = rest_dict
    select_scan_params.inputs.resource = "scan_parameters"
    wf.connect(inputnode, 'scan', select_scan_params, 'scan')

    # if the scan parameters file is on AWS S3, download it
    s3_scan_params = pe.Node(util.Function(input_names=['file_path',
                                                        'creds_path',
                                                        'dl_dir',
                                                        'img_type'],
                                           output_names=['local_path'],
                                           function=check_for_s3),
                             name='s3_scan_params')

    wf.connect(select_scan_params, 'file_path', s3_scan_params, 'file_path')
    wf.connect(inputnode, 'creds_path', s3_scan_params, 'creds_path')
    wf.connect(inputnode, 'dl_dir', s3_scan_params, 'dl_dir')
    wf.connect(s3_scan_params, 'local_path', outputnode, 'scan_params')

    # field map phase file, for field map distortion correction
    select_fmap_phase = pe.Node(util.Function(input_names=['scan',
                                                           'rest_dict',
                                                           'resource'],
                                              output_names=['file_path'],
                                              function=get_rest),
                                name='select_fmap_phase')
    select_fmap_phase.inputs.rest_dict = rest_dict
    select_fmap_phase.inputs.resource = "fmap_phase"
    wf.connect(inputnode, 'scan', select_fmap_phase, 'scan')

    s3_fmap_phase = pe.Node(util.Function(input_names=['file_path',
                                                       'creds_path',
                                                       'dl_dir',
                                                       'img_type'],
                                          output_names=['local_path'],
                                          function=check_for_s3),
                            name='s3_fmap_phase')
    s3_fmap_phase.inputs.img_type = "other"
    wf.connect(select_fmap_phase, 'file_path', s3_fmap_phase, 'file_path')
    wf.connect(inputnode, 'creds_path', s3_fmap_phase, 'creds_path')
    wf.connect(inputnode, 'dl_dir', s3_fmap_phase, 'dl_dir')
    wf.connect(s3_fmap_phase, 'local_path', outputnode, 'phase_diff')

    # field map magnitude file, for field map distortion correction
    select_fmap_mag = pe.Node(util.Function(input_names=['scan',
                                                         'rest_dict',
                                                         'resource'],
                                            output_names=['file_path'],
                                            function=get_rest),
                              name='select_fmap_mag')
    select_fmap_mag.inputs.rest_dict = rest_dict
    select_fmap_mag.inputs.resource = "fmap_mag"
    wf.connect(inputnode, 'scan', select_fmap_mag, 'scan')

    s3_fmap_mag = pe.Node(util.Function(input_names=['file_path',
                                                     'creds_path',
                                                     'dl_dir',
                                                     'img_type'],
                                        output_names=['local_path'],
                                        function=check_for_s3),
                          name='s3_fmap_mag')
    s3_fmap_mag.inputs.img_type = "other"
    wf.connect(select_fmap_mag, 'file_path', s3_fmap_mag, 'file_path')
    wf.connect(inputnode, 'creds_path', s3_fmap_mag, 'creds_path')
    wf.connect(inputnode, 'dl_dir', s3_fmap_mag, 'dl_dir')
    wf.connect(s3_fmap_mag, 'local_path', outputnode, 'magnitude')

    return wf


# Check if passed-in file is on S3
def check_for_s3(file_path, creds_path, dl_dir=None, img_type='anat'):

    # Import packages
    import os
    import nibabel as nib
    import botocore.exceptions

    from indi_aws import fetch_creds

    # Init variables
    s3_str = 's3://'
    if creds_path:
        if "None" in creds_path or "none" in creds_path or \
                "null" in creds_path:
            creds_path = None
    if dl_dir is None:
        dl_dir = os.getcwd()

    if file_path is None:
        # in case it's something like scan parameters or field map files, but
        # we don't have any
        local_path = file_path
        return local_path

    # TODO: remove this once scan parameter input as dictionary is phased out
    if isinstance(file_path, dict):
        # if this is a dictionary, just skip altogether
        local_path = file_path
        return local_path

    # Explicitly lower-case the "s3"
    if file_path.lower().startswith(s3_str):
        file_path_sp = file_path.split('/')
        file_path_sp[0] = file_path_sp[0].lower()
        file_path = '/'.join(file_path_sp)

    # Check for s3 string in filepaths
    if file_path.startswith(s3_str):
        # Get bucket name and bucket object
        bucket_name = file_path.replace(s3_str, '').split('/')[0]
        bucket = fetch_creds.return_bucket(creds_path, bucket_name)

        # Extract relative key path from bucket and local path
        s3_prefix = os.path.join(s3_str, bucket_name)
        s3_key = file_path.replace(s3_prefix, '').lstrip('/')
        local_path = os.path.join(dl_dir, s3_key)

        # Get local directory and create folders if they dont exist
        local_dir = os.path.dirname(local_path)
        if not os.path.exists(local_dir):
            os.makedirs(local_dir)

        # Download file
        try:
            print("Attempting to download from AWS S3: {0}".format(file_path))
            bucket.download_file(Key=s3_key, Filename=local_path)
        except botocore.exceptions.ClientError as exc:
            error_code = int(exc.response['Error']['Code'])
            if error_code == 403:
                err_msg = 'Access to bucket: "%s" is denied; using credentials '\
                          'in subject list: "%s"; cannot access the file "%s"'\
                          % (bucket_name, creds_path, file_path)
                raise Exception(err_msg)
            elif error_code == 404:
                err_msg = 'Bucket: "%s" does not exist; check spelling and try '\
                          'again' % bucket_name
                raise Exception(err_msg)
            else:
                err_msg = 'Unable to connect to bucket: "%s". Error message:\n%s'\
                          % (bucket_name, exc)
        except Exception as exc:
            err_msg = 'Unable to connect to bucket: "%s". Error message:\n%s'\
                      % (bucket_name, exc)
            raise Exception(err_msg)

    # Otherwise just return what was passed in
    else:
        local_path = file_path

    # Check image dimensionality
    if '.nii' in local_path:
        try:
            img_nii = nib.load(local_path)
        except Exception as e:
            # TODO: come up with a better option for handling rogue S3 files
            # TODO: that Nibabel chokes on
            print(str(e))
            return local_path

        if img_type == 'anat':
            if len(img_nii.shape) != 3:
                raise IOError('File: %s must be an anatomical image with 3 '\
                              'dimensions but %d dimensions found!'
                              % (local_path, len(img_nii.shape)))
        elif img_type == 'func':
            if len(img_nii.shape) != 4:
                raise IOError('File: %s must be a functional image with 4 '\
                              'dimensions but %d dimensions found!'
                              % (local_path, len(img_nii.shape)))
        elif img_type == "other":
            pass

    # Return the local path
    return local_path


def create_anat_datasource(wf_name='anat_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'anat', 'creds_path',
                                        'dl_dir'],
                                mandatory_inputs=True),
                        name='inputnode')

    check_s3_node = pe.Node(util.Function(input_names=['file_path',
                                                       'creds_path',
                                                       'dl_dir',
                                                       'img_type'],
                                          output_names=['local_path'],
                                          function=check_for_s3),
                            name='check_for_s3')

    wf.connect(inputnode, 'anat', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')
    wf.connect(inputnode, 'dl_dir', check_s3_node, 'dl_dir')
    check_s3_node.inputs.img_type = 'anat'

    outputnode = pe.Node(util.IdentityInterface(fields=['subject',
                                                        'anat']),
                         name='outputspec')

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(check_s3_node, 'local_path', outputnode, 'anat')

    # Return the workflow
    return wf


def create_roi_mask_dataflow(masks, wf_name='datasource_roi_mask'):

    import nipype.interfaces.io as nio
    import os

    wf = pe.Workflow(name=wf_name)  

    mask_dict = {}

    for mask_file in masks:

        mask_file = mask_file.rstrip('\r\n')

        if not os.path.exists(mask_file):
            err = '\n\n[!] CPAC says: One of your ROI/mask specification ' \
                  'files (under ROI TSE Options) does not have a correct ' \
                  'path or does not exist.\nTip: If all the paths are okay, '\
                  'then ensure there are no whitespaces or blank lines in ' \
                  'your ROI specification file.\n\n'
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
                  "files (under ROI TSE options) does not have '.nii' or " \
                  "'.nii.gz' as an extension.\n\nMask file: %s\n\n" \
                  % mask_file
            raise Exception(err)

        if not (base_name in mask_dict):
            mask_dict[base_name] = mask_file
        else:
            err = "\n\n[!] CPAC says: You have two or more ROI/mask files " \
            "with the same name - please make sure these files are named " \
            "differently.\n\nDuplicate name: %s\n\n" % mask_file
            raise Exception(err)

    inputnode = pe.Node(util.IdentityInterface(fields=['mask'],
                                               mandatory_inputs=True),
                        name='inputspec')

    inputnode.iterables = [('mask', mask_dict.keys())]

    selectmask = pe.Node(util.Function(input_names=['map', 'map_dct'],
                                       output_names=['out_file'],
                                       function=get_map),
                         name='select_mask')
    selectmask.inputs.map_dct = mask_dict

    outputnode = pe.Node(util.IdentityInterface(fields=['out_file']),
                         name='outputspec')

    wf.connect(inputnode, 'mask', selectmask, 'map')
    wf.connect(selectmask, 'out_file', outputnode, 'out_file')

    return wf


def create_spatial_map_dataflow(spatial_maps, wf_name='datasource_maps'):

    import os

    wf = pe.Workflow(name=wf_name)
    
    spatial_map_dict = {}
    
    for spatial_map_file in spatial_maps:

        spatial_map_file = spatial_map_file.rstrip('\r\n')

        if not os.path.exists(spatial_map_file):
            print "\n\nERROR: One of your spatial map files (under " \
                  "Spatial Regression options) does not have a correct " \
                  "path or does not exist. \nTip: If all the paths are " \
                  "okay, then ensure there are no whitespaces or blank " \
                  "lines in your spatial map specification file.\n\nError " \
                  "name: datasource_0001\n\n"
            raise Exception

        base_file = os.path.basename(spatial_map_file)
        base_name = ''
        try:
            if base_file.endswith('.nii'):
                base_name = os.path.splitext(base_file)[0]
            elif base_file.endswith('.nii.gz'):
                base_name = \
                    os.path.splitext(os.path.splitext(base_file)[0])[0]
            else:
                raise Exception("File extension not in .nii and .nii.gz "
                                "File: %s" % spatial_map_file)
        except Exception, e:
            print('error in spatial_map_dataflow: ', e)

        if not (base_name in spatial_map_dict):
            spatial_map_dict[base_name] = spatial_map_file
        else:
            raise ValueError('Files with same name not allowed %s %s'
                             % (spatial_map_file,
                                spatial_map_dict[base_name]))

    inputnode = pe.Node(util.IdentityInterface(fields=['spatial_map'],
                                               mandatory_inputs=True),
                        name='inputspec')

    inputnode.iterables = [('spatial_map', spatial_map_dict.keys())]

    select_spatial_map = pe.Node(util.Function(input_names=['map',
                                                            'map_dct'],
                                               output_names=['out_file'],
                                               function=get_map),
                                 name='select_spatial_map')
    select_spatial_map.inputs.map_dct = spatial_map_dict

    wf.connect(inputnode, 'spatial_map', select_spatial_map, 'map')

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


