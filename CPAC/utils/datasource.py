import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.utils import function


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
        scan_resources.keys()
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
    check_scan = pe.Node(function.Function(input_names=['func_scan_dct',
                                                        'scan'],
                                           output_names=[],
                                           function=check_func_scan,
                                           as_module=True),
                         name='check_func_scan')

    check_scan.inputs.func_scan_dct = rest_dict
    wf.connect(inputnode, 'scan', check_scan, 'scan')

    # get the functional scan itself
    selectrest = pe.Node(function.Function(input_names=['scan',
                                                        'rest_dict',
                                                        'resource'],
                                           output_names=['file_path'],
                                           function=get_rest,
                                           as_module=True),
                         name='selectrest')
    selectrest.inputs.rest_dict = rest_dict
    selectrest.inputs.resource = "scan"
    wf.connect(inputnode, 'scan', selectrest, 'scan')

    # check to see if it's on an Amazon AWS S3 bucket, and download it, if it
    # is - otherwise, just return the local file path
    check_s3_node = pe.Node(function.Function(input_names=['file_path',
                                                           'creds_path',
                                                           'dl_dir',
                                                           'img_type'],
                                              output_names=['local_path'],
                                              function=check_for_s3,
                                              as_module=True),
                            name='check_for_s3')

    wf.connect(selectrest, 'file_path', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')
    wf.connect(inputnode, 'dl_dir', check_s3_node, 'dl_dir')
    check_s3_node.inputs.img_type = 'func'

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(check_s3_node, 'local_path', outputnode, 'rest')
    wf.connect(inputnode, 'scan', outputnode, 'scan')

    # scan parameters CSV
    select_scan_params = pe.Node(function.Function(input_names=['scan',
                                                                'rest_dict',
                                                                'resource'],
                                                   output_names=['file_path'],
                                                   function=get_rest,
                                                   as_module=True),
                                 name='select_scan_params')
    select_scan_params.inputs.rest_dict = rest_dict
    select_scan_params.inputs.resource = "scan_parameters"
    wf.connect(inputnode, 'scan', select_scan_params, 'scan')

    # if the scan parameters file is on AWS S3, download it
    s3_scan_params = pe.Node(function.Function(input_names=['file_path',
                                                            'creds_path',
                                                            'dl_dir',
                                                            'img_type'],
                                               output_names=['local_path'],
                                               function=check_for_s3,
                                               as_module=True),
                             name='s3_scan_params')

    wf.connect(select_scan_params, 'file_path', s3_scan_params, 'file_path')
    wf.connect(inputnode, 'creds_path', s3_scan_params, 'creds_path')
    wf.connect(inputnode, 'dl_dir', s3_scan_params, 'dl_dir')
    wf.connect(s3_scan_params, 'local_path', outputnode, 'scan_params')

    # field map phase file, for field map distortion correction
    select_fmap_phase = pe.Node(function.Function(input_names=['scan',
                                                               'rest_dict',
                                                               'resource'],
                                                  output_names=['file_path'],
                                                  function=get_rest,
                                                  as_module=True),
                                name='select_fmap_phase')
    select_fmap_phase.inputs.rest_dict = rest_dict
    select_fmap_phase.inputs.resource = "fmap_phase"
    wf.connect(inputnode, 'scan', select_fmap_phase, 'scan')

    s3_fmap_phase = pe.Node(function.Function(input_names=['file_path',
                                                           'creds_path',
                                                           'dl_dir',
                                                           'img_type'],
                                              output_names=['local_path'],
                                              function=check_for_s3,
                                              as_module=True),
                            name='s3_fmap_phase')
    s3_fmap_phase.inputs.img_type = "other"
    wf.connect(select_fmap_phase, 'file_path', s3_fmap_phase, 'file_path')
    wf.connect(inputnode, 'creds_path', s3_fmap_phase, 'creds_path')
    wf.connect(inputnode, 'dl_dir', s3_fmap_phase, 'dl_dir')
    wf.connect(s3_fmap_phase, 'local_path', outputnode, 'phase_diff')

    # field map magnitude file, for field map distortion correction
    select_fmap_mag = pe.Node(function.Function(input_names=['scan',
                                                             'rest_dict',
                                                             'resource'],
                                                output_names=['file_path'],
                                                function=get_rest,
                                                as_module=True),
                              name='select_fmap_mag')
    select_fmap_mag.inputs.rest_dict = rest_dict
    select_fmap_mag.inputs.resource = "fmap_mag"
    wf.connect(inputnode, 'scan', select_fmap_mag, 'scan')

    s3_fmap_mag = pe.Node(function.Function(input_names=['file_path',
                                                         'creds_path',
                                                         'dl_dir',
                                                         'img_type'],
                                            output_names=['local_path'],
                                            function=check_for_s3,
                                            as_module=True),
                          name='s3_fmap_mag')
    s3_fmap_mag.inputs.img_type = "other"
    wf.connect(select_fmap_mag, 'file_path', s3_fmap_mag, 'file_path')
    wf.connect(inputnode, 'creds_path', s3_fmap_mag, 'creds_path')
    wf.connect(inputnode, 'dl_dir', s3_fmap_mag, 'dl_dir')
    wf.connect(s3_fmap_mag, 'local_path', outputnode, 'magnitude')

    return wf


def create_check_for_s3_node(name, file_path, img_type='other', creds_path=None, dl_dir=None):

    check_s3_node = pe.Node(function.Function(input_names=['file_path',
                                                           'creds_path',
                                                           'dl_dir',
                                                           'img_type'],
                                              output_names=['local_path'],
                                              function=check_for_s3,
                                              as_module=True),
                            name='check_for_s3_%s' % name)

    check_s3_node.inputs.set(
        file_path=file_path,
        creds_path=creds_path,
        dl_dir=dl_dir,
        img_type=img_type
    )

    return check_s3_node


# Check if passed-in file is on S3
def check_for_s3(file_path, creds_path=None, dl_dir=None, img_type='other'):

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
        return None

    # TODO: remove this once scan parameter input as dictionary is phased out
    if isinstance(file_path, dict):
        # if this is a dictionary, just skip altogether
        local_path = file_path
        return local_path

    # Explicitly lower-case the "s3"
    if file_path.lower().startswith(s3_str):
        
        file_path = s3_str + file_path[len(s3_str):]

        # Get bucket name and bucket object
        bucket_name = file_path[len(s3_str):].split('/')[0]
        bucket = fetch_creds.return_bucket(creds_path, bucket_name)

        # Extract relative key path from bucket and local path
        s3_prefix = s3_str + bucket_name
        s3_key = file_path[len(s3_prefix) + 1:]
        local_path = os.path.join(dl_dir, bucket_name, s3_key)

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

            err_msg = str(exc)
            if error_code == 403:
                err_msg = 'Access to bucket: "%s" is denied; using credentials '\
                          'in subject list: "%s"; cannot access the file "%s"'\
                          % (bucket_name, creds_path, file_path)
            elif error_code == 404:
                err_msg = 'File: {0} does not exist; check spelling and try '\
                          'again'.format(os.path.join(bucket_name, s3_key))
            else:
                err_msg = 'Unable to connect to bucket: "%s". Error message:\n%s'\
                          % (bucket_name, exc)
            
            raise Exception(err_msg)

        except Exception as exc:
            err_msg = 'Unable to connect to bucket: "%s". Error message:\n%s'\
                      % (bucket_name, exc)
            raise Exception(err_msg)

    # Otherwise just return what was passed in
    else:
        local_path = file_path

    # Check if it exists or it is sucessfuly downloaded
    if not os.path.exists(local_path):
        raise IOError('File %s does not exists!' % (local_path))

    # Check image dimensionality
    if local_path.endswith('.nii') or local_path.endswith('.nii.gz'):
        img_nii = nib.load(local_path)

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

    check_s3_node = pe.Node(function.Function(input_names=['file_path',
                                                           'creds_path',
                                                           'dl_dir',
                                                           'img_type'],
                                              output_names=['local_path'],
                                              function=check_for_s3,
                                              as_module=True),
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

    import os

    mask_dict = {}

    for mask_file in masks:

        mask_file = mask_file.rstrip('\r\n')

        if mask_file.strip() == '' or mask_file.startswith('#'):
            continue

        base_file = os.path.basename(mask_file)

        try:
            valid_extensions = ['.nii', '.nii.gz']

            base_name = [
                base_file[:-len(ext)]
                for ext in valid_extensions
                if base_file.endswith(ext)
            ][0]

            if base_name in mask_dict:
                raise ValueError(
                    'Files with same name not allowed: %s %s' % (
                        mask_file,
                        mask_dict[base_name]
                    )
                )

            mask_dict[base_name] = mask_file

        except IndexError as e:
            raise Exception('Error in spatial_map_dataflow: '
                            'File extension not in .nii and .nii.gz')

        except Exception as e:
            raise e


    wf = pe.Workflow(name=wf_name)  

    inputnode = pe.Node(util.IdentityInterface(fields=['mask',
                                                       'mask_file',
                                                       'creds_path',
                                                       'dl_dir'],
                                               mandatory_inputs=True),
                        name='inputspec')

    mask_keys, mask_values = \
        zip(*mask_dict.items())

    inputnode.synchronize = True
    inputnode.iterables = [
        ('mask', mask_keys),
        ('mask_file', mask_values),
    ]

    check_s3_node = pe.Node(function.Function(input_names=['file_path',
                                                           'creds_path',
                                                           'dl_dir',
                                                           'img_type'],
                                              output_names=['local_path'],
                                              function=check_for_s3,
                                              as_module=True),
                            name='check_for_s3')

    wf.connect(inputnode, 'mask_file', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')
    wf.connect(inputnode, 'dl_dir', check_s3_node, 'dl_dir')
    check_s3_node.inputs.img_type = 'mask'

    outputnode = pe.Node(util.IdentityInterface(fields=['out_file']),
                         name='outputspec')

    wf.connect(check_s3_node, 'local_path', outputnode, 'out_file')

    return wf


def create_spatial_map_dataflow(spatial_maps, wf_name='datasource_maps'):

    import os

    wf = pe.Workflow(name=wf_name)
    
    spatial_map_dict = {}
    
    for spatial_map_file in spatial_maps:

        spatial_map_file = spatial_map_file.rstrip('\r\n')
        base_file = os.path.basename(spatial_map_file)

        try:
            valid_extensions = ['.nii', '.nii.gz']

            base_name = [
                base_file[:-len(ext)]
                for ext in valid_extensions
                if base_file.endswith(ext)
            ][0]

            if base_name in spatial_map_dict:
                raise ValueError(
                    'Files with same name not allowed: %s %s' % (
                        spatial_map_file,
                        spatial_map_dict[base_name]
                    )
                )
        
            spatial_map_dict[base_name] = spatial_map_file
        
        except IndexError as e:
            raise Exception('Error in spatial_map_dataflow: '
                            'File extension not in .nii and .nii.gz')

    inputnode = pe.Node(util.IdentityInterface(fields=['spatial_map',
                                                       'spatial_map_file',
                                                       'creds_path',
                                                       'dl_dir'],
                                               mandatory_inputs=True),
                        name='inputspec')

    spatial_map_keys, spatial_map_values = \
        zip(*spatial_map_dict.items())

    inputnode.synchronize = True
    inputnode.iterables = [
        ('spatial_map', spatial_map_keys),
        ('spatial_map_file', spatial_map_values),
    ]

    check_s3_node = pe.Node(function.Function(input_names=['file_path',
                                                           'creds_path',
                                                           'dl_dir',
                                                           'img_type'],
                                              output_names=['local_path'],
                                              function=check_for_s3,
                                              as_module=True),
                            name='check_for_s3')

    wf.connect(inputnode, 'spatial_map_file', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')
    wf.connect(inputnode, 'dl_dir', check_s3_node, 'dl_dir')
    check_s3_node.inputs.img_type = 'mask'

    select_spatial_map = pe.Node(util.IdentityInterface(fields=['out_file'],
                                                        mandatory_inputs=True),
                                 name='select_spatial_map')

    wf.connect(check_s3_node, 'local_path', select_spatial_map, 'out_file')

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

    selectmodel = pe.Node(function.Function(input_names=['model',
                                                         'ftest',
                                                         'model_name'],
                                            output_names=['fts_file',
                                                          'con_file',
                                                          'grp_file',
                                                          'mat_file'],
                                            function=select_model_files,
                                            as_module=True),
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
