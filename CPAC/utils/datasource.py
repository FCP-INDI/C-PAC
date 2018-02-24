import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def get_rest(scan, rest_dict):
    # return the time-series NIFTI file of the chosen series/scan
    return rest_dict[scan]["scan"]


def extract_scan_params_dct(scan_params_dct):
    return scan_params_dct


def create_func_datasource(rest_dict, wf_name='func_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'scan', 'creds_path'],
                                mandatory_inputs=True),
                        name='inputnode')

    outputnode = pe.Node(util.IdentityInterface(fields=['subject', 'rest',
                                                        'scan', 'scan_params',
                                                        'phase_diff',
                                                        'magnitude']),
                         name='outputspec')

    scan_names = rest_dict.keys()
    inputnode.iterables = [('scan', scan_names)]

    for scan in scan_names:

        scan_resources = rest_dict[scan]

        # actual 4D time series file
        if "scan" not in scan_resources.keys():
            err = "\n\n[!] The {0} scan is missing its actual time-series " \
                  "scan file, which should be a filepath labeled with the " \
                  "'scan' key.\n\n".format(scan)
            raise Exception(err)

        scan_file = scan_resources["scan"]

        if '.' in scan_file or '+' in scan_file or '*' in scan_file:
            raise Exception('\n\n[!] Scan names cannot contain any special '
                            'characters (., +, *, etc.). Please update this '
                            'and try again.\n\nScan: {0}'
                            '\n\n'.format(scan_file))

        # scan parameters CSV
        if "scan_parameters" in scan_resources.keys():
            if isinstance(scan_resources["scan_parameters"], str):
                if "s3://" in scan_resources["scan_parameters"]:
                    # if the scan parameters file is on AWS S3, download it
                    s3_scan_params = \
                        pe.Node(util.Function(input_names=['file_path',
                                                           'creds_path',
                                                           'img_type'],
                                              output_names=['local_path'],
                                              function=check_for_s3),
                                name='s3_scan_params')

                    s3_scan_params.inputs.file_path = \
                        scan_resources["scan_parameters"]

                    wf.connect(inputnode, 'creds_path',
                               s3_scan_params, 'creds_path')
                    wf.connect(s3_scan_params, 'local_path',
                               outputnode, 'scan_params')

            elif isinstance(scan_resources["scan_parameters"], dict):
                get_scan_params_dct = \
                        pe.Node(util.Function(input_names=['scan_params_dct'],
                                              output_names=['scan_params_dct'],
                                              function=extract_scan_params_dct),
                                name='get_scan_params_dct')
                get_scan_params_dct.inputs.scan_params_dct = \
                    scan_resources["scan_parameters"]
                wf.connect(get_scan_params_dct, 'scan_params_dct',
                           outputnode, 'scan_params')

        # field map files (if applicable)
        if "fmap_phase" in scan_resources.keys():

            fmap_phase = scan_resources["fmap_phase"]

            if "fmap_mag" not in scan_resources.keys():
                err = "\n\n[!] The field map phase difference file has " \
                      "been listed for scan {0}, but there is no field " \
                      "map magnitude file.\n\nPhase difference file " \
                      "listed: {1}\n\n".format(scan, fmap_phase)
                raise Exception(err)

            s3_fmap_phase = pe.Node(util.Function(input_names=['file_path',
                                                               'creds_path',
                                                               'img_type'],
                                                  output_names=['local_path'],
                                                  function=check_for_s3),
                                    name='s3_fmap_phase')
            s3_fmap_phase.inputs.file_path = fmap_phase
            s3_fmap_phase.inputs.img_type = "other"
            wf.connect(inputnode, 'creds_path', s3_fmap_phase, 'creds_path')
            wf.connect(s3_fmap_phase, 'local_path', outputnode, 'phase_diff')

        if "fmap_mag" in scan_resources.keys():

            fmap_mag = scan_resources["fmap_mag"]

            if "fmap_phase" not in scan_resources.keys():
                err = "\n\n[!] The field map magnitude file has been" \
                      "listed for scan {0}, but there is no field map phase" \
                      "difference file.\n\nPhase magnitude file " \
                      "listed: {1}\n\n".format(scan, fmap_mag)
                raise Exception(err)

            s3_fmap_mag = pe.Node(util.Function(input_names=['file_path',
                                                             'creds_path',
                                                             'img_type'],
                                                output_names=['local_path'],
                                                function=check_for_s3),
                                  name='s3_fmap_mag')
            s3_fmap_mag.inputs.file_path = fmap_mag
            s3_fmap_mag.inputs.img_type = "other"
            wf.connect(inputnode, 'creds_path', s3_fmap_mag, 'creds_path')
            wf.connect(s3_fmap_mag, 'local_path', outputnode, 'magnitude')

    selectrest = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['rest'],
                                       function=get_rest),
                         name='selectrest')
    selectrest.inputs.rest_dict = rest_dict

    check_s3_node = pe.Node(util.Function(input_names=['file_path',
                                                       'creds_path',
                                                       'img_type'],
                                          output_names=['local_path'],
                                          function=check_for_s3),
                            name='check_for_s3')

    wf.connect(selectrest, 'rest', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')
    check_s3_node.inputs.img_type = 'func'

    wf.connect(inputnode, 'scan', selectrest, 'scan')

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(check_s3_node, 'local_path', outputnode, 'rest')
    wf.connect(inputnode, 'scan', outputnode, 'scan')

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
    if dl_dir is None:
        dl_dir = os.getcwd()

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
        local_path = os.path.join(dl_dir, os.path.basename(s3_key))

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
                                fields=['subject', 'anat', 'creds_path'],
                                mandatory_inputs=True),
                        name='inputnode')

    check_s3_node = pe.Node(util.Function(input_names=['file_path',
                                                       'creds_path',
                                                       'img_type'],
                                          output_names=['local_path'],
                                          function=check_for_s3),
                            name='check_for_s3')

    wf.connect(inputnode, 'anat', check_s3_node, 'file_path')
    wf.connect(inputnode, 'creds_path', check_s3_node, 'creds_path')
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

    select_spatial_map = pe.Node(util.Function(input_names=['scan',
                                                            'rest_dict'],
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


