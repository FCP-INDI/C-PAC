# CPAC/utils/test_init.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module contains functions that assist in initializing CPAC
tests resources
'''

# Return tests data config file
def populate_template_config(config_type):
    '''
    Function to read in a template config file from the
    CPAC_RESOURCE_DIR and populate it with actual filepaths

    Parameters
    ----------
    config_type : string
        config file to populate; accepts 'data_config' and
        'pipeline_config'

    Returns
    -------
    config_test : string
        filepath to the newly written config file for testing
    '''

    # Import packages
    import os

    # Init variables
    resource_dir = return_resource_dir()
    templates_dir = return_resource_subfolder('templates')
    yamls = ['data_config', 'pipeline_config']
    texts = ['centrality_spec', 'map_spec', 'mask_spec',
             'roi_spec', 'seed_spec', 'spatial_maps_spec']

    # Check config type and build path
    if config_type in yamls:
        ext = '.yml'
        out_name = 'configs'
    elif config_type in texts:
        ext = '.txt'
        out_name = 'resources'
    else:
        # Check if it's supported, otherwise raise an Exception
        err_msg = 'config_type parameter: %s is unsupported' % config_type
        raise Exception(err_msg)

    # Get template and output paths
    template_path = os.path.join(templates_dir, config_type + ext)
    output_dir = return_resource_subfolder(out_name)
    output_path = os.path.join(output_dir, config_type + ext)

    # Open the files
    tmp_f = open(template_path, 'r')
    out_f = open(output_path, 'w')

    # Replace 'RESOURCE_DIR' string with actual directory
    for line in tmp_f:
        out_f.write(line.replace('RESOURCE_DIR', resource_dir))

    # Close file objects
    tmp_f.close()
    out_f.close()

    # Return filepath
    return output_path


# Populate all of the template paths
def populate_all_templates():
    '''
    Function to populate all of the template files

    Parameters
    ----------
    None

    Returns
    -------
    None
    '''

    # Import packages

    # Init variables
    outputs = []
    config_types = ['data_config', 'pipeline_config', 'centrality_spec',
                    'map_spec', 'mask_spec', 'roi_spec', 'seed_spec',
                    'spatial_maps_spec']

    # Populate all of the config templates with actual paths
    for config_type in config_types:
        output = populate_template_config(config_type)
        outputs.append(output)

    # Check that they all returned a value
    if len(outputs) == len(config_types):
        print 'Successfully populated and saved templates!'
    else:
        err_msg = 'Something went wrong during template population'
        raise Exception(err_msg)


# Get the AWS credentials
def return_aws_creds():
    '''
    Function to return the AWS credentials files located in the
    CPAC_RESOURCE_DIR

    Parameters
    ----------
    None

    Returns
    -------
    aws_creds : string
        filepath to the AWS credentials with access key id and secret
        access key
    db_creds: string
        filepath to the AWS database credentials with the username
        and password along with other details to connect to database
    '''

    # Import packages
    import os

    # Init variables
    creds_dir = return_resource_subfolder('creds')

    # Get credentials path
    aws_creds = os.path.join(creds_dir, 'aws_creds.csv')
    db_creds = os.path.join(creds_dir, 'db_creds.csv')

    # Return the aws creds and databse creds
    return aws_creds, db_creds


# Get the default test bucket name
def return_bucket_name():
    '''
    Function to return the default S3 bucket name used in test suite

    Parameters
    ----------
    None

    Returns
    -------
    bucket_name : string
        default S3 bucket name for testing
    '''

    # Set default bucket name
    bucket_name = 'fcp-indi'

    # Return bucket name
    return bucket_name


# Grab all nifti files within directory
def return_all_niis(base_dir):
    '''
    Function to walk through a base directory and all subsequent files
    and return the filepaths of all nifti files found

    Parameters
    ----------
    base_dir : string
        filepath to the base directory to search through

    Returns
    -------
    nii_list : list
        a list of filepath strings of the nifti files found in base_dir
    '''

    # Import packages
    import os

    # Init variables
    nii_list = []

    # Collect computed outputs
    for root, dirs, files in os.walk(base_dir):
        if files:
            nii_list.extend([os.path.join(root, file) for file in files \
                           if file.endswith('.nii.gz')])

    # Return the list of files
    return nii_list


# Look for CPAC_RESOURCE_DIR to be in environment
def return_resource_dir():
    '''
    Function to return the filepath of the CPAC_RESOURCE_DIR; note the
    CPAC_RESOURCE_DIR environment variable must be set

    Parameters
    ----------
    None

    Returns
    -------
    resource_dir : string
        the file path on disk where the cpac resources folder is
    '''
    # Import packages
    import os

    # Init variables
    resource_dir = os.getenv('CPAC_RESOURCE_DIR')

    # Check if set
    if not resource_dir:
        err_msg = 'CPAC_RESOURCE_DIR environment variable not set!\n' \
                  'Set this to the directory of the cpac_resources folder'
        raise Exception(err_msg)
    else:
        return resource_dir


# Return any subfolder of the resource directory
def return_resource_subfolder(subfolder):
    '''
    Funnction to return subfolders of the CPAC_RESOURCE_DIR

    Parameters
    ----------
    subfolder : string
        subfolder name to return path of

    Returns
    -------
    resource_subfolder : string
        filepath to the resource subfolder
    '''

    # Import packages
    import os

    # Init variables
    resource_dir = return_resource_dir()
    in_settings = ['configs', 'creds', 'resources',
                   'subject_lists', 'templates']

    # Check if its a sub-subfolder
    if subfolder in in_settings:
        resource_subfolder = os.path.join(resource_dir, 'settings', subfolder)
    else:
        resource_subfolder = os.path.join(resource_dir, subfolder)

    # Return subfolder
    return resource_subfolder


# Return test strategies obj file
def return_strats_obj():
    '''
    Function to return the file path of the strategies obj file from
    the CPAC_RESOURCE_DIR

    Parameters
    ----------
    None

    Returns
    -------
    strats_obj : string
        filepath to the strategies obj file
    '''

    # Import packages
    import os

    # Init variables
    settings_dir = return_resource_subfolder('resources')

    # Get strategies obj
    strats_obj = os.path.join(settings_dir, 'strategies_test.obj')

    # Return filepath
    return strats_obj


# Return tests subject list
def return_subject_list():
    '''
    Function to return the file path of the subject list file from
    the CPAC_RESOURCE_DIR

    Parameters
    ----------
    None

    Returns
    -------
    subject_list : string
        filepath to the subject list yaml file
    '''

    # Import packages
    import os

    # Init variables
    config_dir = return_resource_subfolder('subject_lists')

    # Get sublist
    subject_list = os.path.join(config_dir, 'CPAC_subject_list_test.yml')

    # Return filepath
    return subject_list


# Return the test subjects measure directories
def return_subj_measure_dirs(measure):
    '''
    Function to grab the base directories of the test subject's output
    files for a given measure or workflow

    Parameters
    ----------
    measure : string
        the measure or workflow or derivative of interest to parse for;
        this must be the folder name where all of the subject's test
        outputs are located (e.g. 'network_centrality')

    Returns
    -------
    subj_measure_dirs : list
        a list of strings of the base directories for each instance of
        the desired measure folder within the test subjects outputs
    '''

    # Import packages
    import glob
    import os

    # Init variables
    test_subj = return_test_subj()
    outputs_dir = return_resource_subfolder('output')

    # Root directories (cpac_resources/output/reg/subj_sess/scan/measure/..)
    subj_measure_dirs = \
        glob.glob(os.path.join(outputs_dir, '*', '%s*' % test_subj,
                               '*', measure))

    # Check to see if the directories exist
    if len(subj_measure_dirs) == 0:
        err_msg = 'Unable to find any subject directories for the %s measure.' \
                  % measure
        raise Exception(err_msg)

    # Return base directories for test measures outputs
    return subj_measure_dirs


# Get subject for individual tests
def return_test_subj():
    '''
    Function to return the subject id; note the
    CPAC_RESOURCE_DIR environment variable must be set

    Parameters
    ----------
    None

    Returns
    -------
    resource_dir : string
        the file path on disk where the cpac resources folder is
    '''

    # Import packages
    import os

    # Init variables
    test_subj = os.getenv('CPAC_TEST_SUBJ')

    # Get cpac resource directory and get a list of subject folders
    input_dir = return_resource_subfolder('input')
    site_dir = os.path.join(input_dir, 'site_1')

    # Get list of subject directories
    subs = os.listdir(site_dir)

    # Check if set and exists
    if not test_subj:
        err_msg = 'CPAC_TEST_SUBJ environment variable not set!\n' \
                  'Set this to the subject id of the desired subject to test' \
                  'from the cpac_resources folder.'
        raise Exception(err_msg)
    elif test_subj not in subs:
        err_msg = 'Test subject %s is not in the cpac_resources subject ' \
                  'directory %s. Please specify different CPAC_TEST_SUBJ.' \
                  %(test_subj, site_dir)
    else:
        return test_subj


# Smooth nifti file
def smooth_nii_file(self, nii_file, fwhm, mask_file=None):
    '''
    Function to Gaussian smooth nifti files and optionally using a mask
    on the smoothed data

    Parameters
    ----------
    nii_file : string
        filepath to the nifti file to smooth
    fwhm : float
        FWHM for Gaussian smoothing kernel, in mm
    mask_file : string (optional); default=None
        filepath to the mask file to use

    Returns
    -------
    smooth_arr : numpy.ndarray
        smoothed nifti image as a numpy array
    '''

    # Import packages
    import nibabel as nib
    import numpy as np
    import scipy.ndimage

    # Init variables
    raw_nii = nib.load(nii_file)
    raw_arr = raw_nii.get_data()

    # Check parameters
    if mask_file:
        mask_arr = nib.load(mask_file).get_data()
        # Check the mask shape matches the raw nifti
        if mask_arr.shape != raw_arr.shape:
            err_msg = 'Mask file has different dimensions than nifti.\n' \
                      'Check the paths are correct and try again.'
            raise Exception(err_msg)

    # Calculate sigma for smoothing
    mm_res = np.abs(raw_nii.affine[0][0])
    sigma = fwhm/2.3548/mm_res

    # Smooth input
    smooth_arr = scipy.ndimage.gaussian_filter(raw_arr, sigma, order=0)

    # And mask if using one (this writes it to a 1d array)
    if mask_arr:
        smooth_out = smooth_arr[mask_arr.astype('bool')]
        smooth_arr = np.zeros(mask_arr.shape, dtype=float)

        # Get mask coordinates and populate smoothed image
        coords = np.argwhere(mask_arr)
        for idx, xyz in enumerate(coords):
            x, y, z = xyz
            smooth_arr[x, y, z] = smooth_out[idx]

    # Return the smoothed array
    return smooth_arr

# Setup log file
def setup_test_logger(logger_name, log_file, level, to_screen=False):
    '''
    Function to initialize and configure a logger that can write to file
    and (optionally) the screen.

    Parameters
    ----------
    logger_name : string
        name of the logger
    log_file : string
        file path to the log file on disk
    level : integer
        indicates the level at which the logger should log; this is
        controlled by integers that come with the python logging
        package. (e.g. logging.INFO=20, logging.DEBUG=10)
    to_screen : boolean (optional)
        flag to indicate whether to enable logging to the screen

    Returns
    -------
    logger : logging.Logger object
        Python logging.Logger object which is capable of logging run-
        time information about the program to file and/or screen
    '''

    # Import packages
    import logging

    # Init logger, formatter, filehandler, streamhandler
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    formatter = logging.Formatter('%(asctime)s : %(message)s')

    # Write logs to file
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Write to screen, if desired
    if to_screen:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    # Return the logger
    return logger